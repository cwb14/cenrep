"""S2: satellite family library from TRF seeds (rotation-invariant, multimer-aware)."""
import random
from pathlib import Path

import numpy as np
import pandas as pd
import pyhmmer
import pysam
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextMSA, TextSequence

from .annotate import run_significant
from .params import (BLAST_BATCH, FAMILY_IDENT, FULL_MONOMER_COV, FULL_PERIOD_TOL, LEN_TOL, MSA_MAX_GAP,
                     NHMMER_REPORT_E, OVERHANG_MIN_BP, RDNA_FLANK_EVALUE, RDNA_TAG_EVALUE, RELATED_LEN_TOL)
from .util import (CenrepError, canonical_rotation, log, merge_intervals, parse_fasta_text, revcomp, run, union_bp,
                   write_fasta)

ABC = Alphabet.dna()
RRNA_CM = Path(__file__).resolve().parent / "data" / "rfam_rrna.cm"   # Rfam 5S/5.8S/SSU/LSU (package data, data/README.md)
BLAST_FMT = "6 qseqid sseqid pident length qstart qend sstart send qlen slen bitscore"
FAM_COLS = ["fam_id", "rep_seed", "base_seed", "consensus_len", "n_seeds", "max_genome_span",
            "genome_spans", "hmm_mode", "n_hmm_seqs", "superfamily", "rdna", "chain_group"]


def select_seeds(trf: dict, p) -> pd.DataFrame:
    parts = []
    for g, df in trf.items():
        d = df[(df.period >= p.min_period) & (df.period <= p.max_period)
               & ((df.end - df.start) >= p.min_seed_span)].copy()
        d["genome"] = g
        parts.append(d)
    s = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if s.empty:
        return s
    s["span"] = s.end - s.start
    s = s.sort_values(["span", "genome", "chrom", "start"],
                      ascending=[False, True, True, True]).reset_index(drop=True)
    s["seed_id"] = [f"s{i}" for i in range(len(s))]
    return s


def base_unit(mem: pd.DataFrame, fastas: dict):
    """The sequence a family model is founded on: one whole monomer. Multimers have periods near k x the monomer
    (k >= 2), so the monomer class is every seed with a period below 1.5 x the shortest; within it the largest-span
    seed whose consensus spans its whole period is used. TRF drops alignment columns at indels, so a consensus can
    be shorter than its period (A. arenosa pAa: period 179, consensus 170) and a model founded on it misses part of
    every monomer; if no monomer-class seed has a whole consensus, one period of the largest array is taken from
    the genome. mem is sorted by span, largest first. Returns (seed row, sequence)."""
    mono = mem[mem.period < 1.5 * mem.period.min()]
    full = mono[(mono.consensus.str.len() - mono.period).abs() <= np.maximum(1, FULL_PERIOD_TOL * mono.period)]
    if len(full):
        b = full.iloc[0]
        return b, b.consensus.upper()
    b = mono.iloc[0]
    with pysam.FastaFile(str(fastas[b.genome])) as fa:
        return b, fa.fetch(b.chrom, int(b.start), int(b.start + b.period)).upper()


def seed_match(slen: int, llen: int, hsps, t_same: float, t_multi: float) -> bool:
    """Shorter seed (length slen, BLAST query) vs the longer seed's dimer (2*llen, subject).
    hsps: [(pident, qstart, qend, sstart, send), ...], BLAST 1-based coordinates."""
    ratio = llen / slen
    k = max(1, round(ratio))
    if abs(ratio - k) > LEN_TOL * k:
        return False
    t = 100 * (t_same if k == 1 else t_multi)
    good = [h for h in hsps if h[0] >= t]
    if not any(h[2] - h[1] + 1 >= 0.8 * slen for h in good):
        return False
    cov = np.zeros(llen, bool)
    for h in good:
        a, b = min(h[3], h[4]) - 1, max(h[3], h[4])
        cov[np.arange(a, b) % llen] = True          # circular coverage of the longer monomer
    return cov.mean() >= 0.8


def blast_pairs(queries: dict, subjects: dict, workdir: Path, threads: int) -> dict:
    """BLAST query monomers against subject dimers. Returns {(qid, sid): [hsp, ...]}."""
    workdir.mkdir(parents=True, exist_ok=True)
    q, s, db, out = workdir / "q.fa", workdir / "s.fa", workdir / "sdb", workdir / "hits.tsv"
    write_fasta(q, queries.items())
    write_fasta(s, ((k, v + v) for k, v in subjects.items()))
    run(["makeblastdb", "-in", s, "-dbtype", "nucl", "-out", db])
    # thread mode only (hits are identical): split by query when queries outnumber subjects (many seeds vs a few
    # representatives); BLAST's default otherwise -- splitting 68 queries over 25 threads kept 2 busy (Avena)
    mt = ["-mt_mode", 1] if len(queries) > len(subjects) else []
    run(["blastn", "-task", "blastn", "-query", q, "-db", db, "-word_size", 11, "-evalue", "1e-5",
         "-dust", "no", "-soft_masking", "false", "-max_target_seqs", 100000,
         "-num_threads", threads, *mt, "-outfmt", BLAST_FMT, "-out", out])
    pairs = {}
    with open(out) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            pairs.setdefault((x[0], x[1]), []).append(
                (float(x[2]), int(x[4]), int(x[5]), int(x[6]), int(x[7])))
    return pairs


def _related(a, b, seqs, pairs, t_same, t_multi) -> bool:
    """Candidate a vs representative b. The lenient multimer identity (t_multi) applies only when a is the
    longer (a multimer joining a monomer representative). A monomer joining a multimer representative must
    meet t_same: TRF often calls the largest array at a multimer period, and at t_multi a related but distinct
    subfamily (CentBd/CentBs, CEN180/pAa) would join it depending on seed order (review C1)."""
    la, lb = len(seqs[a]), len(seqs[b])
    if la <= lb:
        return seed_match(la, lb, pairs.get((a, b), []), t_same, t_same)
    return seed_match(lb, la, pairs.get((b, a), []), t_same, t_multi)


def cluster_seeds(seqs: dict, order: list, workdir: Path, threads: int, t_same: float, t_multi: float) -> dict:
    """Greedy centroid clustering in priority order (no chaining). Returns {id: representative id}.
    Works in batches of BLAST_BATCH candidate representatives so BLAST output stays bounded."""
    assign, rest, rnd = {}, list(order), 0
    while rest:
        rnd += 1
        cand, others = rest[:BLAST_BATCH], rest[BLAST_BATCH:]
        sub = {i: seqs[i] for i in cand}
        pairs = blast_pairs(sub, sub, workdir / f"r{rnd}c", threads)
        reps = []
        for c in cand:
            hit = next((r for r in reps if _related(c, r, seqs, pairs, t_same, t_multi)), None)
            assign[c] = c if hit is None else hit
            if hit is None:
                reps.append(c)
        if others:
            o, r = {i: seqs[i] for i in others}, {i: seqs[i] for i in reps}
            pairs = blast_pairs(o, r, workdir / f"r{rnd}o1", threads)
            for k, v in blast_pairs(r, o, workdir / f"r{rnd}o2", threads).items():
                pairs.setdefault(k, []).extend(v)
            for x in others:
                hit = next((rp for rp in reps if _related(x, rp, seqs, pairs, t_same, t_multi)), None)
                if hit is not None:
                    assign[x] = hit
        rest = [x for x in others if x not in assign]
        log.debug(f"families: round {rnd}: {len(reps)} new representatives, {len(rest)} seeds left")
    return assign


def absorb_multimers(cons: dict, workdir: Path, threads: int, ident: float) -> dict:
    """Map each family whose consensus is a k-fold multimer (k >= 2) of a shorter family's consensus
    (>= ident per HSP, >= 80% circular tiling; seed_match rules) to that monomer family, following chains
    to the root. TRF reports HOR-rich arrays at the HOR period; if some monomers of the HOR copy are too
    divergent for the 80% multimer rule of the clustering, the multimer would otherwise become its own
    family and its longer HMM hits would displace monomer hits at annotation. Returns {child: root}."""
    if len(cons) < 2:
        return {}
    pairs = blast_pairs(cons, cons, workdir, threads)
    parent = {}
    for f, sf in cons.items():
        best = None
        for g, sg in cons.items():
            if g == f or len(sg) >= len(sf):
                continue
            ratio = len(sf) / len(sg)
            k = round(ratio)
            if k < 2 or abs(ratio - k) > LEN_TOL * k:
                continue
            hsps = pairs.get((g, f), [])
            if seed_match(len(sg), len(sf), hsps, ident, ident):
                score = float(np.mean([h[0] for h in hsps if h[0] >= 100 * ident]))
                if best is None or score > best[0]:
                    best = (score, g)
        if best:
            parent[f] = best[1]
    root = {}
    for f in parent:
        r = f
        while r in parent:
            r = parent[r]
        root[f] = r
    return root


def rdna_tags(cons: dict, workdir: Path, threads: int, cm=RRNA_CM, evalue: float = RDNA_TAG_EVALUE) -> dict:
    """Tag families whose consensus (searched as a dimer, so a gene spanning the rotation point is found)
    contains an rRNA gene (Rfam 5S/5.8S/SSU/LSU, cmsearch). rDNA arrays (e.g. pericentromeric 5S rDNA)
    look like satellites but are never centromeric (review C2). Returns {fam: rRNA model name or '.'}."""
    tags = {f: "." for f in cons}
    if not cons:
        return tags
    if not Path(cm).exists():
        raise CenrepError(f"rRNA models missing: {cm} (part of the cenrep package; reinstall it)")
    workdir.mkdir(parents=True, exist_ok=True)
    q, tbl = workdir / "cons_dimers.fa", workdir / "rrna.tbl"
    write_fasta(q, ((f, s + s) for f, s in cons.items()))
    run(["cmsearch", "--cpu", threads, "-E", evalue, "--noali", "--tblout", tbl, cm, q])
    best = {}
    for line in open(tbl):
        if line.startswith("#"):
            continue
        x = line.split()
        if x[16] != "!":
            continue
        if x[0] not in best or float(x[15]) < best[x[0]][0]:
            best[x[0]] = (float(x[15]), x[2])
    for f, (_, name) in best.items():
        tags[f] = name
    return tags


RRNA_45S = ("SSU_rRNA_eukarya", "5_8S_rRNA", "LSU_rRNA_eukarya")


def rdna_neighbourhood(fasta, chroms: dict, arrays: pd.DataFrame, fams: list, workdir: Path, threads: int,
                       flank: int = 15_000, n_arrays: int = 10, cm=RRNA_CM, evalue: float = RDNA_FLANK_EVALUE) -> dict:
    """45S rDNA intergenic-spacer subrepeats form satellite-like arrays whose consensus holds no rRNA gene, so
    rdna_tags misses them (P. infirma chr1A NOR: Poa_sat18_125 was called centromeric). Search the flanks
    (+/- flank) of each family's n_arrays largest arrays for 45S rRNA genes (Rfam SSU/5.8S/LSU; cmsearch
    --hmmonly, ~70x faster than the CM and as sensitive for these conserved genes).
    Returns {fam: (arrays tested, arrays with a 45S rRNA gene in a flank)}."""
    out = {f: (0, 0) for f in fams}
    if not fams:
        return out
    if not Path(cm).exists():
        raise CenrepError(f"rRNA models missing: {cm} (part of the cenrep package; reinstall it)")
    workdir.mkdir(parents=True, exist_ok=True)
    q, tbl = workdir / "flanks.fa", workdir / "flanks.tbl"
    tested = {}
    with pysam.FastaFile(str(fasta)) as fa, open(q, "w") as fh:
        for f in fams:
            a = arrays[arrays.family == f]
            a = a.assign(_len=a.end - a.start).sort_values("_len", ascending=False).head(n_arrays)
            tested[f] = len(a)
            for i, r in enumerate(a.itertuples(index=False)):
                L = chroms[r.chrom]
                for side, (s, e) in (("L", (max(0, r.start - flank), r.start)), ("R", (r.end, min(L, r.end + flank)))):
                    if e - s >= 100:
                        fh.write(f">{f}|{i}|{side}\n{fa.fetch(r.chrom, s, e)}\n")
    if not any(tested.values()):
        return out
    run(["cmsearch", "--hmmonly", "--cpu", threads, "-E", evalue, "--noali", "--tblout", tbl, cm, q])
    hit = set()
    for line in open(tbl):
        if line.startswith("#"):
            continue
        x = line.split()
        if x[2] in RRNA_45S:
            f, i, _ = x[0].rsplit("|", 2)
            hit.add((f, i))
    return {f: (tested[f], sum(1 for g, _ in hit if g == f)) for f in fams}


def trim_overhang(seq: str, workdir: Path, threads: int, min_len: int = OVERHANG_MIN_BP,
                  ident: float = FAMILY_IDENT) -> str:
    """Trim a terminally redundant consensus to its true period. A consensus that is one unit plus the start of the
    next unit (Msookensis: Cent728 + 68 bp = 796 bp) has its prefix recurring at position p and running to the end;
    its HMM then outscores the true monomer model on every array. Internal motif repeats (not reaching the end)
    and short chance matches (< min_len) are left alone."""
    seq = seq.upper()
    n = len(seq)
    if n < 2 * min_len:
        return seq
    pairs = blast_pairs({"q": seq}, {"s": seq}, workdir, threads)          # subject searched as seq+seq
    best = None
    for pid, qs, qe, ss, se in pairs.get(("q", "s"), []):
        p = ss - qs
        if (se > ss and qs <= 5 and se >= n - 5 and pid >= 100 * ident and qe - qs + 1 >= min_len
                and n / 2 <= p <= n - min_len and (best is None or p < best)):
            best = p
    return seq[:best] if best else seq


def frame_rotation(seq: str, refs: dict, workdir: Path, threads: int, ident: float, _flipped: bool = False) -> str:
    """Rotate (and orient) a monomer consensus into the frame of the most similar related reference family
    (length within 15%, HSP >= ident over >= 80% of seq); canonical rotation if none is related. Independent
    lexicographic rotations put near-identical subfamilies out of phase (110/178 bp), so their hits
    half-overlapped and a tenth or more of the array bp stayed unannotated."""
    seq = seq.upper()
    cand = {k: v for k, v in refs.items() if abs(len(v) / len(seq) - 1) <= RELATED_LEN_TOL}
    if cand:
        pairs = blast_pairs({"q": seq}, cand, workdir, threads)
        best = None
        for k in cand:
            for pid, qs, qe, ss, se in pairs.get(("q", k), []):
                if pid >= 100 * ident and qe - qs + 1 >= 0.8 * len(seq) and (best is None or pid > best[0]):
                    best = (pid, qs, ss, se, len(cand[k]))
        if best:
            _, qs, ss, se, L = best
            if se < ss:                                   # opposite orientation: flip once and realign
                if not _flipped:
                    return frame_rotation(revcomp(seq), refs, workdir, threads, ident, True)
            else:
                r = ((qs - 1) - (ss - 1) % L) % len(seq)
                return seq[r:] + seq[:r]
    return canonical_rotation(seq)


def dedup_families(cons: dict, order: list, workdir: Path, threads: int, ident: float) -> dict:
    """Map near-identical families (length within 5%, HSP >= ident over >= 90% of the shorter, any rotation)
    to the earliest such family in order. TRF-seed clustering at 90% can yield families whose refined consensus
    sequences are 98-100% identical; their HMMs would then split monomers between them arbitrarily."""
    if len(cons) < 2:
        return {}
    pairs = blast_pairs(cons, cons, workdir, threads)
    merged = {}
    for i, f in enumerate(order):
        if f in merged:
            continue
        for g in order[i + 1:]:
            if g in merged or abs(len(cons[g]) / len(cons[f]) - 1) > LEN_TOL:
                continue
            short = min(len(cons[f]), len(cons[g]))
            if any(h[0] >= 100 * ident and h[2] - h[1] + 1 >= 0.9 * short for h in pairs.get((g, f), [])):
                merged[g] = f
    return merged


def chain_groups(cons: dict, order: list, workdir: Path, threads: int, ident: float) -> dict:
    """Connected components of close subfamilies (HSP >= ident over >= 80% of the shorter, length within 15%);
    arrays are chained per group so interleaved variant monomers (CentO-165 among CentO-155) are not dropped as
    per-family singletons. Returns {fam: group representative = earliest member in order}."""
    parent = {f: f for f in cons}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    if len(cons) > 1:
        pairs = blast_pairs(cons, cons, workdir, threads)
        for (a, b), hs in pairs.items():
            la, lb = len(cons[a]), len(cons[b])
            if a == b or abs(la / lb - 1) > RELATED_LEN_TOL:
                continue
            if any(h[0] >= 100 * ident and h[2] - h[1] + 1 >= 0.8 * min(la, lb) for h in hs):
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[max(ra, rb, key=order.index)] = min(ra, rb, key=order.index)
    return {f: find(f) for f in cons}


def msa_consensus(rows: list, max_gap: float = MSA_MAX_GAP) -> str:
    arr = np.array([list(r) for r in rows])
    out = []
    for col in arr.T:
        if (col == "-").mean() >= max_gap:
            continue
        vals, cnt = np.unique(col[col != "-"], return_counts=True)
        out.append(vals[np.argmax(cnt)])
    return "".join(out)


def single_hmm(name: str, seq: str):
    bld, bg = pyhmmer.plan7.Builder(ABC), pyhmmer.plan7.Background(ABC)
    hmm, _, _ = bld.build(TextSequence(name=name, sequence=seq).digitize(ABC), bg)
    hmm.name = name
    return hmm


def msa_hmm(name: str, rows: list):
    msa = TextMSA(name=name, sequences=[TextSequence(name=f"m{i}", sequence=r) for i, r in enumerate(rows)])
    bld, bg = pyhmmer.plan7.Builder(ABC), pyhmmer.plan7.Background(ABC)
    hmm, _, _ = bld.build_msa(msa.digitize(ABC), bg)
    hmm.name = name
    return hmm


def mafft(seqs: list, workdir: Path, threads: int) -> list:
    workdir.mkdir(parents=True, exist_ok=True)
    inp = workdir / "mafft_in.fa"
    write_fasta(inp, ((f"m{i}", s) for i, s in enumerate(seqs)))
    out = run(["mafft", "--auto", "--thread", threads, "--quiet", inp]).stdout
    return [s.upper() for _, s in parse_fasta_text(out)]


def seed_intervals(members: pd.DataFrame, p) -> list:
    """Largest merged seed intervals per genome, up to p.hmm_sample_bp per genome."""
    out = []
    for g, d in members.groupby("genome"):
        iv = [(c, s, e) for c, x in d.groupby("chrom") for s, e in merge_intervals(zip(x.start, x.end))]
        iv.sort(key=lambda t: t[1] - t[2])                  # longest first
        budget = p.hmm_sample_bp
        for c, s, e in iv:
            if budget <= 0:
                break
            e = min(e, s + budget)
            out.append((g, c, s, e))
            budget -= e - s
    return out


def pad_hit(lo: int, hi: int, strand: str, hmm_from: int, hmm_to: int, M: int, seqlen: int):
    """Extend a local nhmmer hit (1-based lo..hi) to the full model length, clipped to the sequence.
    Local alignments trim a few model positions at monomer edges; without padding, sampled monomers
    are short and the MSA-derived model loses its end columns (150-bp satellite -> 146-bp model)."""
    left, right = (hmm_from - 1, M - hmm_to) if strand == "+" else (M - hmm_to, hmm_from - 1)
    return max(1, lo - left), min(seqlen, hi + right)


def sample_monomers(hmm, intervals: list, fastas: dict, p, threads: int, rng: random.Random) -> list:
    """nhmmer hmm vs the intervals; up to p.hmm_sample monomers with >= 90% model coverage, significant alone or in
    a run of copies (annotate.run_significant; E-values are for all intervals together), in consensus orientation."""
    raw = []
    for g in dict.fromkeys(i[0] for i in intervals):
        with pysam.FastaFile(str(fastas[g])) as fa:
            raw += [fa.fetch(c, s, e).upper() for gg, c, s, e in intervals if gg == g]
    if not raw:
        return []
    blk = DigitalSequenceBlock(ABC, [TextSequence(name=f"i{i}", sequence=x).digitize(ABC)
                                     for i, x in enumerate(raw)])
    rows = []
    for hits in pyhmmer.hmmer.nhmmer([hmm], blk, cpus=threads, E=NHMMER_REPORT_E):
        for h in hits:
            if h.evalue > NHMMER_REPORT_E:
                continue
            a = h.best_domain.alignment
            lo, hi = sorted((a.target_from, a.target_to))
            rows.append((int(h.name[1:]), lo, hi, h.best_domain.strand, float(h.evalue), a.hmm_from, a.hmm_to))
    hits = pd.DataFrame(rows, columns=["seq", "start", "end", "strand", "evalue", "hmm_from", "hmm_to"])
    hits = hits.assign(family=hmm.name, model_len=hmm.M)
    total = sum(len(x) for x in raw)
    ok = np.zeros(len(hits), bool)
    for _, d in hits.groupby("seq"):
        ok[d.index.to_numpy()] = run_significant(d, total, p.nhmmer_evalue, p.array_gap_factor)
    out = []
    for r in hits[ok & (hits.hmm_to - hits.hmm_from + 1 >= FULL_MONOMER_COV * hmm.M)].itertuples(index=False):
        src = raw[r.seq]
        lo, hi = pad_hit(r.start, r.end, r.strand, r.hmm_from, r.hmm_to, hmm.M, len(src))
        sub = src[lo - 1:hi]
        out.append(sub if r.strand == "+" else revcomp(sub))
    rng.shuffle(out)
    return out[:p.hmm_sample]


def build_hmm(fam_id: str, rep_seq: str, intervals, fastas, workdir: Path, threads: int, p, rng, refs=None):
    """Two-pass family model in the monomer frame of the most similar earlier family (refs: {fam: consensus}).
    Returns (hmm, consensus, mode, n_sampled)."""
    refs = refs or {}
    c0 = frame_rotation(trim_overhang(rep_seq, workdir / "trim0", threads), refs, workdir / "frame0", threads,
                        p.superfam_ident)
    h = single_hmm(fam_id, c0)
    s1 = sample_monomers(h, intervals, fastas, p, threads, rng)
    if len(s1) < p.min_hmm_seqs:
        log.warning(f"families: {fam_id}: only {len(s1)} monomers sampled; single-sequence HMM")
        return h, c0, "single", len(s1)
    rows1 = mafft(s1, workdir, threads)
    c1 = frame_rotation(trim_overhang(msa_consensus(rows1), workdir / "trim1", threads), refs,
                        workdir / "frame1", threads, p.superfam_ident)
    h = single_hmm(fam_id, c1)
    s2 = sample_monomers(h, intervals, fastas, p, threads, rng)
    if len(s2) < p.min_hmm_seqs:                     # a refinement never replaces its model with a worse one:
        log.warning(f"families: {fam_id}: only {len(s2)} monomers in pass 2; profile from pass 1")   # (Arachis:
        h = msa_hmm(fam_id, rows1)                   # pass-1 consensus collapsed to 20 bp, pass 2 found nothing)
        return h, h.consensus.upper(), "profile", len(s1)
    h = msa_hmm(fam_id, mafft(s2, workdir, threads))
    return h, h.consensus.upper(), "profile", len(s2)


def build_families(trf: dict, genomes: list, workdir: Path, threads: int, p) -> pd.DataFrame:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    seeds = select_seeds(trf, p)
    fams, hmms = [], []
    if seeds.empty:
        log.warning("families: no TRF seeds; no satellite families")
    else:
        seqs = dict(zip(seeds.seed_id, seeds.consensus))
        assign = cluster_seeds(seqs, list(seeds.seed_id), workdir / "blast", threads, p.fam_ident, p.multi_ident)
        seeds["rep"] = seeds.seed_id.map(assign)
        seeds["fam_id"] = "."
        fastas = {g.name: g.fasta for g in genomes}
        bases = {r: base_unit(d, fastas)[1] for r, d in seeds.groupby("rep", sort=False)}
        absorbed = absorb_multimers(bases, workdir / "blast_absorb", threads, p.absorb_ident)
        if absorbed:
            seeds["rep"] = seeds.rep.map(lambda r: absorbed.get(r, r))
            log.info(f"families: {len(absorbed)} multimer cluster(s) absorbed into their monomer family")
        spans = {}
        for (rep, g), d in seeds.groupby(["rep", "genome"]):
            spans.setdefault(rep, {})[g] = sum(union_bp(zip(x.start, x.end)) for _, x in d.groupby("chrom"))
        cand = sorted(((max(v.values()), r) for r, v in spans.items() if max(v.values()) >= p.min_fam_span),
                      reverse=True)
        top, dropped = cand[:p.max_families], cand[p.max_families:]
        log.info(f"families: {len(seeds)} seeds -> {len(set(assign.values()))} clusters; "
                 f"{len(top)} carried forward" + (f"; {len(dropped)} more >= min_fam_span dropped by max_families "
                                                  f"(largest {dropped[0][0]} bp)" if dropped else ""))
        refs = {}                                        # frame references: final consensus of earlier families
        for i, (mx, rep) in enumerate(top, 1):
            fid = f"fam{i:03d}"
            mem = seeds[seeds.rep == rep]
            base, unit = base_unit(mem, fastas)
            hmm, cons, mode, n = build_hmm(fid, unit, seed_intervals(mem, p), fastas,     # per-family seed: a model
                                           workdir / fid, threads, p, random.Random(i), refs)   # never depends on others
            refs[fid] = cons
            seeds.loc[mem.index, "fam_id"] = fid
            fams.append(dict(fam_id=fid, rep_seed=rep, base_seed=base.seed_id, consensus_len=hmm.M,
                             consensus=cons, n_seeds=len(mem), max_genome_span=mx,
                             genome_spans=",".join(f"{g}:{v}" for g, v in sorted(spans[rep].items())),
                             hmm_mode=mode, n_hmm_seqs=n))
            hmms.append(hmm)
            log.info(f"families: {fid} {hmm.M} bp, {len(mem)} seeds, max span {mx} bp, hmm={mode} (n={n})")
        if fams:
            merged = dedup_families({f["fam_id"]: f["consensus"] for f in fams}, [f["fam_id"] for f in fams],
                                    workdir / "blast_dedup", threads, p.dedup_ident)
            if merged:
                seeds["fam_id"] = seeds.fam_id.map(lambda x: merged.get(x, x))
                keep = [i for i, f in enumerate(fams) if f["fam_id"] not in merged]
                fams, hmms = [fams[i] for i in keep], [hmms[i] for i in keep]
                for f in fams:
                    mem = seeds[seeds.fam_id == f["fam_id"]]
                    sp = {g: sum(union_bp(zip(x.start, x.end)) for _, x in d.groupby("chrom"))
                          for g, d in mem.groupby("genome")}
                    f.update(n_seeds=len(mem), max_genome_span=max(sp.values()),
                             genome_spans=",".join(f"{g}:{v}" for g, v in sorted(sp.items())))
                log.info(f"families: merged near-identical families: "
                         f"{', '.join(f'{a}->{b}' for a, b in merged.items())}")
            cs = {f["fam_id"]: f["consensus"] for f in fams}
            sf = cluster_seeds(cs, list(cs), workdir / "blast_sf", threads, p.superfam_ident, p.superfam_ident)
            rd = rdna_tags(cs, workdir / "rdna", threads)
            cg = chain_groups(cs, list(cs), workdir / "blast_chain", threads, p.chain_ident)
            for f in fams:
                f["superfamily"], f["rdna"], f["chain_group"] = sf[f["fam_id"]], rd[f["fam_id"]], cg[f["fam_id"]]
            tagged = [f"{f['fam_id']}={f['rdna']}" for f in fams if f["rdna"] != "."]
            log.info(f"families: rDNA-tagged: {', '.join(tagged) or 'none'}")
        seeds.to_csv(workdir / "seeds.tsv.gz", sep="\t", index=False)
    table = pd.DataFrame(fams, columns=FAM_COLS)
    table.to_csv(workdir / "families.tsv", sep="\t", index=False)
    write_fasta(workdir / "families.fa", ((f["fam_id"], f["consensus"]) for f in fams))
    with open(workdir / "families.hmm", "wb") as fh:
        for h in hmms:
            h.write(fh)
    return table
