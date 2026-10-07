"""S3: genome-wide nhmmer monomer calls (chromosome chunks on a thread pool), overlap resolution, array chaining."""
import threading
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import pyhmmer
import pysam
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence

from .params import ANNOT_CHUNK, FULL_MONOMER_COV, MULTIMER_MIN_COV, MULTIMER_RATIO, NHMMER_REPORT_E, OVERLAP_MAX_COV
from .util import log, union_bp

ABC = Alphabet.dna()
MON_COLS = ["chrom", "start", "end", "strand", "family", "bitscore", "evalue",
            "hmm_from", "hmm_to", "model_len", "identity"]
ARR_COLS = ["array_id", "chrom", "start", "end", "strand", "family", "n_monomers", "mean_identity", "monomer_len"]


def _identity(a) -> float:
    """Exact matches / alignment columns (HMM consensus vs target, case-insensitive)."""
    h = np.frombuffer(a.hmm_sequence.upper().encode(), np.uint8)
    t = np.frombuffer(a.target_sequence.upper().encode(), np.uint8)
    return float((h == t).mean())


def run_significant(df: pd.DataFrame, L: int, evalue: float, gap_factor: float) -> np.ndarray:
    """Which hits are significant alone or as part of a run of whole copies. A copy of a short unit cannot be
    significant alone (a 40-bp model scores <= ~20 bits; E 1e-5 on a 500-Mb chromosome needs ~47), but a tandem run
    of copies is no chance event: after one hit, a random hit starts within the array gap W = gap_factor x unit + 1
    with probability ~W x E_i / L. A run of r whole copies (>= FULL_MONOMER_COV of the model; same family and strand;
    gaps <= gap_factor x unit, as chain_arrays joins monomers) therefore has E = prod(E_i) x (W / L)^(r - 1).
    Fragments count only on their own E: weak partial hits of long models cluster in repeats (not independent), and
    an array is made of whole copies. df: one sequence of length L, E-values for that length (evalue), hmm_from,
    hmm_to, unit in model_len. Returns a mask aligned to df: own E or run E <= evalue."""
    if df.empty:
        return np.zeros(0, bool)
    alone = (np.log10(np.maximum(df.evalue.to_numpy(), 1e-300)) <= np.log10(evalue))
    whole = (df.hmm_to - df.hmm_from + 1).to_numpy() >= FULL_MONOMER_COV * df.model_len.to_numpy()
    out = alone.copy()
    if not whole.any():
        return out
    d = df.loc[whole, ["family", "strand", "start", "end", "evalue", "model_len"]].reset_index(drop=True)
    d = d.assign(_i=np.flatnonzero(whole)).sort_values(["family", "strand", "start"], kind="stable")
    s, e = d.start.to_numpy(), d.end.to_numpy()
    unit = d.model_len.to_numpy().astype(float)
    key = (d.family.astype(str) + d.strand.astype(str)).to_numpy()
    same = np.concatenate([[False], key[1:] == key[:-1]])
    reach = pd.Series(e).groupby(np.cumsum(~same)).cummax().to_numpy()         # furthest end so far in each key
    brk = ~same | np.concatenate([[True], s[1:] - reach[:-1] > gap_factor * unit[1:]])
    run = np.cumsum(brk)
    loge = np.log10(np.maximum(d.evalue.to_numpy(), 1e-300))
    g = pd.DataFrame({"run": run, "loge": loge, "w": np.log10((gap_factor * unit + 1) / L)})
    tot = g.groupby("run").agg(s=("loge", "sum"), n=("loge", "size"), w=("w", "first"))
    run_e = (tot.s + (tot.n - 1) * tot.w).to_numpy()[run - 1]
    out[d._i.to_numpy()] |= run_e <= np.log10(evalue)
    return out


def _rows(hits, mlen: dict, chrom: str, offset: int, scale: float) -> list:
    """MON_COLS rows of the hits HMMER would report (E <= NHMMER_REPORT_E per search; pyhmmer's TopHits also holds
    weaker ones); E-values scaled to the chromosome."""
    rows = []
    for hs in hits:
        fam = hs.query.name
        for h in hs:
            if h.evalue > NHMMER_REPORT_E:
                continue
            d = h.best_domain
            a = d.alignment
            lo, hi = sorted((a.target_from, a.target_to))
            rows.append((chrom, offset + lo - 1, offset + hi, d.strand, fam, float(h.score), float(h.evalue * scale),
                         a.hmm_from, a.hmm_to, mlen[fam], _identity(a)))
    return rows


def significant(rows: list, L: int, evalue: float, gap_factor: float) -> pd.DataFrame:
    """Hits of one sequence (E-values for length L) that are significant alone or in a run; chunk-overlap duplicates
    are dropped first so that no copy counts twice."""
    df = pd.DataFrame(rows, columns=MON_COLS).drop_duplicates(["start", "end", "strand", "family"])
    return df[run_significant(df, L, evalue, gap_factor)].reset_index(drop=True)


def search_chrom(hmms, chrom: str, seq: str, threads: int, evalue: float, gap_factor: float = 2.0) -> list:
    if not hmms:
        return []
    blk = DigitalSequenceBlock(ABC, [TextSequence(name=chrom, sequence=seq.upper()).digitize(ABC)])
    rows = _rows(pyhmmer.hmmer.nhmmer(hmms, blk, cpus=threads, E=NHMMER_REPORT_E), {h.name: h.M for h in hmms},
                 chrom, 0, 1.0)
    return list(significant(rows, len(seq), evalue, gap_factor).itertuples(index=False, name=None))


def _chunks(chroms: dict, size: int, overlap: int):
    for c, n in chroms.items():
        s = 0
        while s < n:
            e = min(n, s + size)
            yield c, s, e
            if e == n:
                break
            s = e - overlap


def search_genome(hmms, genome, threads: int, evalue: float, chunk: int = ANNOT_CHUNK, overlap: int = 10_000,
                  ref_len: int = 0, gap_factor: float = 2.0) -> dict:
    """{chrom: overlap-resolved monomer hits}. pyhmmer's nhmmer parallelises over query models only, so a few
    long models would leave most threads idle; instead every (chunk x model) job runs on a thread pool (pyhmmer
    releases the GIL). One job per model, not per chunk: a satellite-dense 2-Mb chunk searched with all ~40 models
    took > 40 min on one thread while the other 24 idled (Arachis); hits and E-values are the same, since each
    model's search is independent. Chunk E-values are rescaled by chromosome/chunk length, so the cutoff equals
    that of one whole-chromosome search. overlap must exceed the longest model; duplicate and edge-truncated
    hits in overlaps are removed by resolve_overlaps. Sequences shorter than ref_len are scaled as if ref_len long,
    so short contigs get the cutoff of a chromosome instead of a looser one. Hits are kept if significant alone or
    as part of a run of copies (run_significant), so short units are found too."""
    out = {c: [] for c in genome.chroms}
    if hmms:
        mlen = {h.name: h.M for h in hmms}
        local = threading.local()

        def work(job):
            c, s, e, h = job
            if not hasattr(local, "fa"):
                local.fa, local.key = pysam.FastaFile(str(genome.fasta)), None
            if local.key != (c, s, e):                    # the thread's last chunk is reused when it comes again
                local.blk = DigitalSequenceBlock(ABC, [TextSequence(name=c, sequence=local.fa.fetch(c, s, e).upper())
                                                       .digitize(ABC)])
                local.key = (c, s, e)
            scale = max(genome.chroms[c], ref_len) / (e - s)
            return c, _rows(pyhmmer.hmmer.nhmmer([h], local.blk, cpus=1, E=NHMMER_REPORT_E), mlen, c, s, scale)

        overlap = max(overlap, 2 * max(mlen.values()))
        jobs = ((c, s, e, h) for c, s, e in _chunks(genome.chroms, chunk, overlap) for h in hmms)
        with ThreadPoolExecutor(max(1, threads)) as ex:
            for c, rows in ex.map(work, jobs):
                out[c].extend(rows)
    return {c: resolve_overlaps(drop_multimer_hits(significant(r, max(genome.chroms[c], ref_len), evalue, gap_factor),
                                                   evalue=evalue))
            for c, r in out.items()}


def _union(s: np.ndarray, e: np.ndarray):
    """Merged union of intervals sorted by start: (starts, ends, bp before each interval)."""
    if not len(s):
        return s, e, np.zeros(1, np.int64)
    run = np.maximum.accumulate(e)
    first = np.concatenate([[True], s[1:] > run[:-1]])
    us = s[first]
    ue = np.maximum.reduceat(e, np.flatnonzero(first))
    return us, ue, np.concatenate([[0], np.cumsum(ue - us)])


def _covered_below(x: np.ndarray, us, ue, cum) -> np.ndarray:
    """bp of the union below each position x."""
    k = np.searchsorted(us, x, side="left")
    over = np.where(k > 0, np.maximum(ue[np.maximum(k - 1, 0)] - x, 0), 0)
    return cum[k] - over


def drop_multimer_hits(df: pd.DataFrame, ratio: float = MULTIMER_RATIO, min_cov: float = MULTIMER_MIN_COV,
                       evalue: float = None) -> pd.DataFrame:
    """Describe satellites by their shortest repeat unit: drop a hit when hits of families whose model is
    >= ratio x shorter cover >= min_cov of it. Composite or non-integer multimer families (CentBd x7 + spacer,
    CentO x2.3) otherwise outscore the monomer hits on the monomer arrays, because bitscore grows with length.
    With evalue, only shorter hits significant on their own count: copies significant only as a run
    (run_significant) fill sequence nothing else explains but never erase a hit that is significant alone (Avena:
    run-only 38-41-bp copies erased 87/130-bp arrays and then lost the overlap contest, leaving sequence unannotated).
    Vectorised: per distinct long-model length, the union of the shorter-model hits is built once and each hit's
    coverage read from prefix sums (Avena C genomes: tens of millions of hits)."""
    if df.empty:
        return df
    df = df.sort_values(["start", "end"]).reset_index(drop=True)
    s, e, m = df.start.to_numpy(), df.end.to_numpy(), df.model_len.to_numpy()
    strong = np.ones(len(df), bool) if evalue is None else df.evalue.to_numpy() <= evalue
    keep = np.ones(len(df), bool)
    cand = m >= ratio * m.min()                          # only hits of long models can be multimers
    for t in np.unique(m[cand]):
        short = (m * ratio <= t) & strong
        if not short.any():
            continue
        us, ue, cum = _union(s[short], e[short])
        idx = np.flatnonzero(cand & (m == t))
        cov = _covered_below(e[idx], us, ue, cum) - _covered_below(s[idx], us, ue, cum)
        keep[idx[(cov > 0) & (cov >= min_cov * (e[idx] - s[idx]))]] = False
    return df[keep].reset_index(drop=True)


def resolve_overlaps(df: pd.DataFrame, max_cov: float = OVERLAP_MAX_COV) -> pd.DataFrame:
    """Greedy by bitscore: accept a hit only if <= max_cov of it is already covered by accepted hits
    (coverage bitmap over the chromosome span). A pairwise ">50% of the shorter" rule let hits of phase-shifted
    subfamilies that overlap each of two accepted neighbours by <= 50% survive (A. thaliana CEN180: +27% monomers)."""
    if df.empty:
        return df
    df = df.sort_values(["start", "end"]).reset_index(drop=True)
    s, e = df.start.to_numpy(), df.end.to_numpy()
    lo = int(s.min())
    occ = np.zeros(int(e.max()) - lo, bool)
    keep = np.zeros(len(df), bool)
    for i in np.argsort(-df.bitscore.to_numpy(), kind="stable"):
        a, b = s[i] - lo, e[i] - lo
        if occ[a:b].sum() <= max_cov * (b - a):
            occ[a:b] = True
            keep[i] = True
    return df[keep].reset_index(drop=True)


def chain_arrays(mon: pd.DataFrame, fam_len: dict, p, groups: dict = None):
    """Monomers on one strand whose families share a chain group (close subfamilies; default: the family itself)
    with gap <= array_gap_factor x the longest member model form an array; arrays need >= min_array_monomers.
    An array's family is its majority family by monomer bp; monomers keep their own family.
    Returns (monomers in arrays with array_id, arrays)."""
    if mon.empty:
        return mon.assign(array_id=pd.Series(dtype="int64")), pd.DataFrame(columns=ARR_COLS)
    groups = groups or {}
    mon = mon.assign(_grp=mon.family.map(lambda f: groups.get(f, f)))
    mon = mon.sort_values(["chrom", "_grp", "strand", "start"]).reset_index(drop=True)
    ids = np.full(len(mon), -1, np.int64)
    rows, aid = [], 0
    ident, fam, blen = mon.identity.to_numpy(), mon.family.to_numpy(), (mon.end - mon.start).to_numpy()
    for (c, _, st), g in mon.groupby(["chrom", "_grp", "strand"], sort=False):
        idx = g.index.to_numpy()
        s, e = g.start.to_numpy(), g.end.to_numpy()
        maxgap = p.array_gap_factor * max(fam_len[f] for f in set(fam[idx]))
        brk = np.where(s[1:] - np.maximum.accumulate(e)[:-1] > maxgap)[0] + 1
        for block in np.split(np.arange(len(idx)), brk):
            if len(block) < p.min_array_monomers:
                continue
            ii = idx[block]
            ids[ii] = aid
            major = pd.Series(blen[ii]).groupby(fam[ii]).sum().idxmax()
            rows.append((aid, c, int(s[block].min()), int(e[block].max()), st, major, len(block),
                         float(ident[ii].mean()), int(fam_len[major])))
            aid += 1
    mon["array_id"] = ids
    return mon[mon.array_id >= 0].drop(columns="_grp").copy(), pd.DataFrame(rows, columns=ARR_COLS)


def annotate_genome(genome, hmms, fam_len: dict, threads: int, p, out_mon, out_arr, groups: dict = None):
    per = search_genome(hmms, genome, threads, p.nhmmer_evalue, ref_len=p.evalue_ref_len, gap_factor=p.array_gap_factor)
    parts = [per[c] for c in genome.chroms if len(per[c])]
    log.debug(f"annotate: {genome.name}: monomer hits on {len(parts)}/{len(genome.chroms)} sequences")
    mon = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=MON_COLS)
    mon, arrays = chain_arrays(mon, fam_len, p, groups)
    order = {c: i for i, c in enumerate(genome.chroms)}
    mon = mon.assign(ci=mon.chrom.map(order)).sort_values(["ci", "start", "end"]).drop(columns="ci")
    mon = mon.reset_index(drop=True)
    mon.insert(0, "mon_idx", np.arange(len(mon), dtype=np.int64))
    arrays = arrays.assign(ci=arrays.chrom.map(order)).sort_values(["ci", "start"]).drop(columns="ci")
    mon.to_csv(out_mon, sep="\t", index=False)
    arrays.to_csv(out_arr, sep="\t", index=False)
    log.info(f"annotate: {genome.name}: {len(mon)} monomers in {len(arrays)} arrays")
    return mon, arrays.reset_index(drop=True)
