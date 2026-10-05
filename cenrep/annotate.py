"""S3: genome-wide nhmmer monomer calls (chromosome chunks on a thread pool), overlap resolution, array chaining."""
import threading
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import pyhmmer
import pysam
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence

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


def search_chrom(hmms, chrom: str, seq: str, threads: int, evalue: float) -> list:
    if not hmms:
        return []
    blk = DigitalSequenceBlock(ABC, [TextSequence(name=chrom, sequence=seq.upper()).digitize(ABC)])
    mlen = {h.name: h.M for h in hmms}
    rows = []
    for hits in pyhmmer.hmmer.nhmmer(hmms, blk, cpus=threads, E=evalue):
        fam = hits.query.name
        for h in hits:
            if h.evalue > evalue:
                continue
            d = h.best_domain
            a = d.alignment
            lo, hi = sorted((a.target_from, a.target_to))
            rows.append((chrom, lo - 1, hi, d.strand, fam, float(h.score), float(h.evalue),
                         a.hmm_from, a.hmm_to, mlen[fam], _identity(a)))
    return rows


def _chunks(chroms: dict, size: int, overlap: int):
    for c, n in chroms.items():
        s = 0
        while s < n:
            e = min(n, s + size)
            yield c, s, e
            if e == n:
                break
            s = e - overlap


def search_genome(hmms, genome, threads: int, evalue: float, chunk: int = 2_000_000, overlap: int = 10_000,
                  ref_len: int = 0) -> dict:
    """{chrom: overlap-resolved monomer hits}. pyhmmer's nhmmer parallelises over query models only, so a few
    long models would leave most threads idle; instead every (chunk x all models) job runs on a thread pool
    (pyhmmer releases the GIL). Chunk E-values are rescaled by chromosome/chunk length, so the cutoff equals
    that of one whole-chromosome search. overlap must exceed the longest model; duplicate and edge-truncated
    hits in overlaps are removed by resolve_overlaps. Sequences shorter than ref_len are scaled as if ref_len long,
    so short contigs get the cutoff of a chromosome instead of a looser one."""
    out = {c: [] for c in genome.chroms}
    if hmms:
        mlen = {h.name: h.M for h in hmms}
        local = threading.local()

        def work(job):
            c, s, e = job
            if not hasattr(local, "fa"):
                local.fa = pysam.FastaFile(str(genome.fasta))
            scale = max(genome.chroms[c], ref_len) / (e - s)
            blk = DigitalSequenceBlock(ABC, [TextSequence(name=c, sequence=local.fa.fetch(c, s, e).upper())
                                             .digitize(ABC)])
            rows = []
            for hits in pyhmmer.hmmer.nhmmer(hmms, blk, cpus=1, E=evalue):
                fam = hits.query.name
                for h in hits:
                    if h.evalue * scale > evalue:
                        continue
                    d = h.best_domain
                    a = d.alignment
                    lo, hi = sorted((a.target_from, a.target_to))
                    rows.append((c, s + lo - 1, s + hi, d.strand, fam, float(h.score), float(h.evalue * scale),
                                 a.hmm_from, a.hmm_to, mlen[fam], _identity(a)))
            return c, rows

        overlap = max(overlap, 2 * max(mlen.values()))
        with ThreadPoolExecutor(max(1, threads)) as ex:
            for c, rows in ex.map(work, _chunks(genome.chroms, chunk, overlap)):
                out[c].extend(rows)
    return {c: resolve_overlaps(drop_multimer_hits(pd.DataFrame(r, columns=MON_COLS))) for c, r in out.items()}


def drop_multimer_hits(df: pd.DataFrame, ratio: float = 1.9, min_cov: float = 0.5) -> pd.DataFrame:
    """Describe satellites by their shortest repeat unit: drop a hit when hits of families whose model is
    >= ratio x shorter cover >= min_cov of it. Composite or non-integer multimer families (CentBd x7 + spacer,
    CentO x2.3) otherwise outscore the monomer hits on the monomer arrays, because bitscore grows with length."""
    if df.empty:
        return df
    df = df.sort_values(["start", "end"]).reset_index(drop=True)
    s, e, m = df.start.to_numpy(), df.end.to_numpy(), df.model_len.to_numpy()
    maxlen = int((e - s).max())
    keep = np.ones(len(df), bool)
    for i in np.where(m >= ratio * m.min())[0]:          # only hits of long models can be multimers
        lo, hi = np.searchsorted(s, s[i] - maxlen), np.searchsorted(s, e[i])
        iv = [(max(s[j], s[i]), min(e[j], e[i])) for j in range(lo, hi)
              if j != i and m[j] * ratio <= m[i] and e[j] > s[i]]
        if iv and union_bp(iv) >= min_cov * (e[i] - s[i]):
            keep[i] = False
    return df[keep].reset_index(drop=True)


def resolve_overlaps(df: pd.DataFrame, max_cov: float = 0.2) -> pd.DataFrame:
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
    per = search_genome(hmms, genome, threads, p.nhmmer_evalue, ref_len=p.evalue_ref_len)
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
