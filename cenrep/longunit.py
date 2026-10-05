"""S5b: tandem arrays with units > 2 kb (beyond TRF's MaxPeriod), from windowed minimap2 self-alignment.
A tandem array of unit d self-aligns with a shift of d; an alignment at least d long means >= 2 adjacent copies."""
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import pysam

from .util import log, merge_intervals, run

LU_COLS = ["chrom", "start", "end", "unit", "copies", "identity"]


def windows(start: int, end: int, size: int, overlap: int):
    s = start
    while s < end:
        e = min(end, s + size)
        yield s, e
        if e == end:
            break
        s = e - overlap


def build_masks(frames) -> dict:
    """Merged intervals to N-mask (satellite arrays, TRF rows): {chrom: (starts, ends)}."""
    iv = {}
    for df in frames:
        for c, s, e in zip(df.chrom, df.start, df.end):
            iv.setdefault(c, []).append((int(s), int(e)))
    out = {}
    for c, x in iv.items():
        m = np.array(merge_intervals(x), dtype=np.int64).reshape(-1, 2)
        out[c] = (m[:, 0], m[:, 1])
    return out


def mask_seq(seq: str, offset: int, starts, ends) -> str:
    b = bytearray(seq.encode())
    n = len(b)
    i0 = int(np.searchsorted(ends, offset, side="right"))
    i1 = int(np.searchsorted(starts, offset + n, side="left"))
    for s, e in zip(starts[i0:i1], ends[i0:i1]):
        a, z = max(int(s) - offset, 0), min(int(e) - offset, n)
        if z > a:
            b[a:z] = b"N" * (z - a)
    return b.decode()


def selfaln_hits(seq: str, offset: int, tmpdir, p) -> list:
    fd, path = tempfile.mkstemp(suffix=".fa", dir=tmpdir)
    with os.fdopen(fd, "w") as f:
        f.write(f">w\n{seq}\n")
    try:
        out = run(["minimap2", "-x", "asm20", "-c", "-DP", "--dual=no", "-t", "1", path, path]).stdout
    finally:
        os.unlink(path)
    hits = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12 or f[4] != "+":
            continue
        qs, qe, ts, te, nm, al = int(f[2]), int(f[3]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
        d = ts - qs
        if d <= 0 or not p.lu_min_unit <= d <= p.lu_max_unit or qe - qs < d or nm / al < p.lu_min_ident:
            continue
        hits.append((offset + qs, offset + te, d, nm / al))
    return hits


def merge_hits(hits, chrom: str) -> list:
    """Overlapping hits form a group; within it each hit joins the smallest base unit it is a multiple of
    (5%), and each base unit's spans are unioned into arrays."""
    out, groups, cur, cur_end = [], [], [], -1
    for h in sorted(hits):
        if cur and h[0] >= cur_end:
            groups.append(cur)
            cur, cur_end = [], -1
        cur.append(h)
        cur_end = max(cur_end, h[1])
    if cur:
        groups.append(cur)
    for g in groups:
        bases = {}
        for h in sorted(g, key=lambda x: x[2]):
            for b in bases:
                k = round(h[2] / b)
                if k >= 1 and abs(h[2] / b - k) <= 0.05 * k:
                    bases[b].append(h)
                    break
            else:
                bases[h[2]] = [h]
        for b, hs in bases.items():
            unit_hits = [h for h in hs if round(h[2] / b) == 1]
            unit = int(np.median([h[2] for h in unit_hits]))
            ident = float(np.mean([h[3] for h in unit_hits]))
            for s, e in merge_intervals((h[0], h[1]) for h in hs):
                out.append(dict(chrom=chrom, start=s, end=e, unit=unit,
                                copies=round((e - s) / unit, 2), identity=round(ident, 4)))
    return out


def detect(fasta, targets, masks: dict, tmpdir, threads: int, p) -> pd.DataFrame:
    jobs = [(c, ws, we) for c, s, e in targets for ws, we in windows(s, e, p.lu_window, p.lu_window_overlap)]
    empty = (np.zeros(0, np.int64), np.zeros(0, np.int64))

    def work(job):
        c, ws, we = job
        with pysam.FastaFile(str(fasta)) as fa:
            seq = mask_seq(fa.fetch(c, ws, we).upper(), ws, *masks.get(c, empty))
        if len(seq) - seq.count("N") < 2 * p.lu_min_unit:
            return c, []
        return c, selfaln_hits(seq, ws, tmpdir, p)

    per = {}
    with ThreadPoolExecutor(max(1, threads)) as ex:
        for c, h in ex.map(work, jobs):
            per.setdefault(c, []).extend(h)
    rows = [r for c, h in per.items() for r in merge_hits(h, c)]
    log.debug(f"longunit: {len(jobs)} windows, {len(rows)} arrays")
    return pd.DataFrame(rows, columns=LU_COLS)
