"""S6: monomer variant classes (cd-hit-est) and higher-order repeats (HORs) via maximal repeats.

Each centromere becomes an integer string (one symbol per monomer class, consensus orientation, unique
negative separators between arrays and for unclassed monomers). A HOR unit is a maximal repeat of >= 2
monomers with >= 2 distinct classes that does not contain two full copies of a shorter period."""
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from .util import log, revcomp, run, write_fasta

MONHOR_COLS = ["mon_idx", "variant_class", "hor_unit", "hor_occ"]
OCC_COLS = ["occ_id", "unit", "chrom", "start", "end", "strand", "array_id", "n_monomers", "tandem_run",
            "tandem_copies"]
UNIT_COLS = ["unit", "family", "chrom", "n_monomers", "occurrences", "max_tandem_copies", "class_string"]


def suffix_array(t) -> np.ndarray:
    """Prefix doubling, O(n log^2 n) with numpy sorts. Shorter suffix sorts first on ties."""
    t = np.asarray(t)
    n = len(t)
    if n == 0:
        return np.zeros(0, np.int64)
    rank = np.unique(t, return_inverse=True)[1].astype(np.int64).ravel()
    k = 1
    while True:
        r2 = np.full(n, -1, np.int64)
        if k < n:
            r2[:n - k] = rank[k:]
        sa = np.lexsort((r2, rank))
        a, b = rank[sa], r2[sa]
        new = np.zeros(n, np.int64)
        new[1:] = np.cumsum((a[1:] != a[:-1]) | (b[1:] != b[:-1]))
        rank = np.empty(n, np.int64)
        rank[sa] = new
        if new[-1] == n - 1:
            return sa
        k *= 2


def lcp_array(t, sa) -> np.ndarray:
    """Kasai: lcp[i] = LCP(suffix sa[i-1], suffix sa[i]); lcp[0] = 0."""
    tl, sal = [int(x) for x in t], [int(x) for x in sa]
    n = len(tl)
    rank = [0] * n
    for i, s in enumerate(sal):
        rank[s] = i
    lcp, h = [0] * n, 0
    for i in range(n):
        r = rank[i]
        if r == 0:
            h = 0
            continue
        j = sal[r - 1]
        while i + h < n and j + h < n and tl[i + h] == tl[j + h]:
            h += 1
        lcp[r] = h
        if h:
            h -= 1
    return np.array(lcp, np.int64)


def maximal_repeats(t, sa, lcp, max_len: int, min_len: int = 2, keep=None) -> list:
    """Bottom-up lcp-interval traversal (Abouelhoda et al. 2004). A repeat is reported if its interval
    is left-diverse. lcp is capped at max_len, so the work is O(n * max_len)."""
    t = np.asarray(t)
    n = len(t)
    if n < 2:
        return []
    lcp = np.minimum(lcp, max_len)
    prev = np.where(sa > 0, t[np.maximum(sa - 1, 0)], -(10 ** 15) - np.arange(n))
    out, stack = [], [[0, 0]]
    for i in range(1, n + 1):
        cur = int(lcp[i]) if i < n else 0
        lb = i - 1
        while cur < stack[-1][0]:
            l, lb = stack.pop()
            if l >= min_len:
                if keep is None or keep(tuple(int(x) for x in t[sa[lb]:sa[lb] + l])):
                    pv = prev[lb:i]
                    if pv.min() != pv.max():
                        out.append((l, sa[lb:i]))
        if cur > stack[-1][0]:
            stack.append([cur, lb])
    return out


def smallest_period(u) -> int:
    n = len(u)
    if n == 0:
        return 0
    f, k = [0] * n, 0
    for i in range(1, n):
        while k and u[i] != u[k]:
            k = f[k - 1]
        if u[i] == u[k]:
            k += 1
        f[i] = k
    return n - f[-1]


def is_hor_unit(u) -> bool:
    return len(u) >= 2 and min(u) >= 0 and len(set(u)) >= 2 and 2 * smallest_period(u) > len(u)


def _greedy_count(occ, l) -> int:
    n, last = 0, -(10 ** 18)
    for q in occ:
        if q >= last + l:
            n, last = n + 1, q
    return n


def hor_candidate(w) -> bool:
    """Keep a maximal repeat if it, or its period unit, is a valid HOR unit."""
    return is_hor_unit(w) or is_hor_unit(w[:smallest_period(w)])


def assign_hors(t, reps) -> list:
    """Candidates are valid maximal repeats plus the period unit of periodic maximal repeats: a tandem HOR
    array ending in a partial copy (abcd x 14 + abc) has no maximal 'abcd' (every copy is followed by 'a');
    its unit is the period of the maximal repeat 'abcdabc', occurring at each full period of each occurrence.
    Rank units by achievable non-overlapping coverage, then length; accept occurrences whose monomers
    are unassigned; keep units with >= 2 accepted occurrences."""
    t = np.asarray(t)
    occs = {}
    for l, occ in reps:
        occ = sorted(int(x) for x in occ)
        w = tuple(int(x) for x in t[occ[0]:occ[0] + l])
        if is_hor_unit(w):
            occs.setdefault(w, set()).update(occ)
        per = smallest_period(w)
        if per < l and is_hor_unit(w[:per]):
            occs.setdefault(w[:per], set()).update(q + k * per for q in occ for k in range(l // per))
    cands = []
    for u, occ in occs.items():
        occ = sorted(occ)
        c = _greedy_count(occ, len(u))
        if c >= 2:
            cands.append((c * len(u), len(u), occ, u))
    cands.sort(key=lambda x: (-x[0], -x[1], x[2][0]))
    used = np.zeros(len(t), bool)
    units = []
    for _, l, occ, u in cands:
        acc = []
        for q in occ:
            if not used[q:q + l].any():
                used[q:q + l] = True
                acc.append(q)
        if len(acc) >= 2:
            units.append((u, acc))
        else:
            for q in acc:
                used[q:q + l] = False
    return units


def tandem_runs(acc, l) -> list:
    runs = [[acc[0]]]
    for q in acc[1:]:
        if q == runs[-1][-1] + l:
            runs[-1].append(q)
        else:
            runs.append([q])
    return runs


def variant_classes(seqs: dict, workdir, threads: int, ident: float) -> dict:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    fa, out = workdir / "mon.fa", workdir / "mon.cdhit"
    write_fasta(fa, seqs.items())
    run(["cd-hit-est", "-i", fa, "-o", out, "-c", ident, "-n", 10 if ident >= 0.95 else 8, "-G", 1,
         "-aS", 0.9, "-aL", 0.9, "-g", 1, "-r", 0, "-d", 0, "-M", 0, "-T", threads])
    cl, res = -1, {}
    with open(str(out) + ".clstr") as f:
        for line in f:
            if line.startswith(">Cluster"):
                cl = int(line.split()[1])
            else:
                res[line.split(">", 1)[1].split("...", 1)[0]] = cl
    return res


def run_hor(genome, mon: pd.DataFrame, arrays: pd.DataFrame, regions: dict, primary: set, workdir,
            threads: int, p):
    workdir = Path(workdir)
    assert (mon.mon_idx.to_numpy() == np.arange(len(mon))).all(), "mon_idx must equal row position"
    cov = (mon.hmm_to - mon.hmm_from + 1) / mon.model_len
    sel = mon[mon.family.isin(primary) & (cov >= p.hor_min_cov)]
    seqs = {}
    with pysam.FastaFile(str(genome.fasta)) as fa:
        for c, g in sel.groupby("chrom", sort=False):
            cs = fa.fetch(c).upper()
            for i, s, e, st in zip(g.mon_idx, g.start, g.end, g.strand):
                x = cs[s:e]
                seqs[int(i)] = x if st == "+" else revcomp(x)
    klass = {}
    for f in sorted(primary):
        ids = sel.mon_idx[sel.family == f].astype(int).tolist()
        if ids:
            for name, k in variant_classes({f"m{i}": seqs[i] for i in ids}, workdir / f"cdhit_{f}",
                                           threads, p.hor_ident).items():
                klass[int(name[1:])] = (f, k)
    code = {key: n for n, key in enumerate(sorted(set(klass.values())))}
    inv = {v: k for k, v in code.items()}
    m_start, m_end = mon.start.to_numpy(), mon.end.to_numpy()
    by_array = {a: g.mon_idx.astype(int).tolist() for a, g in mon.groupby("array_id", sort=False)}
    hor_unit, hor_occ, occ_rows, unit_rows = {}, {}, [], []
    for c, r in regions.items():
        if r["evidence"] == "none":
            continue
        a = arrays[(arrays.chrom == c) & arrays.family.isin(primary)
                   & (arrays.start >= r["cen_start"]) & (arrays.end <= r["cen_end"])].sort_values("start")
        t, pos2mon, pos2arr, sep = [], [], [], -1
        for aid, st in zip(a.array_id, a.strand):
            idxs = by_array.get(aid, [])
            for i in (idxs[::-1] if st == "-" else idxs):
                if i in klass:
                    t.append(code[klass[i]])
                else:
                    t.append(sep)
                    sep -= 1
                pos2mon.append(i)
                pos2arr.append((aid, st))
            t.append(sep)
            sep -= 1
            pos2mon.append(-1)
            pos2arr.append(None)
        if sum(x >= 0 for x in t) < 4:
            continue
        t = np.array(t, np.int64)
        sa = suffix_array(t)
        units = assign_hors(t, maximal_repeats(t, sa, lcp_array(t, sa), p.hor_max_len, keep=hor_candidate))
        units.sort(key=lambda x: (-len(x[0]) * len(x[1]), x[1][0]))
        n_occ = 0
        for k, (u, acc) in enumerate(units, 1):
            uid, l = f"HOR_{c}_{k}", len(u)
            runs = tandem_runs(acc, l)
            for ri, rr in enumerate(runs, 1):
                for q in rr:
                    n_occ += 1
                    oid = f"hor_{c}_{n_occ}"
                    mids = pos2mon[q:q + l]
                    aid, st = pos2arr[q]
                    occ_rows.append((oid, uid, c, int(m_start[mids].min()), int(m_end[mids].max()), st, aid, l,
                                     f"{uid}.r{ri}", len(rr)))
                    for mi in mids:
                        hor_unit[mi], hor_occ[mi] = uid, oid
            unit_rows.append((uid, inv[u[0]][0], c, l, len(acc), max(len(rr) for rr in runs),
                              "-".join(str(inv[x][1]) for x in u)))
        log.debug(f"hor: {genome.name} {c}: {len(t)} symbols, {len(units)} HOR units")
    mh = pd.DataFrame([(i, k[1], hor_unit.get(i, "."), hor_occ.get(i, ".")) for i, k in sorted(klass.items())],
                      columns=MONHOR_COLS)
    return mh, pd.DataFrame(occ_rows, columns=OCC_COLS), pd.DataFrame(unit_rows, columns=UNIT_COLS)
