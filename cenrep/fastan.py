"""S1 engine: FasTAN (Myers) tandem arrays, with each array's interval, mean unit, divergence of adjacent copies and
the cuts between copies. FasTAN reports no unit sequence, so a seed's consensus is the partial-order alignment (abPOA)
of its copies."""
import os
import random
import shutil
import tempfile

import numpy as np
import pyabpoa

from .params import CONS_UNITS, FASTAN_ARGS, LEN_TOL, RELATED_LEN_TOL, SEED_MAX_PERIOD
from .util import CenrepError, log, run

TOOLS = ["FAtoGDB", "FasTAN", "ONEview"]
DUP_TOL = LEN_TOL   # a unit is a mean over the copies aligned, so the same array seen from two chunks can differ
KEY = ["fastan", 1, FASTAN_ARGS, CONS_UNITS, RELATED_LEN_TOL, SEED_MAX_PERIOD, DUP_TOL]   # cache key; bump 1 if rows change
NAN = float("nan")
_POA = None


def parse_ano(text: str) -> list:
    """ONEview text of a FasTAN .1ano -> [(beg, end, unit, pct_diff, cuts)], 0-based half-open; cuts are the copy
    boundaries relative to beg (P line), none for a near tandem (two copies with a short gap)."""
    rows = []
    for line in text.splitlines():
        t = line.split()
        if not t:
            continue
        try:
            if t[0] == "M":
                rows.append([int(t[2]), int(t[3]), 0, 0, []])
            elif rows and t[0] == "L":
                rows[-1][2] = int(t[2])
            elif rows and t[0] == "X":
                rows[-1][3] = int(t[1])
            elif rows and t[0] == "P":
                rows[-1][4] = [int(x) for x in t[2:2 + int(t[1])]]
        except (ValueError, IndexError):                 # a format change, not a bad record: stop
            raise CenrepError(f"FasTAN output: unparsable .1ano line {line[:80]!r} (FasTAN/ONEview version?)") from None
    return [tuple(r) for r in rows]


def run_fastan(seq: str, tmpdir) -> str:
    """FasTAN on one sequence -> the ONEview text of its .1ano (FasTAN masks only from a GDB, hence FAtoGDB)."""
    work = tempfile.mkdtemp(dir=tmpdir)
    try:
        with open(os.path.join(work, "s.fa"), "w") as f:
            f.write(f">s\n{seq}\n")
        run(["FAtoGDB", "s.fa"], cwd=work)
        run(["FasTAN", *FASTAN_ARGS, "s.1gdb"], cwd=work)
        return run(["ONEview", "-h", "s.1ano"], cwd=work).stdout
    finally:
        shutil.rmtree(work, ignore_errors=True)


def trims_last_unit(tmpdir) -> bool:
    """True if the installed FasTAN ends each interval where the lagging copy of its first-wave alignment ends, one
    unit before the array does: the 0.8 release does (fixed upstream on 2026-08-09). Probe: ten exact copies of a
    150-bp unit at [2000, 3500) between random flanks."""
    r = random.Random(1)
    rs = lambda n: "".join(r.choice("ACGT") for _ in range(n))
    unit = rs(150)
    ends = [e for _, e, u, _, _ in parse_ano(run_fastan(rs(2000) + unit * 10 + rs(2000), tmpdir)) if u == 150]
    if not ends:
        raise CenrepError("FasTAN self-test found no array (10 copies of a 150-bp unit); check FasTAN, FAtoGDB and "
                          "ONEview (FASTGA)")
    return max(ends) < 3500 - 75


def setup(tmpdir) -> bool:
    trims = trims_last_unit(tmpdir)
    if trims:
        log.info("seeds: this FasTAN build ends each array one copy early (0.8 release); ends are extended by one unit")
    return trims


def unit_consensus(units: list, period: int) -> str:
    """Majority unit of an array: abPOA consensus of up to CONS_UNITS copies spread along it. Copies off the period by
    more than RELATED_LEN_TOL are misparsed (a missed or extra cut) and left out. One usable copy is returned as is;
    none gives ''."""
    global _POA
    ok = [u for u in units if abs(len(u) - period) <= RELATED_LEN_TOL * period]
    if len(ok) > CONS_UNITS:
        ok = [ok[i] for i in np.linspace(0, len(ok) - 1, CONS_UNITS).round().astype(int)]
    if len(ok) < 2:
        return ok[0] if ok else ""
    if _POA is None:                                 # most frequent base per column (gap-majority columns dropped):
        _POA = pyabpoa.msa_aligner(cons_algrm="MF")  # the majority rule of families.msa_consensus
    return _POA.msa(ok, out_cons=True, out_msa=False).cons_seq[0]


def chunk_rows(chrom: str, seq: str, offset: int, tmpdir, trims: bool) -> list:
    """Seed-table rows (seeds.COLS) of the FasTAN arrays of one chunk at offset, units up to SEED_MAX_PERIOD (longer
    ones are the long-unit detector's). pct_match = 100 - FasTAN's % differences between adjacent copies; score = the
    matched bp; FasTAN measures no indel share or entropy (NaN). trims: see trims_last_unit."""
    seq = seq.upper()
    if not any(b in seq for b in "ACGT"):
        return []                                    # FAtoGDB rejects a sequence without bases
    out = []
    for beg, end, unit, diff, cuts in parse_ano(run_fastan(seq, tmpdir)):
        if unit > SEED_MAX_PERIOD:
            continue
        if trims:                                    # add the last copy; an array never crosses an assembly gap
            gap = seq.find("N", end, end + unit)
            end = min(len(seq), end + unit) if gap < 0 else gap
        cons = unit_consensus([seq[beg + a:beg + b] for a, b in zip(cuts, cuts[1:])], unit) or seq[beg:beg + unit]
        span = end - beg
        out.append((chrom, offset + beg, offset + end, unit, round(span / unit, 1), len(cons), 100 - diff, NAN,
                    round(span * (100 - diff) / 100), NAN, cons))
    return out
