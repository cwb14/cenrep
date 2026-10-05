import time

import numpy as np
import pandas as pd
import pytest

from cenrep.annotate import ARR_COLS, MON_COLS
from cenrep.genomes import Genome
from cenrep.hor import (assign_hors, hor_candidate, is_hor_unit, lcp_array, maximal_repeats, run_hor, smallest_period,
                        suffix_array, tandem_runs, variant_classes)
from cenrep.params import Params
from conftest import have
from synth import mutate, rseq


def brute_sa(t):
    return sorted(range(len(t)), key=lambda i: list(t[i:]))


def brute_maximal(t, min_len=2):
    n, out = len(t), set()
    for l in range(min_len, n + 1):
        occ = {}
        for i in range(n - l + 1):
            occ.setdefault(tuple(t[i:i + l]), []).append(i)
        for ps in occ.values():
            if len(ps) < 2:
                continue
            left = {t[q - 1] if q > 0 else ("start", q) for q in ps}
            right = {t[q + l] if q + l < n else ("end", q) for q in ps}
            if len(left) > 1 and len(right) > 1:
                out.add((l, tuple(ps)))
    return out


def test_suffix_array_and_lcp_match_brute_force():
    rng = np.random.default_rng(1)
    for _ in range(60):
        t = rng.integers(0, 3, size=int(rng.integers(1, 40)))
        if rng.random() < 0.5:
            t[rng.integers(0, len(t))] = -1
        sa = suffix_array(t)
        assert sa.tolist() == brute_sa(t.tolist())
        lcp = lcp_array(t, sa)
        for i in range(1, len(t)):
            a, b = t[sa[i - 1]:].tolist(), t[sa[i]:].tolist()
            k = 0
            while k < min(len(a), len(b)) and a[k] == b[k]:
                k += 1
            assert lcp[i] == k


def test_maximal_repeats_match_brute_force():
    rng = np.random.default_rng(2)
    for _ in range(60):
        t = rng.integers(0, 3, size=int(rng.integers(2, 30)))
        sa = suffix_array(t)
        got = {(l, tuple(sorted(int(x) for x in occ)))
               for l, occ in maximal_repeats(t, sa, lcp_array(t, sa), max_len=10 ** 6)}
        assert got == brute_maximal(t.tolist())


def test_smallest_period_and_unit_rules():
    assert smallest_period((1, 2, 3, 1, 2, 3)) == 3 and smallest_period((1, 2, 1)) == 2
    assert is_hor_unit((1, 2, 3, 4)) and is_hor_unit((1, 2, 1)) and is_hor_unit((1, 2, 3, 4, 1, 2, 3))
    assert not is_hor_unit((1, 2, 1, 2)) and not is_hor_unit((5, 5, 5)) and not is_hor_unit((7,))


def _hors(t, max_len=200):
    t = np.array(t)
    sa = suffix_array(t)
    return assign_hors(t, maximal_repeats(t, sa, lcp_array(t, sa), max_len, keep=hor_candidate))


def test_tandem_hor_prefers_full_unit_over_partial_repeat():
    units = _hors([-1] + [0, 1, 2, 3] * 14 + [0, 1, 2] + [-2])
    assert len(units) == 1
    u, acc = units[0]
    assert u == (0, 1, 2, 3) and len(acc) == 14 and tandem_runs(acc, 4) == [acc]


def test_exact_tandem_hor_without_partial_copy():
    units = _hors([-1] + [0, 1, 2, 3] * 14 + [-2])
    assert units == [((0, 1, 2, 3), list(range(1, 57, 4)))]


def test_non_tandem_duplication():
    assert _hors([-1, 9, 0, 1, 10, 11, 0, 1, 12, -2]) == [((0, 1), [2, 6])]


def test_homogeneous_run_is_not_a_hor_and_is_fast():
    t0 = time.time()
    assert _hors([-1] + [0] * 20000 + [-2]) == []
    assert time.time() - t0 < 10


@pytest.mark.skipif(not have("cd-hit-est"), reason="cd-hit-est missing")
def test_variant_classes_groups_near_identical(tmp_path):
    rng = np.random.default_rng(4)
    base = rseq(rng, 150)
    v = [mutate(rng, base, 0.08) for _ in range(3)]
    cl = variant_classes({f"m{i}": mutate(rng, v[i % 3], 0.002) for i in range(12)}, tmp_path, 2, 0.97)
    groups = {}
    for name, k in cl.items():
        groups.setdefault(k, set()).add(int(name[1:]) % 3)
    assert len(groups) == 3 and all(len(g) == 1 for g in groups.values())


@pytest.mark.skipif(not have("cd-hit-est"), reason="cd-hit-est missing")
def test_run_hor_tiny_region(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">c\n" + "ACGT" * 250 + "\n")
    import pysam
    pysam.faidx(str(fa))
    g = Genome("g", fa, {"c": 1000})
    mon = pd.DataFrame([(0, "c", 0, 150, "+", "f", 90.0, 1e-20, 1, 150, 150, 0.9, 0)],
                       columns=["mon_idx"] + MON_COLS + ["array_id"])
    arrays = pd.DataFrame([(0, "c", 0, 150, "+", "f", 1, 0.9, 150)], columns=ARR_COLS)
    regions = {"c": dict(evidence="satellite", cen_start=0, cen_end=150)}
    mh, occ, units = run_hor(g, mon, arrays, regions, {"f"}, tmp_path / "w", 1, Params())
    assert occ.empty and units.empty and len(mh) == 1
