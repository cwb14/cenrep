import numpy as np
import pytest

from cenrep import fastan
from conftest import have

need_fastan = pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview"), reason="FasTAN tools missing")

ANO = """g
S 1 c
C 5500

M 0 1999 3352
L 3 150
X 0
P 10 70 220 370 520 670 820 970 1120 1270 1420
M 0 4000 5700
L 3 848
X 2
"""


def rand(rng, n):
    return "".join(rng.choice(list("ACGT"), n))


def indel_copy(rng, s, rate):
    """A copy of s with substitutions, deletions and insertions (70/15/15) at total rate."""
    out = []
    for ch in s:
        x = rng.random()
        if x < 0.7 * rate:
            out.append(rng.choice([b for b in "ACGT" if b != ch]))
        elif x < 0.85 * rate:
            continue
        elif x < rate:
            out.append(ch + rng.choice(list("ACGT")))
        else:
            out.append(ch)
    return "".join(out)


def edit_in(a, t):
    """Edit distance of all of a against the best-matching substring of t (test-only reference)."""
    prev = [0] * (len(t) + 1)
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(t)
        for j, ct in enumerate(t, 1):
            cur[j] = min(prev[j - 1] + (ca != ct), prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    return min(prev)


def test_parse_ano_reads_intervals_units_divergence_and_cuts():
    assert fastan.parse_ano(ANO) == [
        (1999, 3352, 150, 0, [70, 220, 370, 520, 670, 820, 970, 1120, 1270, 1420]),
        (4000, 5700, 848, 2, []),                 # a near tandem: no partition line, no cuts of its own
    ]


def test_unparsable_fastan_output_fails_clearly():
    from cenrep.util import CenrepError
    with pytest.raises(CenrepError, match="M 0 1999"):
        fastan.parse_ano("M 0 1999\nL 3 150\n")


def test_consensus_recovers_the_unit_from_divergent_copies():
    rng = np.random.default_rng(3)
    unit = rand(rng, 171)
    copies = [indel_copy(rng, unit, 0.10) for _ in range(40)]
    assert min(edit_in(c, unit) for c in copies[:5]) >= 5          # no single copy is the unit
    assert edit_in(fastan.unit_consensus(copies, 171), unit) <= 1


def test_consensus_leaves_out_misparsed_copies():
    """A missed cut gives a double copy; it must not stretch the consensus to two units."""
    rng = np.random.default_rng(4)
    unit = rand(rng, 120)
    good = [indel_copy(rng, unit, 0.05) for _ in range(10)]
    double = [indel_copy(rng, unit, 0.05) + indel_copy(rng, unit, 0.05) for _ in range(25)]
    c = fastan.unit_consensus(good + double, 120)
    assert abs(len(c) - 120) <= 3 and edit_in(c, unit) <= 2


def test_consensus_of_one_copy_is_that_copy_and_of_none_is_empty():
    assert fastan.unit_consensus(["ACGTTGCA", "A" * 30], 8) == "ACGTTGCA"
    assert fastan.unit_consensus(["A" * 30], 8) == ""


@need_fastan
def test_planted_array_has_its_coordinates_period_and_unit(tmp_path):
    """Whatever the FasTAN build (the 0.8 release ends intervals one unit early), the row spans the planted array in
    chromosome coordinates and its consensus is the planted unit."""
    rng = np.random.default_rng(5)
    unit = rand(rng, 171)
    seq = rand(rng, 3000) + "".join(indel_copy(rng, unit, 0.05) for _ in range(30)) + rand(rng, 3000)
    end = len(seq) - 3000
    rows = fastan.chunk_rows("c1", seq, 10_000, tmp_path, fastan.trims_last_unit(tmp_path))
    hit = [r for r in rows if abs(r[3] - 171) <= 5]
    assert len(hit) == 1
    chrom, s, e, period, copies, csize, pct, _, _, _, cons = hit[0]
    assert chrom == "c1" and abs(s - 13_000) <= 20 and abs(e - (10_000 + end)) <= 20
    assert abs(copies - (e - s) / period) < 0.1 and csize == len(cons) and pct >= 85
    assert edit_in(cons, unit + unit) <= 3


@need_fastan
def test_added_last_copy_never_crosses_an_assembly_gap(tmp_path):
    """The last copies (90 bp) are shorter than the mean unit (~98 bp), so adding one mean unit to the reported end
    overshoots the array; the array ends where the gap starts."""
    rng = np.random.default_rng(6)
    unit = rand(rng, 90)
    left = rand(rng, 2000) + "".join(unit[:45] + rand(rng, 10) + unit[45:] for _ in range(15)) + unit * 5
    seq = left + "N" * 100 + rand(rng, 2000)
    hit = [r for r in fastan.chunk_rows("c1", seq, 0, tmp_path, True) if 85 <= r[3] <= 105]
    assert hit and max(r[2] for r in hit) == len(left)


@need_fastan
def test_units_beyond_the_seed_range_are_left_to_the_long_unit_detector(tmp_path):
    rng = np.random.default_rng(7)
    unit = rand(rng, 2500)
    seq = rand(rng, 3000) + "".join(indel_copy(rng, unit, 0.02) for _ in range(4)) + rand(rng, 3000)
    raw = fastan.parse_ano(fastan.run_fastan(seq, tmp_path))
    assert any(abs(u - 2500) <= 50 for _, _, u, _, _ in raw)                  # FasTAN reports it ...
    assert all(r[3] <= 2000 for r in fastan.chunk_rows("c1", seq, 0, tmp_path, False))   # ... cenrep drops it


@need_fastan
def test_sequence_without_bases_gives_no_rows(tmp_path):
    assert fastan.chunk_rows("c1", "N" * 5000, 0, tmp_path, False) == []
    assert fastan.chunk_rows("c1", "", 0, tmp_path, False) == []
