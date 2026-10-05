import pytest

from cenrep.util import (CenrepError, canonical_rotation, clipped_bp, merge_intervals, parse_fasta_text,
                         revcomp, run, union_bp)


def test_canonical_rotation_is_rotation_invariant():
    s = "GATTACAGATTTCC"
    rots = {canonical_rotation(s[i:] + s[:i]) for i in range(len(s))}
    assert len(rots) == 1 and len(next(iter(rots))) == len(s)


def test_merge_intervals_gap():
    assert merge_intervals([(10, 20), (0, 5), (22, 30)], gap=2) == [[0, 5], [10, 30]]
    assert merge_intervals([(0, 5), (5, 8)]) == [[0, 8]]


def test_union_and_clipped_bp():
    iv = [(0, 10), (5, 15), (20, 30)]
    assert union_bp(iv) == 25
    assert clipped_bp(iv, 8, 25) == 12


def test_revcomp_and_fasta_parse():
    assert revcomp("AACGTn") == "nACGTT"
    assert parse_fasta_text(">a x\nAC\nGT\n>b\nTT\n") == [("a", "ACGT"), ("b", "TT")]


def test_run_raises_on_failure():
    with pytest.raises(CenrepError):
        run(["false"])
