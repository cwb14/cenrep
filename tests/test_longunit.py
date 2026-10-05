import numpy as np
import pandas as pd
import pysam
import pytest

from cenrep.longunit import build_masks, detect, mask_seq, merge_hits, windows
from conftest import have
from synth import mutate, rseq


def test_windows():
    assert list(windows(0, 250, 100, 20)) == [(0, 100), (80, 180), (160, 250)]


def test_mask_seq_clips_to_window():
    starts, ends = np.array([5, 20]), np.array([8, 30])
    assert mask_seq("A" * 25, 3, starts, ends) == "AA" + "NNN" + "A" * 12 + "N" * 8


def test_build_masks_merges():
    df = pd.DataFrame({"chrom": ["c", "c", "d"], "start": [0, 5, 1], "end": [10, 20, 2]})
    m = build_masks([df])
    assert m["c"][0].tolist() == [0] and m["c"][1].tolist() == [20] and m["d"][1].tolist() == [2]


def test_merge_hits_multiples_and_distinct_arrays():
    hits = [(1000, 16000, 3000, 0.97), (1000, 16000, 6000, 0.96), (1100, 15800, 9000, 0.95),
            (15500, 30000, 4100, 0.90)]
    out = merge_hits(hits, "c")
    assert sorted((r["unit"], r["start"], r["end"]) for r in out) == [(3000, 1000, 16000), (4100, 15500, 30000)]


@pytest.mark.skipif(not have("minimap2"), reason="minimap2 missing")
def test_detect_finds_3kb_unit(tmp_path, test_params):
    rng = np.random.default_rng(9)
    u = rseq(rng, 3000)
    seq = rseq(rng, 30000) + "".join(mutate(rng, u, 0.03) for _ in range(6)) + rseq(rng, 30000)
    fa = tmp_path / "x.fa"
    fa.write_text(">c\n" + seq + "\n")
    pysam.faidx(str(fa))
    lu = detect(fa, [("c", 0, len(seq))], {}, tmp_path, 2, test_params)
    assert len(lu) == 1 and abs(lu.unit.iloc[0] - 3000) <= 150
    assert abs(lu.start.iloc[0] - 30000) <= 300 and abs(lu.end.iloc[0] - 48000) <= 300
