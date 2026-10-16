import pytest

from cenrep.trf import parse_ngs
from cenrep.util import union_bp
from conftest import have

NGS = ("@c1\n"
       "11 60 10 5.0 10 95 0 90 25 25 25 25 2.0 ACGTACGTAC ACGTACGTACACGT FLANKLEFT FLANKRIGHT\n"
       "garbage line\n")


def test_parse_ngs_offsets_and_skips_bad_lines():
    assert parse_ngs(NGS, "c1", 1000) == [("c1", 1010, 1060, 10, 5.0, 10, 95, 0, 90, 2.0, "ACGTACGTAC")]


@pytest.mark.skipif(not have("trf"), reason="trf missing")
def test_trf_finds_planted_satellite_across_chunk_boundary(synth_trf, synth_complex):
    cdir, truth = synth_complex
    df = synth_trf["Pa_A"]
    s, e = truth["genomes"]["Pa_A"]["Pa_A01"]["cen"]      # 130,000-169,000 crosses the 145-150 kb chunk overlap
    k = (df.period / 150).round()                          # TRF may report the HOR block at k x 150 (e.g. 600)
    mult = (k >= 1) & ((df.period / 150 - k).abs() <= 0.03 * k)
    hit = df[(df.chrom == "Pa_A01") & mult & (df.start < e) & (df.end > s)]
    assert union_bp(zip(hit.start, hit.end)) > 0.7 * (e - s)
    for p in hit.itertuples():                            # no chunk duplicates left
        dup = hit[(hit.index != p.Index) & ((hit.period - p.period).abs() <= 1)
                  & (hit.start < p.end) & (hit.end > p.start)]
        for q in dup.itertuples():
            ov = min(p.end, q.end) - max(p.start, q.start)
            assert ov < 0.5 * min(p.end - p.start, q.end - q.start)
