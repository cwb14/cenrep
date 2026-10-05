import pandas as pd
import pysam

from cenrep.annotate import ARR_COLS
from cenrep.centromere import finalize, load_regions, longunit_call, satellite_call, write_regions
from cenrep.longunit import LU_COLS
from cenrep.params import Params


def test_satellite_call_main_cluster_and_competing():
    a = pd.DataFrame([(0, "c", 5_000_000, 5_200_000, "+", "A", 1100, 0.9, 180),
                      (1, "c", 5_300_000, 5_350_000, "+", "A", 270, 0.9, 180),
                      (2, "c", 8_000_000, 8_120_000, "+", "A", 660, 0.9, 180)], columns=ARR_COLS)
    call = satellite_call(a, {"A"}, 10_000_000, Params())
    assert (call["cen_start"], call["cen_end"]) == (5_000_000, 5_350_000)
    assert call["sat_bp"] == 250_000 and call["second_cluster_bp"] == 120_000 and call["evidence"] == "satellite"


def test_satellite_call_weak_and_absent():
    p = Params()
    small = pd.DataFrame([(0, "c", 100_000, 108_000, "+", "A", 44, 0.9, 180)], columns=ARR_COLS)
    assert satellite_call(small, {"A"}, 10_000_000, p)["evidence"] == "weak_satellite"
    assert satellite_call(small, {"B"}, 10_000_000, p) is None


def test_finalize_flags_gaps_and_near_end(tmp_path):
    seq = "A" * 1000 + "N" * 20 + "A" * 1000 + "N" * 5 + "A" * 2000
    fa = tmp_path / "x.fa"
    fa.write_text(">c\n" + seq + "\n")
    pysam.faidx(str(fa))
    p = Params(subtel_min_bp=500)
    call = dict(cen_start=100, cen_end=3000, evidence="satellite", main_family="A", main_family_bp=2000,
                sat_bp=2000, second_cluster_bp=0, cen_families="A")
    with pysam.FastaFile(str(fa)) as h:
        row = finalize("c", len(seq), call, h, None, None, p)
    assert row["n_gaps"] == 1 and row["gap_bp"] == 20
    assert row["flags"] == "low_satellite,near_end,gaps:1" and row["confidence"] == "low"


def test_none_row_when_no_call():
    row = finalize("c", 1000, None, None, None, None, Params())
    assert row["evidence"] == "none" and row["confidence"] == "none" and row["cen_start"] == "."


def test_longunit_call_requires_concentration():
    p = Params()
    lu = pd.DataFrame([("c", 4_000_000, 4_600_000, 2700, 222.2, 0.95)], columns=LU_COLS)
    call = longunit_call(lu, 10_000_000, None, None, p)
    assert call["evidence"] == "long_unit_tandem" and call["cen_start"] == 4_000_000
    scattered = pd.DataFrame([("c", s, s + 30_000, 2700, 11.1, 0.9) for s in range(1_000_000, 9_000_000, 1_000_000)],
                             columns=LU_COLS)
    assert longunit_call(scattered, 10_000_000, None, None, p) is None


def test_regions_roundtrip(tmp_path):
    rows = {"c1": finalize("c1", 1000, None, None, None, None, Params())}
    rows["c2"] = dict(rows["c1"], chrom="c2", cen_start=10, cen_end=500, evidence="satellite", confidence="high")
    write_regions(rows, tmp_path / "r.tsv")
    back = load_regions(tmp_path / "r.tsv")
    assert back["c1"]["evidence"] == "none" and back["c2"]["cen_end"] == 500 and back["c2"]["chrom_len"] == 1000


def test_diffuse_region_is_flagged_and_not_high(tmp_path):
    """Real data (Coffea): arrays scattered through pericentromeric heterochromatin chain into 5-22 Mb regions
    that are only 6-9% satellite; such a region must not be reported as a high-confidence centromere."""
    fa = tmp_path / "x.fa"
    fa.write_text(">c\n" + "A" * 20_000_000 + "\n")
    pysam.faidx(str(fa))
    call = dict(cen_start=5_000_000, cen_end=12_000_000, evidence="satellite", main_family="A",
                main_family_bp=700_000, sat_bp=700_000, second_cluster_bp=0, cen_families="A")
    with pysam.FastaFile(str(fa)) as h:
        row = finalize("c", 20_000_000, call, h, None, None, Params())
    assert "diffuse" in row["flags"] and row["confidence"] == "medium"


def _contigs(tmp_path):
    from cenrep.genomes import Genome
    fa = tmp_path / "ctg.fa"
    fa.write_text("".join(f">{c}\n{'A' * 200_000}\n" for c in ("c1", "c2", "c3")))
    pysam.faidx(str(fa))
    arrays = pd.DataFrame([(0, "c1", 0, 60_000, "+", "A", 333, 0.9, 180),
                           (1, "c2", 50_000, 58_000, "+", "A", 44, 0.9, 180)], columns=ARR_COLS)
    return Genome("ctg", fa, {"c1": 200_000, "c2": 200_000, "c3": 200_000}), arrays


def test_borrowed_calls_need_satellite_evidence(tmp_path, monkeypatch):
    """On a borrowed (contig) genome only full satellite calls count: a contig without the centromeric
    satellite is no evidence of a satellite-free centromere, so no weak-satellite or long-unit fallback, and
    the long-unit scan covers call flanks only."""
    from cenrep import centromere, longunit
    from cenrep.trf import COLS
    g, arrays = _contigs(tmp_path)
    seen = []
    lu = pd.DataFrame([("c3", 80_000, 140_000, 3000, 20.0, 0.95)], columns=LU_COLS)
    monkeypatch.setattr(longunit, "detect", lambda fa, targets, *a: seen.append(targets) or lu)
    p = Params(subtel_min_bp=10_000)
    trf = pd.DataFrame(columns=COLS)
    rows, _ = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t1", 1, p)
    assert [rows[c]["evidence"] for c in ("c1", "c2", "c3")] == ["satellite", "weak_satellite", "long_unit_tandem"]
    rows, _ = centromere.call_regions(g, arrays, {"A"}, None, trf, tmp_path / "t2", 1, p, borrowed=True)
    assert [rows[c]["evidence"] for c in ("c1", "c2", "c3")] == ["satellite", "none", "none"]
    assert seen[-1] == [("c1", 0, 200_000)]


def test_borrowed_contig_end_flag_does_not_lower_confidence(tmp_path):
    g, _ = _contigs(tmp_path)
    call = dict(cen_start=0, cen_end=60_000, evidence="satellite", main_family="A", main_family_bp=60_000,
                sat_bp=60_000, second_cluster_bp=0, cen_families="A")
    with pysam.FastaFile(str(g.fasta)) as h:
        own = finalize("c1", 200_000, call, h, None, None, Params())
        bor = finalize("c1", 200_000, call, h, None, None, Params(), borrowed=True)
    assert own["flags"] == "near_end" and own["confidence"] == "medium"
    assert bor["flags"] == "contig_end" and bor["confidence"] == "high"
