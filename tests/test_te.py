import math

import numpy as np
import pytest

from cenrep.te import agrees, binom_sf, centromeric_clade, core, load_te_gff, marker

L = 10_000_000


def test_binom_sf_exact_values():
    assert binom_sf(10, 10, 0.5) == pytest.approx(1 / 1024)
    assert binom_sf(0, 10, 0.3) == 1.0
    assert binom_sf(3, 4, 0.5) == pytest.approx(5 / 16)


def test_load_te_gff_reads_clades_from_ltrquest_and_edta_style_lines(tmp_path):
    gff = tmp_path / "g.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "c1\tLTRquest\tLTR_retrotransposon\t101\t200\t.\t+\t.\tID=a;superfamily=Gypsy;clade=CRM;family=f1\n"
        "c1\tLTRquest\tlong_terminal_repeat\t101\t120\t.\t+\t.\tParent=a\n"
        "c1\tEDTA\tGypsy_LTR_retrotransposon\t1001\t2000\t.\t-\t.\tID=b;Clade=CRM\n"
        "c1\tEDTA\tCopia_LTR_retrotransposon\t5001\t6000\t.\t-\t.\tID=c;Clade=Ale\n"
        "c1\tLTRquest\tLTR_retrotransposon\t7001\t8000\t.\t+\t.\tID=d;clade=unknown\n"
        "c1\tLTRquest\tLTR_retrotransposon\t9001\t9100\t.\t+\t.\tID=e\n"
        "c9\tLTRquest\tLTR_retrotransposon\t1\t100\t.\t+\t.\tID=f;clade=CRM\n")
    t = load_te_gff(gff, {"c1": 20_000})
    assert set(t) == {"CRM", "Ale"} and t["CRM"]["c1"].tolist() == [150, 1500] and t["Ale"]["c1"].tolist() == [5500]


def test_core_finds_a_concentrated_cluster_and_ignores_uniform_spread():
    rng = np.random.default_rng(1)
    clustered = np.sort(np.concatenate([rng.integers(4_000_000, 4_400_000, 60), rng.integers(0, L, 40)]))
    c = core(clustered, L)
    assert c["significant"] and c["start"] <= 4_200_000 <= c["end"] and c["end"] - c["start"] <= 1_000_000
    assert c["k"] >= 55
    assert not core(np.sort(rng.integers(0, L, 100)), L)["significant"]
    assert core(clustered[:4], L) is None                               # too few elements to estimate a peak
    # Sphagnum-like: 8 of 9 elements of a clade in one 0.2-Mb cluster (moss genomes carry few intact LTR-RTs)
    sparse = np.sort(np.concatenate([rng.integers(6_000_000, 6_200_000, 8), [1_000_000]]))
    assert core(sparse, L)["significant"]


def test_centromeric_clade_is_the_concentrated_one_on_most_chromosomes():
    rng = np.random.default_rng(2)
    chroms = {f"c{i}": L for i in range(5)}
    centro = {c: np.sort(np.concatenate([rng.integers(i * 1_000_000 + 2_000_000, i * 1_000_000 + 2_300_000, 50),
                                         rng.integers(0, L, 30)])) for i, c in enumerate(chroms)}
    spread = {c: np.sort(rng.integers(0, L, 200)) for c in chroms}
    clade, cores = centromeric_clade({"CRM": centro, "Ale": spread}, chroms)
    assert clade == "CRM" and all(cores[c]["significant"] for c in chroms)
    assert centromeric_clade({"Ale": spread}, chroms) == (None, {})
    two = {c: v for c, v in list(centro.items())[:2]}                     # concentrated on only 2 of 5 chromosomes
    assert centromeric_clade({"CRM": two}, chroms)[0] is None


def test_agrees_within_the_wider_of_call_and_core():
    c = dict(start=1_000_000, end=1_500_000)                             # 0.5-Mb core
    assert agrees(1_200_000, 1_300_000, c)                               # inside
    assert agrees(1_900_000, 2_000_000, c)                               # 0.4 Mb away: within the core width
    assert agrees(3_000_000, 6_000_000, c)                               # 1.5 Mb away: within the 3-Mb call width
    assert not agrees(2_100_000, 2_200_000, c)                           # 0.6 Mb away, both narrower


def test_marker_prefers_the_centromere_specific_crm_lineage():
    """CRM (CR chromoviruses: CRM maize, CRR rice, CRB Brassica, CRS sugarcane) is centromere-specific in plants;
    other clades are often pericentromeric. A qualifying CRM is the marker even if another clade is more concentrated,
    and only CRM evidence is decisive."""
    a, b = {"c1": dict(start=1, end=2)}, {"c1": dict(start=5, end=6)}
    assert marker([("Ale", a), ("CRM", b)]) == ("CRM", b, True)
    assert marker([("Retand", a), ("Athila", b)]) == ("Retand", a, False)
    assert marker([]) is None


def test_ltr_overlap_reports_tandem_repeats_inside_retrotransposons(tmp_path):
    """Musa: a 62-bp repeat lies inside SIRE elements; a true satellite array lies outside them."""
    import pandas as pd
    from cenrep.te import load_te_intervals, ltr_overlap
    gff = tmp_path / "g.gff3"
    gff.write_text("##gff-version 3\n"
                   "c1\tx\tLTR_retrotransposon\t1001\t9000\t.\t+\t.\tID=a;clade=SIRE\n"
                   "c1\tx\tLTR_retrotransposon\t20001\t28000\t.\t+\t.\tID=b;clade=SIRE\n"
                   "c1\tx\tLTR_retrotransposon\t24001\t26000\t.\t+\t.\tID=c\n")       # nested, no clade
    iv = load_te_intervals(gff, {"c1": 100_000})
    assert iv["c1"][2] == ["SIRE", "SIRE", "unknown"]
    arrays = pd.DataFrame(dict(chrom=["c1"] * 3, start=[2000, 21000, 50000], end=[2500, 21500, 60000],
                               family=["T", "T", "S"]))
    ov = ltr_overlap(arrays, iv, longest=8000)
    assert ov["T"] == (1.0, "SIRE") and ov["S"] == (0.0, ".")
