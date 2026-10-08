import numpy as np
import pandas as pd

from cenrep.annotate import ARR_COLS
from cenrep.classify import assign_cen_associated, choose_centromeric, family_metrics, make_clusters
from cenrep.params import Params

P = Params()
CH = {f"c{i}": 10_000_000 for i in range(1, 5)}


def arr(fam, chrom, start, n, mlen=180):
    return (0, chrom, start, start + n * mlen, "+", fam, n, 0.9, mlen)


def test_make_clusters_bp_and_order():
    cl = make_clusters([(0, 100), (50, 200), (10_000, 10_050), (600_000, 600_300)], 500)
    assert [(c.start, c.end, c.bp) for c in cl] == [(600_000, 600_300, 300), (0, 200, 200), (10_000, 10_050, 50)]


def test_chromosome_specific_and_layered_families():
    rows = [arr("A", "c1", 5_000_000, 1000), arr("A", "c2", 5_000_000, 1000),
            arr("B", "c3", 4_000_000, 800), arr("B", "c4", 4_000_000, 800),
            arr("C", "c1", 5_300_000, 200),                                      # layered next to A on c1
            arr("D", "c1", 100_000, 300), arr("D", "c1", 9_700_000, 300),
            arr("D", "c3", 50_000, 300), arr("D", "c3", 9_800_000, 300)]         # subtelomeric: both ends
    rows += [arr("S", f"c{1 + i % 4}", 300_000 + (i // 4) * 900_000, 120) for i in range(40)]   # scattered
    info = {"A": dict(superfamily="A", length=180, rdna="."), "B": dict(superfamily="B", length=180, rdna="."),
            "C": dict(superfamily="A", length=175, rdna="."), "D": dict(superfamily="D", length=90, rdna="."),
            "S": dict(superfamily="S", length=120, rdna=".")}
    fams, mains = family_metrics(pd.DataFrame(rows, columns=ARR_COLS), CH, None, P, info)
    fams, per_chrom = choose_centromeric(fams, mains, CH, P, info)
    f = fams.set_index("family")
    assert f.loc["A", "cls"] == f.loc["B", "cls"] == f.loc["C", "cls"] == "centromeric"
    assert f.loc["D", "cls"] == "other" and "subtel_frac" in f.loc["D", "reason"]
    assert f.loc["S", "cls"] == "other" and "concentration" in f.loc["S", "reason"]
    assert per_chrom["c1"] == ["A", "C"] and per_chrom["c3"] == ["B"]
    assert f.loc["A", "strength"] == "strong" and f.loc["C", "strength"] == "weak"


def test_gene_rich_cluster_is_vetoed():
    genes = {c: np.arange(0, 10_000_000, 5_000, dtype=np.int64) for c in ("c1", "c2")}
    fams, _ = family_metrics(pd.DataFrame([arr("G", "c1", 5_000_000, 500)], columns=ARR_COLS), CH, genes, P)
    assert not fams.eligible.iloc[0] and "gene_ratio" in fams.reason.iloc[0]


def test_empty_arrays():
    fams, mains = family_metrics(pd.DataFrame(columns=ARR_COLS), CH, None, P)
    fams, per_chrom = choose_centromeric(fams, mains, CH, P)
    assert fams.empty and per_chrom == {}


def test_assign_cen_associated():
    rows = [arr("A", "c1", 5_000_000, 1000), arr("X", "c1", 5_050_000, 100), arr("X", "c2", 1_500_000, 20)]
    arrays = pd.DataFrame(rows, columns=ARR_COLS)
    fams, mains = family_metrics(arrays, CH, None, P)
    fams, _ = choose_centromeric(fams, mains, CH, P)
    regions = {"c1": dict(evidence="satellite", cen_start=5_000_000, cen_end=5_180_000),
               "c2": dict(evidence="none")}
    out = assign_cen_associated(fams, arrays, regions, P).set_index("family")
    assert out.loc["X", "cls"] == "cen_associated"                            # 18 of 21.6 kb inside


def test_unrelated_adjacent_family_is_not_co_centromeric():
    """Review C2: an eligible 5S-rDNA-like 500-bp family next to the CEN180 cluster must not become primary;
    a related layer (same superfamily, or period within 15%) still does."""
    rows = [arr("A", "c1", 5_000_000, 1000), arr("A", "c2", 5_000_000, 1000),
            arr("R", "c1", 5_200_000, 60, 500),                       # abutting, unrelated 500-bp family
            arr("L", "c2", 5_250_000, 200, 159)]                      # 159 vs 180: within 15% (AthCEN159-like)
    info = {"A": dict(superfamily="A", length=180, rdna="."), "R": dict(superfamily="R", length=500, rdna="."),
            "L": dict(superfamily="L", length=159, rdna=".")}
    fams, mains = family_metrics(pd.DataFrame(rows, columns=ARR_COLS), CH, None, P, info)
    fams, per_chrom = choose_centromeric(fams, mains, CH, P, info)
    assert per_chrom["c1"] == ["A"] and per_chrom["c2"] == ["A", "L"]
    assert fams.set_index("family").loc["R", "cls"] == "other"


def test_rdna_family_is_never_eligible():
    rows = [arr("R", "c1", 5_000_000, 200, 500), arr("R", "c2", 5_000_000, 200, 500)]
    info = {"R": dict(superfamily="R", length=500, rdna="5S_rRNA")}
    fams, _ = family_metrics(pd.DataFrame(rows, columns=ARR_COLS), CH, None, P, info)
    assert not fams.eligible.iloc[0] and "rDNA (5S_rRNA)" in fams.reason.iloc[0]


def test_strong_family_takes_precedence_over_larger_unrelated_satellite():
    """Real data (B. juncea A09): the genome's centromeric satellite (strong on most chromosomes) has a 69 kb
    cluster, but a pericentromeric 238-bp satellite forms a larger 162 kb cluster nearby. The strong family's
    superfamily must define the centromere where it has a >= chrom_min_bp cluster; the largest-cluster rule is
    only the fallback (chromosome without such a cluster)."""
    info = {"A": dict(superfamily="A", length=180, rdna="."), "T": dict(superfamily="T", length=238, rdna=".")}
    rows = [arr("A", c, 5_000_000, 1000) for c in ("c1", "c2", "c3")]          # strong: 3/4 chromosomes, 180 kb
    rows += [arr("A", "c4", 5_000_000, 390), arr("T", "c4", 4_600_000, 670, 238)]   # c4: A 70 kb, T 160 kb
    fams, mains = family_metrics(pd.DataFrame(rows, columns=ARR_COLS), CH, None, P, info)
    fams, per_chrom = choose_centromeric(fams, mains, CH, P, info)
    assert per_chrom["c4"] == ["A"]
    rows = rows[:3] + [arr("T", "c4", 4_600_000, 670, 238)]                      # c4 without A: fallback to T
    fams, mains = family_metrics(pd.DataFrame(rows, columns=ARR_COLS), CH, None, P, info)
    fams, per_chrom = choose_centromeric(fams, mains, CH, P, info)
    assert per_chrom["c4"] == ["T"]


def test_widespread_pericentromeric_satellite_gets_no_precedence():
    """Real data (B. napus A06, v0.3.6): a 228-bp pericentromeric satellite with >= 100 kb clusters on half the
    chromosomes passed the 'strong' prevalence test, and where its cluster was larger than CentBr's it defined the
    centromere. Precedence belongs only to families whose kind is also the largest cluster on >= strong_prevalence
    of chromosomes; two such unrelated families (each winning half) still compete by cluster size."""
    info = {"A": dict(superfamily="A", length=180, rdna="."), "T": dict(superfamily="T", length=238, rdna="."),
            "B": dict(superfamily="B", length=120, rdna=".")}
    rows = [arr("A", c, 5_000_000, 1000) for c in ("c1", "c2", "c3")]              # A largest on c1-c3 (180 kb)
    rows += [arr("T", c, 4_500_000, 500, 238) for c in ("c1", "c2")]                # T 119 kb beside A
    rows += [arr("A", "c4", 5_000_000, 390), arr("T", "c4", 4_600_000, 670, 238)]  # c4: A 70 kb, T 160 kb
    fams, mains = family_metrics(pd.DataFrame(rows, columns=ARR_COLS), CH, None, P, info)
    assert fams.set_index("family").loc["T", "prevalence"] >= P.strong_prevalence
    fams, per_chrom = choose_centromeric(fams, mains, CH, P, info)
    assert per_chrom["c4"] == ["A"]
    assert fams.set_index("family").loc["T", "cls"] == "other"
    rows = [arr("A", c, 5_000_000, 1000) for c in ("c1", "c2")] + [arr("B", c, 5_000_000, 1200, 120)
                                                                   for c in ("c3", "c4")]
    rows += [arr("B", "c1", 4_500_000, 900, 120)]                                   # B 108 kb beside A 180 kb
    fams, mains = family_metrics(pd.DataFrame(rows, columns=ARR_COLS), CH, None, P, info)
    fams, per_chrom = choose_centromeric(fams, mains, CH, P, info)
    assert per_chrom["c1"] == ["A"] and per_chrom["c3"] == ["B"]


def test_borrowed_classification_uses_source_families():
    """Contig assemblies: per-chromosome tests are meaningless on contigs (a contig holds only part of a centromere),
    so a target genome takes the union of its sources' centromeric families; its own metrics are reported."""
    from cenrep.classify import borrow_centromeric
    ch = {f"t{i}": 1_000_000 for i in range(1, 4)}
    rows = [arr("A", "t1", 900_000, 500), arr("A", "t2", 0, 500), arr("B", "t3", 100_000, 200),
            arr("D", "t1", 0, 30)]
    fams, mains = family_metrics(pd.DataFrame(rows, columns=ARR_COLS), ch, None, P)
    src = {"S1": pd.DataFrame([dict(family="A", cls="centromeric", strength="weak"),
                               dict(family="B", cls="other", strength="."),
                               dict(family="Z", cls="centromeric", strength="strong")]),
           "S2": pd.DataFrame([dict(family="A", cls="centromeric", strength="strong")])}
    out, per_chrom = borrow_centromeric(fams, mains, src, P)
    f = out.set_index("family")
    assert len(out) == len(fams) and f.index.is_unique and "Z" not in f.index
    assert f.loc["A", "cls"] == "centromeric" and f.loc["A", "strength"] == "strong"
    assert f.loc["A", "defines"] == "t1,t2" and f.loc["A", "reason"].startswith("borrowed from S1+S2")
    assert "own contig metrics: " in f.loc["A", "reason"]
    assert f.loc["B", "cls"] == "other" and f.loc["B", "reason"].startswith("not centromeric in S1+S2")
    assert f.loc["D", "defines"] == "." and per_chrom == {"t1": ["A"], "t2": ["A"]}


def test_rdna_associated_family_becomes_ineligible():
    from cenrep.classify import apply_rdna
    fams = pd.DataFrame(dict(family=["igs", "sat", "few"], eligible=[True, True, True],
                             reason=["eligible"] * 3))
    out = apply_rdna(fams, {"igs": (10, 9), "sat": (10, 0), "few": (2, 1)}, Params()).set_index("family")
    assert not out.loc["igs", "eligible"] and out.loc["igs", "reason"].startswith("rDNA-associated")
    assert out.loc["sat", "eligible"] and out.loc["few", "eligible"] == (1 / 2 < Params().rdna_array_frac)


def test_te_rejected_needs_a_majority_of_disagreeing_chromosomes():
    from cenrep.classify import Cluster, te_rejected
    core = dict(start=1_000_000, end=1_500_000)
    cores = {c: core for c in ("c1", "c2", "c3")}
    at, away = Cluster(1_100_000, 1_400_000, 200_000), Cluster(8_000_000, 8_300_000, 200_000)
    fams = pd.DataFrame([dict(family="F", cls="centromeric", defines="c1,c2,c3"),
                         dict(family="H", cls="other", defines=".")])
    m = lambda x: {("F", "c1"): x[0], ("F", "c2"): x[1], ("F", "c3"): x[2]}
    assert te_rejected(fams, m([at, at, away]), cores) == {}                    # 1 of 3 disagree
    assert te_rejected(fams, m([away, away, at]), cores) == {"F": (2, 3)}
    assert te_rejected(fams, m([away, away, away]), {}) == {}                   # no core: no test


def test_te_vote_is_per_family_kind_so_subfamilies_are_not_rejected_alone():
    """A. thaliana CEN180 subfamily / B. juncea CentBr2: a subfamily defining few chromosomes, away from the core
    there, belongs to a kind that agrees overall; an unrelated family must still be judged on its own."""
    from cenrep.classify import Cluster, te_rejected
    core = dict(start=1_000_000, end=1_500_000)
    cores = {c: core for c in ("c1", "c2", "c3", "c4")}
    at, away = Cluster(1_100_000, 1_400_000, 200_000), Cluster(8_000_000, 8_300_000, 200_000)
    fams = pd.DataFrame([dict(family="A178", cls="centromeric", defines="c1,c2,c3"),
                         dict(family="A176", cls="centromeric", defines="c4"),
                         dict(family="X90", cls="centromeric", defines="c4")])
    info = {"A178": dict(superfamily="S1", length=178), "A176": dict(superfamily="S1", length=176),
            "X90": dict(superfamily="S9", length=90)}
    mains = {("A178", "c1"): at, ("A178", "c2"): at, ("A178", "c3"): at, ("A176", "c4"): away, ("X90", "c4"): away}
    assert te_rejected(fams, mains, cores, info) == {"X90": (1, 1)}



def test_only_decisive_retrotransposon_evidence_rejects():
    """B. juncea A: only pericentromeric clades (Retand, Athila) qualify; their cores sit beside CentBr. Evidence that
    is not the centromere-specific CRM lineage flags disagreement but never rejects the satellite."""
    from cenrep.classify import Cluster, choose_with_te
    fams = pd.DataFrame([dict(family="CentBr", bp=500_000, n_arrays=20, prevalence=1.0, concentration=0.9,
                              subtel_frac=0.0, gene_ratio=0.1, median_main_bp=250_000, eligible=True, cls="other",
                              strength=".", defines=".", dominant=False, reason="eligible")])
    chroms = {"c1": 20_000_000, "c2": 20_000_000}
    mains = {("CentBr", c): Cluster(1_100_000, 1_400_000, 250_000) for c in chroms}
    far = {c: dict(start=9_000_000, end=9_500_000) for c in chroms}
    p = Params()
    out, _, rej = choose_with_te(fams, mains, chroms, p, None, ("Retand", far, False))
    assert rej == {} and set(out.cls) == {"centromeric"}
    out, _, rej = choose_with_te(fams, mains, chroms, p, None, ("CRM", far, True))
    assert rej == {"CentBr": (2, 2)} and set(out.cls) == {"other"}


def test_end_zone_counts_only_two_ended_arrays():
    """Subtelomeric repeats cap both ends; a centromere is one locus even next to an end. Fragaria x ananassa Chr1D:
    its CENH3-bound Type II array lies 1.0 Mb from the end of a 27-Mb chromosome and was vetoed as subtelomeric by
    the old one-ended zone with a 1-Mb floor."""
    from cenrep.classify import end_zone, end_zone_bp
    L = 27_106_190
    assert end_zone(L, P) == int(0.05 * L) and end_zone(2_000_000, P) == 100_000      # no floor: small chromosomes
    acro = [(25_974_330, 26_073_822)]                                                  # one end only
    assert end_zone_bp(acro, L, P) == 0
    both = [(0, 40_000), (L - 30_000, L)]
    assert end_zone_bp(both, L, P) == 70_000


def test_acrocentric_satellite_can_be_centromeric():
    rows = [arr("A", "c1", 5_000_000, 1000), arr("A", "c2", 5_000_000, 1000), arr("A", "c3", 5_000_000, 1000),
            arr("T", "c4", 9_700_000, 400)]                                            # c4: centromere near one end
    info = {"A": dict(superfamily="A", length=180, rdna="."), "T": dict(superfamily="T", length=146, rdna=".")}
    fams, mains = family_metrics(pd.DataFrame(rows, columns=ARR_COLS), CH, None, P, info)
    fams, per_chrom = choose_centromeric(fams, mains, CH, P, info)
    assert fams.set_index("family").loc["T", "eligible"] and per_chrom["c4"] == ["T"]
