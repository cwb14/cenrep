import numpy as np
import pandas as pd
import pysam

from cenrep.annotate import ARR_COLS, MON_COLS
from cenrep.centromere import finalize
from cenrep.classify import CLS_COLS
from cenrep.genomes import Genome
from cenrep.hor import MONHOR_COLS, OCC_COLS, UNIT_COLS
from cenrep.longunit import LU_COLS
from cenrep.output import attrs, final_names, tandem_loci, write_genome
from cenrep.params import Params
from cenrep.seeds import COLS as SEED_COLS
from gffcheck import check_gff


def test_attrs_escaping():
    assert attrs({"ID": "a", "flags": "gaps:1,x=y;z", "none": None}) == "ID=a;flags=gaps:1,x%3Dy%3Bz"


def test_tandem_loci_merge_and_top_row():
    rows = pd.DataFrame([("c", 0, 100, 2, 50.0, 2, 90, 0, 50, 1.0, "AT"),
                         ("c", 50, 300, 4, 60.0, 4, 90, 0, 80, 1.5, "ATGC"),
                         ("c", 400, 500, 10, 10.0, 10, 90, 0, 40, 1.9, "ACGTACGTAA")], columns=SEED_COLS)
    loci = tandem_loci(rows)
    assert [(x["start"], x["end"], x["period"]) for x in loci] == [(0, 300, 4), (400, 500, 10)]


def test_final_names_rank_by_bp():
    fam = pd.DataFrame({"fam_id": ["fam001", "fam002"], "consensus_len": [178, 90]})
    arrays = {"g": pd.DataFrame([(0, "c", 0, 1000, "+", "fam002", 11, .9, 90),
                                 (1, "c", 0, 500, "+", "fam001", 3, .9, 178)], columns=ARR_COLS)}
    assert final_names(fam, arrays, "Cx") == {"fam002": "Cx_sat01_90", "fam001": "Cx_sat02_178"}


def _fams(rows):
    base = dict(bp=1, n_arrays=1, prevalence=1.0, concentration=1.0, subtel_frac=0.0, gene_ratio=np.nan,
                median_main_bp=1.0, eligible=True, strength="strong", defines="Sp_A01", reason="r")
    return pd.DataFrame([dict(base, family=f, cls=c) for f, c in rows], columns=CLS_COLS)


def _data(tmp_path, with_region=True):
    fa = tmp_path / "Sp_A.fa"
    fa.write_text(">Sp_A01\n" + "A" * 20000 + "\n>Sp_A02\n" + "C" * 5000 + "\n")
    pysam.faidx(str(fa))
    g = Genome("Sp_A", fa, {"Sp_A01": 20000, "Sp_A02": 5000})
    arrays = pd.DataFrame([(0, "Sp_A01", 5000, 5600, "+", "fam001", 4, 0.9, 150),
                           (1, "Sp_A01", 5600, 6050, "-", "fam001", 3, 0.88, 150),
                           (3, "Sp_A01", 6100, 6400, "+", "fam002", 3, 0.9, 100),
                           (2, "Sp_A01", 15000, 15450, "+", "fam001", 3, 0.8, 150)], columns=ARR_COLS)
    mons = [("Sp_A01", 5000 + 150 * i, 5150 + 150 * i, "+", "fam001", 90.0, 1e-20, 1, 150, 150, 0.9, 0) for i in range(4)]
    mons += [("Sp_A01", 5600 + 150 * i, 5750 + 150 * i, "-", "fam001", 80.0, 1e-18, 1, 150, 150, 0.88, 1) for i in range(3)]
    mons += [("Sp_A01", 6100 + 100 * i, 6200 + 100 * i, "+", "fam002", 50.0, 1e-10, 1, 100, 100, 0.9, 3) for i in range(3)]
    mons += [("Sp_A01", 15000 + 150 * i, 15150 + 150 * i, "+", "fam001", 70.0, 1e-15, 1, 150, 150, 0.8, 2) for i in range(3)]
    mon = pd.DataFrame(mons, columns=MON_COLS + ["array_id"])
    mon.insert(0, "mon_idx", np.arange(len(mon)))
    rows = {"Sp_A02": finalize("Sp_A02", 5000, None, None, None, None, Params())}
    if with_region:
        rows["Sp_A01"] = dict(finalize("Sp_A01", 20000, None, None, None, None, Params()), cen_start=5000,
                              cen_end=7500, evidence="satellite", confidence="high", flags="gaps:1",
                              main_family="fam001", main_family_bp=1050, n_gaps=1, gap_bp=20)
    else:
        rows["Sp_A01"] = finalize("Sp_A01", 20000, None, None, None, None, Params())
    seeds = pd.DataFrame([("Sp_A01", 5000, 5600, 150, 4.0, 150, 90, 1, 500, 1.9, "A" * 150),
                          ("Sp_A01", 6500, 7100, 2, 300.0, 2, 95, 0, 900, 1.0, "AT")], columns=SEED_COLS)
    mh = pd.DataFrame([(i, i % 2, "HOR_Sp_A01_1" if i in (0, 1) else ".", "hor_Sp_A01_1" if i in (0, 1) else ".")
                       for i in range(7)], columns=MONHOR_COLS)
    occ = pd.DataFrame([("hor_Sp_A01_1", "HOR_Sp_A01_1", "Sp_A01", 5000, 5300, "+", 0, 2, "HOR_Sp_A01_1.r1", 1)],
                       columns=OCC_COLS)
    units = pd.DataFrame([("HOR_Sp_A01_1", "fam001", "Sp_A01", 2, 2, 1, "0-1")], columns=UNIT_COLS)
    return dict(genome=g, arrays=arrays, monomers=mon, fams=_fams([("fam001", "centromeric"), ("fam002", "cen_associated")]),
                regions=rows, longunit=pd.DataFrame(columns=LU_COLS), mon_hor=mh, hor_occ=occ, hor_units=units,
                seeds=seeds, genes=None)


def test_write_genome_produces_valid_gff(tmp_path):
    d = _data(tmp_path)
    names = {"fam001": "Sp_sat01_150", "fam002": "Sp_sat02_100"}
    write_genome(d, names, tmp_path, Params())
    feats = check_gff(tmp_path / "Sp_A_centromeric_repeat.gff", d["genome"].chroms)
    types = pd.Series([f[3] for f in feats]).value_counts().to_dict()
    assert types == {"satellite_monomer": 13, "satellite_DNA": 4, "centromere": 1, "tandem_repeat": 1,
                     "higher_order_repeat": 1}
    outside = [f for f in feats if f[3] == "satellite_DNA" and f[1] == 15001][0]
    assert outside[4]["in_centromere"] == "false" and "Parent" not in outside[4]
    tr = [f for f in feats if f[3] == "tandem_repeat"][0]
    assert tr[4]["period"] == "2" and tr[4]["Parent"] == "cen_Sp_A01" and tr[4]["unit_source"] == "fastan"
    s = pd.read_csv(tmp_path / "Sp_A_centromere_summary.tsv", sep="\t", dtype=str)
    assert s.set_index("chrom").loc["Sp_A02", "confidence"] == "none"
    assert s.set_index("chrom").loc["Sp_A01", "cen_start"] == "5001"


def test_write_genome_without_regions(tmp_path):
    d = _data(tmp_path, with_region=False)
    d["fams"] = _fams([("fam001", "other"), ("fam002", "other")])
    write_genome(d, {"fam001": "Sp_sat01_150", "fam002": "Sp_sat02_100"}, tmp_path, Params())
    assert check_gff(tmp_path / "Sp_A_centromeric_repeat.gff", d["genome"].chroms) == []
