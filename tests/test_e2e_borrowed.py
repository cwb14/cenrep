"""End to end with a contig assembly that borrows its centromeric families (classify_from).

Ctg_A = Pa_A cut into contigs: chromosomes 1-3 at 150 kb, inside the centromere, so each half-centromere lies in
a contig end zone and the genome's own per-chromosome tests would reject the centromeric family; chromosome 4
(no centromere) in three 100-kb pieces. The name sorts before Pa_A, so sources must be classified first.
"""
import pandas as pd
import pysam
import pytest

from cenrep.params import Params
from conftest import TEST_PARAMS, have
from gffcheck import check_gff
from synth import L, build_complex

TOOLS = ("FasTAN", "FAtoGDB", "ONEview", "makeblastdb", "blastn", "mafft", "cd-hit-est", "minimap2", "cmsearch")
pytestmark = pytest.mark.skipif(not have(*TOOLS), reason="external tools missing")


@pytest.fixture(scope="module")
def borrowed_run(tmp_path_factory):
    from cenrep.pipeline import Pipeline
    root = tmp_path_factory.mktemp("e2e_borrowed")
    cdir = root / "Synth"
    truth = build_complex(cdir)
    pieces = {}
    with pysam.FastaFile(str(cdir / "Pa_A.fa")) as fa, open(cdir / "Ctg_A.fa", "w") as out:
        for ci in range(1, 5):
            src = f"Pa_A{ci:02d}"
            cuts = [0, 150_000, L] if ci < 4 else [0, 100_000, 200_000, L]
            for k, (s, e) in enumerate(zip(cuts, cuts[1:])):
                name = f"Ctg_A{ci:02d}{'abc'[k]}"
                out.write(f">{name}\n{fa.fetch(src, s, e)}\n")
                pieces[name] = (src, s, e)
    pysam.faidx(str(cdir / "Ctg_A.fa"))
    p = Params(**dict(TEST_PARAMS, classify_from="Ctg_A=Pa_A", evalue_ref_len=L))
    Pipeline(cdir, 4, root / "cache", p).run()
    return cdir, truth, pieces, p


def _tsv(path, **kw):
    return pd.read_csv(path, sep="\t", **kw)


def _summary(cdir):
    return _tsv(cdir / "centromere" / "Ctg_A_centromere_summary.tsv", dtype=str).set_index("chrom")


def test_borrowed_family_table(borrowed_run):
    cdir, _, _, _ = borrowed_run
    work = cdir / "centromere" / "work"
    fam = _tsv(work / "families" / "families.tsv").set_index("fam_id")
    src = _tsv(work / "classify" / "Pa_A.families.tsv").set_index("family")
    tgt = _tsv(work / "classify" / "Ctg_A.families.tsv").set_index("family")
    assert tgt.index.is_unique
    assert set(tgt.index[tgt.cls == "centromeric"]) == set(src.index[src.cls == "centromeric"]) & set(tgt.index)
    main = tgt[(fam.loc[tgt.index, "consensus_len"] - 150).abs().to_numpy() <= 3].bp.idxmax()
    assert tgt.loc[main, "cls"] == "centromeric"                                    # borrowed, whatever its own metrics
    assert tgt.loc[main, "reason"].startswith("borrowed from Pa_A") and tgt.loc[main, "strength"] == src.loc[main, "strength"]
    decoy = [f for f in tgt.index if abs(fam.loc[f, "consensus_len"] - 90) <= 3]
    assert decoy and set(tgt.loc[decoy, "cls"]) == {"other"}


def test_borrowed_calls_match_truth(borrowed_run):
    cdir, truth, pieces, _ = borrowed_run
    s = _summary(cdir)
    for name, (src, a, b) in pieces.items():
        cen = truth["genomes"]["Pa_A"][src]["cen"]
        if cen is None:
            assert s.loc[name, "evidence"] == "none", name
            continue
        ts, te = max(cen[0], a) - a, min(cen[1], b) - a
        assert s.loc[name, "evidence"] == "satellite", name
        assert abs(int(s.loc[name, "cen_start"]) - 1 - ts) <= 150, name
        assert abs(int(s.loc[name, "cen_end"]) - te) <= 150, name
        assert "contig_end" in s.loc[name, "flags"] and "near_end" not in s.loc[name, "flags"], name
    assert "long_unit_tandem" not in set(s.evidence) and "weak_satellite" not in set(s.evidence)


def test_borrowed_long_unit_scan_stays_near_calls(borrowed_run):
    cdir, _, _, p = borrowed_run
    s = _summary(cdir)
    lu = _tsv(cdir / "centromere" / "work" / "centromeres" / "Ctg_A.longunit.tsv", dtype={"chrom": str})
    for r in lu.itertuples(index=False):
        c = s.loc[r.chrom]
        assert c.evidence != "none"
        assert r.start >= int(c.cen_start) - 1 - p.lu_flank and r.end <= int(c.cen_end) + p.lu_flank
    assert ((lu.chrom == "Ctg_A03b") & ((lu.unit - 3000).abs() <= 150)).any()


def test_borrowed_hor_gap_and_gff(borrowed_run):
    from cenrep.genomes import load_fai
    cdir, _, _, _ = borrowed_run
    out = cdir / "centromere"
    u = _tsv(out / "Ctg_A_hor_units.tsv")
    c1 = u[u.chrom == "Ctg_A01a"]
    assert ((c1.n_monomers == 4) & (c1.max_tandem_copies >= 14)).any()
    assert "gaps:1" in _summary(cdir).loc["Ctg_A03b", "flags"]
    assert check_gff(out / "Ctg_A_centromeric_repeat.gff", load_fai(cdir / "Ctg_A.fa"))


def test_family_table_and_figure_cover_all_genomes(borrowed_run):
    cdir, _, _, _ = borrowed_run
    f = _tsv(cdir / "centromere" / "satellite_families.tsv")
    assert set(f.genome) == {"Ctg_A", "Pa_A", "Pb_B", "Tx_A", "Tx_B"} and len(f) == 5 * f.fam_id.nunique()
    assert (cdir / "centromere" / "centromere_overview.pdf").stat().st_size > 5_000
