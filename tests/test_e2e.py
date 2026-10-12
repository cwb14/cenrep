import shutil
import subprocess
import sys

import pandas as pd
import pytest

from cenrep.params import Params
from cenrep.util import CenrepError
from conftest import TEST_PARAMS, have
from gffcheck import check_gff
from synth import build_complex

TOOLS = ("trf", "makeblastdb", "blastn", "mafft", "cd-hit-est", "minimap2")
pytestmark = pytest.mark.skipif(not have(*TOOLS), reason="external tools missing")


@pytest.fixture(scope="module")
def synth_run(tmp_path_factory):
    from cenrep.pipeline import Pipeline
    root = tmp_path_factory.mktemp("e2e")
    cdir = root / "Synth"
    truth = build_complex(cdir)
    Pipeline(cdir, 4, root / "cache", Params(**TEST_PARAMS)).run()
    return cdir, truth, root


def _summary(cdir, g):
    return pd.read_csv(cdir / "centromere" / f"{g}_centromere_summary.tsv", sep="\t", dtype=str).set_index("chrom")


def _gff(cdir, g):
    from cenrep.genomes import load_fai
    return check_gff(cdir / "centromere" / f"{g}_centromeric_repeat.gff", load_fai(cdir / f"{g}.fa"))


def test_outputs_exist_and_gffs_valid(synth_run):
    cdir, truth, _ = synth_run
    out = cdir / "centromere"
    for f in ("satellite_families.tsv", "satellite_families.fa", "satellite_families.hmm", "centromere_overview.pdf"):
        assert (out / f).stat().st_size > 0
    for g in truth["genomes"]:
        assert _gff(cdir, g)


def test_centromere_boundaries_and_none(synth_run):
    cdir, truth, _ = synth_run
    for g, chroms in truth["genomes"].items():
        s = _summary(cdir, g)
        for c, t in chroms.items():
            if t["cen"] is None:
                assert s.loc[c, "evidence"] == "none"
                continue
            assert abs(int(s.loc[c, "cen_start"]) - 1 - t["cen"][0]) <= 150, (g, c)
            assert abs(int(s.loc[c, "cen_end"]) - t["cen"][1]) <= 150, (g, c)


def test_decoy_family_is_other(synth_run):
    cdir, _, _ = synth_run
    f = pd.read_csv(cdir / "centromere" / "satellite_families.tsv", sep="\t")
    dec = f[(f.consensus_len - 90).abs() <= 3]
    assert len(dec) and set(dec["class"]) == {"other"}
    cen = f[(f.consensus_len - 150).abs() <= 3]
    assert set(cen["class"]) == {"centromeric"}


def test_tandem_hor_chr1_and_duplication_chr2(synth_run):
    cdir, _, _ = synth_run
    u = pd.read_csv(cdir / "centromere" / "Pa_A_hor_units.tsv", sep="\t")
    c1 = u[u.chrom == "Pa_A01"]
    assert ((c1.n_monomers == 4) & (c1.max_tandem_copies >= 14)).any()
    c2 = u[u.chrom == "Pa_A02"]
    assert ((c2.n_monomers >= 2) & (c2.occurrences >= 2) & (c2.max_tandem_copies == 1)).any()


def test_chr3_gap_microsat_and_long_unit(synth_run):
    cdir, _, _ = synth_run
    assert "gaps:1" in _summary(cdir, "Pa_A").loc["Pa_A03", "flags"]
    tr = [a for c, s, e, t, a in _gff(cdir, "Pa_A") if t == "tandem_repeat" and c == "Pa_A03"]
    assert any(a["period"] == "2" and a["unit_source"] == "trf" for a in tr)
    assert any(a["unit_source"] == "selfaln" and abs(int(a["period"]) - 3000) <= 150 for a in tr)


def test_resume_and_param_change(synth_run, monkeypatch):
    """Review I1: a checkpoint is reused only if the code version and the parameters of that stage and all
    earlier stages are unchanged; --force must name the earliest affected stage (or an earlier one)."""
    import cenrep.pipeline as pl
    from cenrep.pipeline import Pipeline
    cdir, _, root = synth_run
    Pipeline(cdir, 2, root / "cache", Params(**TEST_PARAMS)).run()          # all checkpoints: no-op
    changed = Params(**dict(TEST_PARAMS, cluster_gap=60_000))              # a classify parameter
    with pytest.raises(CenrepError, match="--force classify"):
        Pipeline(cdir, 2, root / "cache", changed)
    with pytest.raises(CenrepError, match="--force classify"):
        Pipeline(cdir, 2, root / "cache", changed, force="plot")          # too late a stage
    same = Params(**TEST_PARAMS)
    monkeypatch.setitem(pl.STAGE_VERSION, "classify", "9.9.9")             # classify code changed
    with pytest.raises(CenrepError, match="--force classify"):
        Pipeline(cdir, 2, root / "cache", same)
    monkeypatch.setitem(pl.STAGE_VERSION, "trf", "9.9.9")                  # TRF code changed
    with pytest.raises(CenrepError, match="--force trf"):
        Pipeline(cdir, 2, root / "cache", same)
    monkeypatch.undo()
    Pipeline(cdir, 2, root / "cache", changed, force="classify")          # accepted (deletes classify+ checkpoints)


def test_every_array_is_reported_with_family_and_role(synth_run):
    """<genome>_satellite_arrays.gff holds every array of every kept family, with the satellite_families.tsv name and
    the family's role in that genome."""
    from cenrep.classify import ROLES
    from cenrep.genomes import load_fai
    cdir, truth, _ = synth_run
    out = cdir / "centromere"
    fam = pd.read_csv(out / "satellite_families.tsv", sep="\t")
    assert "role" in fam.columns and set(fam.role) <= set(ROLES) | {"absent", "unassessed"}
    for g in truth["genomes"]:
        gff = out / f"{g}_satellite_arrays.gff"
        assert check_gff(gff, load_fai(cdir / f"{g}.fa"))
        lines = [l.split("\t") for l in gff.read_text().splitlines() if not l.startswith("#")]
        n = len(pd.read_csv(out / "work" / "annotate" / f"{g}.arrays.tsv", sep="\t"))
        assert len(lines) == n > 0 and all(l[2] == "satellite_DNA" for l in lines)
        att = [dict(kv.split("=", 1) for kv in l[8].split(";")) for l in lines]
        assert {a["family"] for a in att} <= set(fam.family)
        role = dict(zip(fam[fam.genome == g].family, fam[fam.genome == g].role))
        assert all(a["role"] == role[a["family"]] for a in att)


def test_rdna_column_in_family_table(synth_run):
    cdir, _, _ = synth_run
    f = pd.read_csv(cdir / "centromere" / "satellite_families.tsv", sep="\t", dtype=str)
    assert "rdna" in f.columns and set(f.rdna) == {"."}


def test_bgzipped_inputs_give_identical_results(synth_run, tmp_path):
    """The cleaned polyploid genomes are bgzip FASTA (.fa.gz + .fai + .gzi): every stage reads them through
    pysam, so results must equal the plain-FASTA run."""
    import pysam
    from cenrep.pipeline import Pipeline
    src, truth, _ = synth_run                     # inputs only: test_resume_and_param_change reruns its outputs
    cdir, gdir = tmp_path / "fa" / "Synth", tmp_path / "gz" / "Synth"     # same name: it prefixes family names
    cdir.mkdir(parents=True)
    gdir.mkdir(parents=True)
    for g in truth["genomes"]:
        for ext in (".fa", ".fa.fai"):
            shutil.copy(src / f"{g}{ext}", cdir / f"{g}{ext}")
        gz = gdir / f"{g}.fa.gz"
        pysam.tabix_compress(str(src / f"{g}.fa"), str(gz))
        pysam.faidx(str(gz))
    for d in (cdir, gdir):
        Pipeline(d, 4, tmp_path / "cache", Params(**TEST_PARAMS)).run()
    for g in truth["genomes"]:
        assert _summary(gdir, g).equals(_summary(cdir, g)), g
        assert (gdir / "centromere" / f"{g}_centromeric_repeat.gff").read_text() == \
            (cdir / "centromere" / f"{g}_centromeric_repeat.gff").read_text(), g
    fam = lambda d: pd.read_csv(d / "centromere" / "satellite_families.tsv", sep="\t")
    assert fam(gdir).equals(fam(cdir))


def test_cli_help_is_short():
    """Thresholds are internal constants (params.py), not options."""
    r = subprocess.run([sys.executable, "-m", "cenrep", "run", "--help"], capture_output=True, text=True)
    assert r.returncode == 0 and "GENOME_DIR" in r.stdout and "--outdir" in r.stdout
    assert "--cluster-gap" not in r.stdout and r.stdout.count("\n  -") <= 8


def test_cli_needs_a_genome_dir_and_fails_clearly(tmp_path):
    run = lambda *a: subprocess.run([sys.executable, "-m", "cenrep", "run", *a], capture_output=True, text=True)
    assert run().returncode == 2
    r = run(str(tmp_path / "empty"))
    assert r.returncode == 1 and "no *.fa genomes" in r.stderr
    r = run("--complex-dir", str(tmp_path / "empty"))                      # hidden alias of GENOME_DIR
    assert r.returncode == 1 and "no *.fa genomes" in r.stderr
