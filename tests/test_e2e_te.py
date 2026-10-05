"""Retrotransposon evidence end to end: GENOME_DIR/tes/<genome>.gff3 with a centromere-concentrated clade (CRM) and a
uniformly spread one (Ale)."""
import numpy as np
import pandas as pd
import pytest

from cenrep.params import Params
from conftest import TEST_PARAMS, have
from synth import L, build_complex

TOOLS = ("trf", "makeblastdb", "blastn", "mafft", "cd-hit-est", "minimap2", "cmsearch")
pytestmark = pytest.mark.skipif(not have(*TOOLS), reason="external tools missing")


def write_tes(cdir, truth, crm_at):
    """Per chromosome: 40 CRM elements inside crm_at (lo, hi), 15 CRM and 60 Ale elements anywhere."""
    (cdir / "tes").mkdir()
    rng = np.random.default_rng(7)
    for g, chroms in truth["genomes"].items():
        lines = ["##gff-version 3"]
        for c in chroms:
            mids = [("CRM", x) for x in rng.integers(crm_at[0], crm_at[1], 40)] + \
                   [("CRM", x) for x in rng.integers(0, L, 15)] + [("Ale", x) for x in rng.integers(0, L, 60)]
            for i, (clade, m) in enumerate(mids):
                s, e = max(1, int(m) - 2000), min(L, int(m) + 2000)
                lines.append(f"{c}\tLTRquest\tLTR_retrotransposon\t{s}\t{e}\t.\t+\t.\tID={c}_{i};clade={clade}")
        (cdir / "tes" / f"{g}.gff3").write_text("\n".join(lines) + "\n")


def _run(root, crm_at):
    from cenrep.pipeline import Pipeline
    cdir = root / "Synth"
    truth = build_complex(cdir)
    write_tes(cdir, truth, crm_at)
    Pipeline(cdir, 4, root / "cache", Params(**TEST_PARAMS)).run()
    summ = {g: pd.read_csv(cdir / "centromere" / f"{g}_centromere_summary.tsv", sep="\t", dtype=str).set_index("chrom")
            for g in truth["genomes"]}
    fams = pd.read_csv(cdir / "centromere" / "satellite_families.tsv", sep="\t", dtype=str)
    return truth, summ, fams


def test_retrotransposon_core_at_the_satellite_centromere_confirms_it(tmp_path):
    truth, summ, fams = _run(tmp_path, (135_000, 165_000))
    for g, s in summ.items():
        for c in (f"{g}01", f"{g}02", f"{g}03"):
            assert s.loc[c, "evidence"] == "satellite" and "te_concordant" in s.loc[c, "flags"], (g, c)
            assert s.loc[c, "confidence"] == "high" and s.loc[c, "te_clade"] == "CRM"
        c4 = s.loc[f"{g}04"]                                   # no satellite: the retrotransposon core places it
        assert c4.evidence == "retrotransposon" and c4.confidence == "medium"
        assert 130_000 <= int(c4.cen_start) and int(c4.cen_end) <= 175_000
    cen = fams[(fams.consensus_len.astype(int) - 150).abs() <= 3]
    assert set(cen["class"]) == {"centromeric"} and set(cen["dominant"]) == {"True"}


def test_satellite_away_from_the_retrotransposon_centromere_is_rejected(tmp_path):
    truth, summ, fams = _run(tmp_path, (30_000, 50_000))
    for g, s in summ.items():
        for c in truth["genomes"][g]:
            r = s.loc[c]
            assert r.evidence == "retrotransposon" and r.te_clade == "CRM", (g, c)
            assert 25_000 <= int(r.cen_start) and int(r.cen_end) <= 55_000, (g, c)
    cen = fams[(fams.consensus_len.astype(int) - 150).abs() <= 3]
    assert set(cen["class"]) == {"other"} and cen.reason.str.contains("retrotransposon").all()
