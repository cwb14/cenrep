import shutil

import pytest

from cenrep.params import Params
from synth import build_complex

TEST_PARAMS = dict(
    trf_chunk=150_000, trf_overlap=5_000, min_fam_span=5_000, hmm_sample=100, hmm_sample_bp=100_000,
    chrom_min_bp=5_000, cluster_gap=50_000, subtel_min_bp=20_000, strong_cluster_bp=20_000,
    low_sat_bp=10_000, fallback_min_bp=2_000, lu_flank=20_000, lu_window=100_000, lu_window_overlap=20_000,
)


def have(*tools) -> bool:
    return all(shutil.which(t) for t in tools)


@pytest.fixture(scope="session")
def test_params():
    return Params(**TEST_PARAMS)


@pytest.fixture(scope="session")
def synth_complex(tmp_path_factory):
    cdir = tmp_path_factory.mktemp("synth") / "Synth"
    truth = build_complex(cdir)
    return cdir, truth


@pytest.fixture(scope="session")
def synth_trf_cache(synth_complex, test_params, tmp_path_factory):
    from cenrep.genomes import discover_genomes
    from cenrep.trf import trf_genomes
    cache = tmp_path_factory.mktemp("trfcache")
    trf_genomes(discover_genomes(synth_complex[0]), cache, 4, test_params)
    return cache


@pytest.fixture(scope="session")
def synth_trf(synth_complex, synth_trf_cache, test_params):
    from cenrep.genomes import discover_genomes
    from cenrep.trf import load_trf, trf_genomes
    paths = trf_genomes(discover_genomes(synth_complex[0]), synth_trf_cache, 4, test_params)
    return {g: load_trf(p) for g, p in paths.items()}


@pytest.fixture(scope="session")
def synth_families(synth_complex, synth_trf, test_params, tmp_path_factory):
    from cenrep.families import build_families
    from cenrep.genomes import discover_genomes
    wd = tmp_path_factory.mktemp("families")
    fam = build_families(synth_trf, discover_genomes(synth_complex[0]), wd, 4, test_params)
    return fam, wd


@pytest.fixture(scope="session")
def synth_annot(synth_complex, synth_families, test_params, tmp_path_factory):
    import pyhmmer
    from cenrep.annotate import annotate_genome
    from cenrep.genomes import discover_genomes
    fam, wd = synth_families
    with pyhmmer.plan7.HMMFile(wd / "families.hmm") as f:
        hmms = list(f)
    g = [x for x in discover_genomes(synth_complex[0]) if x.name == "Pa_A"][0]
    out = tmp_path_factory.mktemp("annot")
    return annotate_genome(g, hmms, dict(zip(fam.fam_id, fam.consensus_len)), 4, test_params,
                           out / "Pa_A.monomers.tsv.gz", out / "Pa_A.arrays.tsv")
