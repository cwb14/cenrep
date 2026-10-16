import numpy as np
import pandas as pd
import pysam
import pytest

from cenrep.genomes import discover_genomes
from cenrep.params import Params
from cenrep.seeds import COLS, cache_path, chunks, dedupe, seed_genomes
from cenrep.util import union_bp
from conftest import have


def test_chunks_cover_with_overlap():
    assert list(chunks({"a": 25, "b": 7}, 10, 3)) == [
        ("a", 0, 10), ("a", 7, 17), ("a", 14, 24), ("a", 21, 25), ("b", 0, 7)]


def test_dedupe_keeps_best_of_overlapping_duplicates():
    r1 = ("c", 100, 1000, 150, 6.0, 150, 90, 1, 500, 1.9, "A")
    r2 = ("c", 120, 1000, 151, 5.8, 151, 90, 1, 600, 1.9, "B")   # same array seen from the next chunk
    r3 = ("c", 100, 1000, 300, 3.0, 300, 90, 1, 400, 1.9, "C")   # different period: kept
    assert sorted(r[-1] for r in dedupe([r1, r2, r3])) == ["B", "C"]


def test_dedupe_period_tolerance_scales_with_the_unit():
    """A FasTAN unit is a mean over the copies aligned, so the same 1-kb array seen from two chunks can differ by a
    few bp; at tol 0 (TRF) only 1 bp is allowed."""
    r1 = ("c", 0, 10_000, 1001, 10.0, 1001, 90, 0, 9000, 0, "A")
    r2 = ("c", 0, 9_000, 1004, 9.0, 1004, 90, 0, 8100, 0, "B")
    assert [r[-1] for r in dedupe([r1, r2], 0.05)] == ["A"]
    assert sorted(r[-1] for r in dedupe([r1, r2])) == ["A", "B"]


def test_trf_cache_name_is_unchanged(tmp_path):
    """TRF caches written by earlier cenrep versions (hours per large genome) must still be found."""
    assert cache_path(tmp_path, "abc", Params(seed_engine="trf")).name == "abc.25429390.trf.tsv.gz"


def test_cache_key_depends_on_engine_and_chunking(tmp_path):
    a = cache_path(tmp_path, "abc", Params())
    assert a == cache_path(tmp_path, "abc", Params()) and a.name.endswith(".fastan.tsv.gz")
    assert a != cache_path(tmp_path, "abc", Params(seed_chunk=1_000_000))
    assert a != cache_path(tmp_path, "abc", Params(seed_engine="trf"))


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview"), reason="FasTAN tools missing")
def test_fastan_finds_planted_satellite_across_chunk_boundary(synth_seeds, synth_complex):
    cdir, truth = synth_complex
    df = synth_seeds["Pa_A"]
    s, e = truth["genomes"]["Pa_A"]["Pa_A01"]["cen"]      # 130,000-169,000 crosses the 145-150 kb chunk overlap
    k = (df.period / 150).round()
    mult = (k >= 1) & ((df.period / 150 - k).abs() <= 0.03 * k)
    hit = df[(df.chrom == "Pa_A01") & mult & (df.start < e) & (df.end > s)]
    assert union_bp(zip(hit.start, hit.end)) > 0.9 * (e - s)
    assert (hit.consensus.str.len() - hit.period).abs().max() <= 3
    for p in hit.itertuples():                            # no chunk duplicates left
        dup = hit[(hit.index != p.Index) & ((hit.period - p.period).abs() <= 0.05 * p.period)
                  & (hit.start < p.end) & (hit.end > p.start)]
        for q in dup.itertuples():
            ov = min(p.end, q.end) - max(p.start, q.start)
            assert ov < 0.5 * min(p.end - p.start, q.end - q.start)


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview"), reason="FasTAN tools missing")
def test_seed_cache_is_reused(synth_complex, synth_seed_cache, test_params):
    cdir, _ = synth_complex
    genomes = discover_genomes(cdir)
    before = {p: p.stat().st_mtime for p in synth_seed_cache.glob("*.fastan.tsv.gz")}
    seed_genomes(genomes, synth_seed_cache, 2, test_params)
    after = {p: p.stat().st_mtime for p in synth_seed_cache.glob("*.fastan.tsv.gz")}
    assert before == after and len(after) == 4


def _one_genome(tmp_path):
    from cenrep.genomes import Genome
    fa = tmp_path / "g.fa"
    fa.write_text(">c\n" + "ACGTTGCA" * 500 + "\n")
    pysam.faidx(str(fa))
    return Genome("g", fa, {"c": 4000})


def test_seeds_wait_for_lock_held_by_live_process(tmp_path, test_params, monkeypatch):
    """A second process must not recompute (or corrupt) a genome another live process is running."""
    import gzip
    import os
    import threading
    import time
    import cenrep.seeds as seeds
    from cenrep.util import md5_file
    monkeypatch.setattr(seeds, "POLL", 0.1)
    g = _one_genome(tmp_path)
    cache = tmp_path / "cache"
    cache.mkdir()
    path = seeds.cache_path(cache, md5_file(g.fasta), test_params)
    lock = seeds.lock_path(path)
    lock.write_text(str(os.getpid()))                       # held by a live process (this one)
    res = {}
    t = threading.Thread(target=lambda: res.update(seeds.seed_genomes([g], cache, 1, test_params)))
    t.start()
    time.sleep(0.5)
    assert t.is_alive()                                      # still waiting on the lock
    with gzip.open(path, "wt") as f:                          # the other process finishes
        f.write("\t".join(seeds.COLS) + "\nFROM_OTHER_PROCESS\n")
    lock.unlink()
    t.join(5)
    assert not t.is_alive() and "FROM_OTHER_PROCESS" in gzip.open(res["g"], "rt").read()


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview"), reason="FasTAN tools missing")
def test_seeds_steal_stale_lock(tmp_path, test_params):
    import socket
    import cenrep.seeds as seeds
    from cenrep.util import md5_file
    g = _one_genome(tmp_path)
    cache = tmp_path / "cache"
    cache.mkdir()
    path = seeds.cache_path(cache, md5_file(g.fasta), test_params)
    seeds.lock_path(path).write_text(f"{socket.gethostname()} 999999999")   # dead pid on this host
    out = seeds.seed_genomes([g], cache, 1, test_params)
    assert out["g"].exists() and not seeds.lock_path(path).exists()


def test_seeds_never_steal_a_lock_held_on_another_host(tmp_path, test_params, monkeypatch):
    """On a shared filesystem (SLURM nodes) a PID from another host cannot be checked here: its lock is never stolen."""
    import threading
    import time
    import cenrep.seeds as seeds
    from cenrep.util import md5_file
    monkeypatch.setattr(seeds, "POLL", 0.1)
    g = _one_genome(tmp_path)
    cache = tmp_path / "cache"
    cache.mkdir()
    path = seeds.cache_path(cache, md5_file(g.fasta), test_params)
    lock = seeds.lock_path(path)
    lock.write_text("some-other-node 999999999")             # dead here, but on another host
    t = threading.Thread(target=lambda: seeds.seed_genomes([g], cache, 1, test_params), daemon=True)
    t.start()
    time.sleep(0.5)
    assert t.is_alive() and lock.read_text() == "some-other-node 999999999"


def test_seed_lock_records_host_and_pid(tmp_path):
    import os
    import socket
    import cenrep.seeds as seeds
    lock = tmp_path / "x.lock"
    assert seeds._acquire(lock) and lock.read_text() == f"{socket.gethostname()} {os.getpid()}"


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview"), reason="FasTAN tools missing")
def test_seeds_cache_each_genome_as_soon_as_it_is_done(tmp_path, test_params):
    """A failure in a later genome must not lose the finished genomes of the same run (hours on real data)."""
    from cenrep.genomes import Genome
    from cenrep.seeds import _compute
    rng = np.random.default_rng(3)
    fa = tmp_path / "a.fa"
    fa.write_text(">c1\n" + "".join(rng.choice(list("ACGT"), 30_000)) + "\n")
    pysam.faidx(str(fa))
    good = Genome("a", fa, {"c1": 30_000})
    bad = Genome("b", tmp_path / "missing.fa", {"c1": 30_000})
    pa, pb = tmp_path / "a.fastan.tsv.gz", tmp_path / "b.fastan.tsv.gz"
    with pytest.raises(Exception):
        _compute({pa: good, pb: bad}, tmp_path, 2, test_params)
    assert pa.exists() and not pb.exists()
    assert list(pd.read_csv(pa, sep="\t").columns) == COLS


def test_genome_without_chunks_is_cached(tmp_path, test_params):
    """Zero-length sequences give no chunks; the cache must still be written, or seed_genomes loops. (The trf
    engine needs no tool when there is nothing to scan.)"""
    from dataclasses import replace
    from cenrep.genomes import Genome
    from cenrep.seeds import _compute
    path = tmp_path / "z.trf.tsv.gz"
    _compute({path: Genome("z", tmp_path / "z.fa", {"c1": 0})}, tmp_path, 1, replace(test_params, seed_engine="trf"))
    assert path.exists() and pd.read_csv(path, sep="\t").empty
