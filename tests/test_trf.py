import numpy as np
import pandas as pd
import pysam
import pytest

from cenrep.genomes import discover_genomes
from cenrep.trf import COLS, chunks, dedupe, load_trf, parse_ngs, trf_genomes
from cenrep.util import union_bp
from conftest import have

NGS = ("@c1\n"
       "11 60 10 5.0 10 95 0 90 25 25 25 25 2.0 ACGTACGTAC ACGTACGTACACGT FLANKLEFT FLANKRIGHT\n"
       "garbage line\n")


def test_parse_ngs_offsets_and_skips_bad_lines():
    assert parse_ngs(NGS, "c1", 1000) == [("c1", 1010, 1060, 10, 5.0, 10, 95, 0, 90, 2.0, "ACGTACGTAC")]


def test_chunks_cover_with_overlap():
    assert list(chunks({"a": 25, "b": 7}, 10, 3)) == [
        ("a", 0, 10), ("a", 7, 17), ("a", 14, 24), ("a", 21, 25), ("b", 0, 7)]


def test_dedupe_keeps_best_of_overlapping_duplicates():
    r1 = ("c", 100, 1000, 150, 6.0, 150, 90, 1, 500, 1.9, "A")
    r2 = ("c", 120, 1000, 151, 5.8, 151, 90, 1, 600, 1.9, "B")   # same array seen from the next chunk
    r3 = ("c", 100, 1000, 300, 3.0, 300, 90, 1, 400, 1.9, "C")   # different period: kept
    assert sorted(r[-1] for r in dedupe([r1, r2, r3])) == ["B", "C"]


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


@pytest.mark.skipif(not have("trf"), reason="trf missing")
def test_trf_cache_is_reused(synth_complex, synth_trf_cache, test_params):
    cdir, _ = synth_complex
    genomes = discover_genomes(cdir)
    before = {p: p.stat().st_mtime for p in synth_trf_cache.glob("*.trf.tsv.gz")}
    trf_genomes(genomes, synth_trf_cache, 2, test_params)
    after = {p: p.stat().st_mtime for p in synth_trf_cache.glob("*.trf.tsv.gz")}
    assert before == after and len(after) == 4


def _one_genome(tmp_path):
    import pysam
    from cenrep.genomes import Genome
    fa = tmp_path / "g.fa"
    fa.write_text(">c\n" + "ACGTTGCA" * 500 + "\n")
    pysam.faidx(str(fa))
    return Genome("g", fa, {"c": 4000})


def test_trf_waits_for_lock_held_by_live_process(tmp_path, test_params, monkeypatch):
    """A second process must not recompute (or corrupt) a genome another live process is running."""
    import gzip
    import os
    import threading
    import time
    import cenrep.trf as trf
    from cenrep.util import md5_file
    monkeypatch.setattr(trf, "POLL", 0.1)
    g = _one_genome(tmp_path)
    cache = tmp_path / "cache"
    cache.mkdir()
    path = trf.cache_path(cache, md5_file(g.fasta), test_params)
    lock = trf.lock_path(path)
    lock.write_text(str(os.getpid()))                       # held by a live process (this one)
    res = {}
    t = threading.Thread(target=lambda: res.update(trf.trf_genomes([g], cache, 1, test_params)))
    t.start()
    time.sleep(0.5)
    assert t.is_alive()                                      # still waiting on the lock
    with gzip.open(path, "wt") as f:                          # the other process finishes
        f.write("\t".join(trf.COLS) + "\nFROM_OTHER_PROCESS\n")
    lock.unlink()
    t.join(5)
    assert not t.is_alive() and "FROM_OTHER_PROCESS" in gzip.open(res["g"], "rt").read()


@pytest.mark.skipif(not have("trf"), reason="trf missing")
def test_trf_steals_stale_lock(tmp_path, test_params):
    import cenrep.trf as trf
    from cenrep.util import md5_file
    g = _one_genome(tmp_path)
    cache = tmp_path / "cache"
    cache.mkdir()
    path = trf.cache_path(cache, md5_file(g.fasta), test_params)
    import socket
    trf.lock_path(path).write_text(f"{socket.gethostname()} 999999999")   # dead pid on this host
    out = trf.trf_genomes([g], cache, 1, test_params)
    assert out["g"].exists() and not trf.lock_path(path).exists()


def test_trf_never_steals_a_lock_held_on_another_host(tmp_path, test_params, monkeypatch):
    """On a shared filesystem (SLURM nodes) a PID from another host cannot be checked here: its lock is never stolen."""
    import threading
    import time
    import cenrep.trf as trf
    from cenrep.util import md5_file
    monkeypatch.setattr(trf, "POLL", 0.1)
    g = _one_genome(tmp_path)
    cache = tmp_path / "cache"
    cache.mkdir()
    path = trf.cache_path(cache, md5_file(g.fasta), test_params)
    lock = trf.lock_path(path)
    lock.write_text("some-other-node 999999999")             # dead here, but on another host
    t = threading.Thread(target=lambda: trf.trf_genomes([g], cache, 1, test_params), daemon=True)
    t.start()
    time.sleep(0.5)
    assert t.is_alive() and lock.read_text() == "some-other-node 999999999"


def test_trf_lock_records_host_and_pid(tmp_path):
    import os
    import socket
    import cenrep.trf as trf
    lock = tmp_path / "x.lock"
    assert trf._acquire(lock) and lock.read_text() == f"{socket.gethostname()} {os.getpid()}"


@pytest.mark.skipif(not have("trf"), reason="trf missing")
def test_trf_caches_each_genome_as_soon_as_it_is_done(tmp_path, test_params):
    """A failure in a later genome must not lose the finished genomes of the same TRF run (hours on real data)."""
    from cenrep.genomes import Genome
    from cenrep.trf import _compute
    rng = np.random.default_rng(3)
    fa = tmp_path / "a.fa"
    fa.write_text(">c1\n" + "".join(rng.choice(list("ACGT"), 30_000)) + "\n")
    pysam.faidx(str(fa))
    good = Genome("a", fa, {"c1": 30_000})
    bad = Genome("b", tmp_path / "missing.fa", {"c1": 30_000})
    pa, pb = tmp_path / "a.trf.tsv.gz", tmp_path / "b.trf.tsv.gz"
    with pytest.raises(Exception):
        _compute({pa: good, pb: bad}, tmp_path, 2, test_params)
    assert pa.exists() and not pb.exists()
    assert list(pd.read_csv(pa, sep="\t").columns) == COLS


def test_trf_genome_without_chunks_is_cached(tmp_path, test_params):
    """Review M4: zero-length sequences give no chunks; the cache must still be written, or trf_genomes loops."""
    from cenrep.genomes import Genome
    from cenrep.trf import _compute
    path = tmp_path / "z.trf.tsv.gz"
    _compute({path: Genome("z", tmp_path / "z.fa", {"c1": 0})}, tmp_path, 1, test_params)
    assert path.exists() and pd.read_csv(path, sep="\t").empty
