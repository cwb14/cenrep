"""S1: TRF tandem-repeat discovery on chunked chromosomes, cached by genome md5."""
import gzip
import hashlib
import json
import os
import shutil
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import pysam

from .util import log, md5_file, run

# Benson's recommended weights (match 2, mismatch 7, indel 7, PM 80, PI 10), minscore 50,
# max period 2000 (TRF maximum), -l 6: chunks are <= 5 Mb so no array is longer.
TRF_ARGS = ["2", "7", "7", "80", "10", "50", "2000", "-l", "6", "-h", "-ngs"]
POLL = 15   # seconds between checks while another process holds a genome's TRF lock
COLS = ["chrom", "start", "end", "period", "copies", "consensus_size",
        "pct_match", "pct_indel", "score", "entropy", "consensus"]   # start/end 0-based half-open


def chunks(chroms: dict, size: int, overlap: int):
    for c, n in chroms.items():
        s = 0
        while s < n:
            e = min(n, s + size)
            yield c, s, e
            if e == n:
                break
            s = e - overlap


def parse_ngs(text: str, chrom: str, offset: int) -> list:
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("@"):
            continue
        p = line.split()
        try:
            rows.append((chrom, int(p[0]) - 1 + offset, int(p[1]) + offset, int(p[2]), float(p[3]),
                         int(p[4]), int(p[5]), int(p[6]), int(p[7]), float(p[12]), p[13].upper()))
        except (ValueError, IndexError):
            log.warning(f"TRF: skipping unparsable line on {chrom}: {line[:80]!r}")
    return rows


def _run_chunk(job):
    fasta, chrom, s, e, tmpdir = job
    with pysam.FastaFile(str(fasta)) as fa:
        seq = fa.fetch(chrom, s, e)
    fd, path = tempfile.mkstemp(suffix=".fa", dir=tmpdir)
    with os.fdopen(fd, "w") as f:
        f.write(f">{chrom}\n{seq}\n")
    try:
        p = run(["trf", path] + TRF_ARGS, cwd=tmpdir)
    finally:
        os.unlink(path)
    return str(fasta), parse_ngs(p.stdout, chrom, s)


def dedupe(rows: list) -> list:
    """Drop chunk-overlap duplicates: same chrom, period within 1 bp, overlap >= 50% of the shorter
    -> keep the higher TRF score."""
    rows = sorted(rows, key=lambda r: (r[0], r[1], r[2], r[3]))
    kept, active, cur = [], [], None
    for r in rows:
        if r[0] != cur:
            cur, active = r[0], []
        active = [i for i in active if kept[i][2] > r[1]]
        dup = None
        for i in active:
            k = kept[i]
            if abs(k[3] - r[3]) <= 1 and \
                    min(k[2], r[2]) - max(k[1], r[1]) >= 0.5 * min(k[2] - k[1], r[2] - r[1]):
                dup = i
                break
        if dup is None:
            kept.append(r)
            active.append(len(kept) - 1)
        elif r[8] > kept[dup][8]:
            kept[dup] = r
    return kept


def cache_path(cache_dir, md5: str, p) -> Path:
    """Cache file for one genome: keyed by the FASTA md5 and by the TRF settings that shape the result."""
    key = hashlib.md5(json.dumps([TRF_ARGS, p.trf_chunk, p.trf_overlap]).encode()).hexdigest()[:8]
    return Path(cache_dir) / f"{md5}.{key}.trf.tsv.gz"


def lock_path(path: Path) -> Path:
    return Path(str(path) + ".lock")


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _acquire(lock: Path) -> bool:
    """Create the lock file atomically; steal it if its owner process is dead."""
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        try:
            pid = int(lock.read_text().strip() or 0)
        except (OSError, ValueError):
            return False                          # being written right now; retry next poll
        if pid and not _alive(pid):
            log.warning(f"trf: removing stale lock {lock.name} (pid {pid} is gone)")
            lock.unlink(missing_ok=True)
            return _acquire(lock)
        return False
    os.write(fd, str(os.getpid()).encode())
    os.close(fd)
    return True


def _write_cache(path: Path, g, rows: list) -> None:
    kept = dedupe(rows)
    part = Path(f"{path}.{os.getpid()}.part")
    with gzip.open(part, "wt") as f:
        f.write("\t".join(COLS) + "\n")
        for r in kept:
            f.write("\t".join(map(str, r)) + "\n")
    part.rename(path)
    log.info(f"trf: {g.name} {len(kept)} TRF rows")


def _compute(todo: dict, cache_dir: Path, threads: int, p) -> None:
    """Jobs run genome by genome and results arrive in job order, so each genome is cached as soon as its
    last chunk is back: a later failure loses only the unfinished genomes."""
    tmp = tempfile.mkdtemp(prefix="trf_tmp_", dir=cache_dir)
    jobs, left, target = [], {}, {}
    for path, g in todo.items():
        mine = [(g.fasta, c, s, e, tmp) for c, s, e in chunks(g.chroms, p.trf_chunk, p.trf_overlap)]
        jobs += mine
        left[str(g.fasta)], target[str(g.fasta)] = len(mine), (path, g)
    log.info(f"trf: {len(todo)} genome(s), {len(jobs)} chunks, {threads} workers")
    rows = {fa: [] for fa in left}
    for fa in [fa for fa, n in left.items() if not n]:          # only zero-length sequences: nothing to run
        _write_cache(*target[fa], rows.pop(fa))
    try:
        with ProcessPoolExecutor(threads) as ex:
            for i, (fa, r) in enumerate(ex.map(_run_chunk, jobs, chunksize=1), 1):
                rows[fa].extend(r)
                left[fa] -= 1
                log.debug(f"trf: chunk {i}/{len(jobs)}")
                if not left[fa]:
                    _write_cache(*target[fa], rows.pop(fa))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def trf_genomes(genomes, cache_dir, threads: int, p) -> dict:
    """Run TRF on every genome not cached yet (md5-keyed). Returns {genome name: cached .tsv.gz path}.
    Safe for concurrent runs sharing the cache: a per-genome lock file makes one process compute a genome
    while others wait for its result (stale locks of dead processes are stolen)."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    out, todo = {}, {}
    for g in genomes:
        path = cache_path(cache_dir, md5_file(g.fasta), p)
        out[g.name] = path
        if path.exists():
            log.info(f"trf: {g.name} cached ({path.name})")
        else:
            todo.setdefault(path, g)          # md5-identical genomes share one run
    waiting_logged = set()
    while True:
        todo = {path: g for path, g in todo.items() if not path.exists()}
        if not todo:
            return out
        mine = {path: g for path, g in todo.items() if _acquire(lock_path(path))}
        if mine:
            try:
                _compute(mine, cache_dir, threads, p)
            finally:
                for path in mine:
                    lock_path(path).unlink(missing_ok=True)
            continue
        for path, g in todo.items():
            if path not in waiting_logged:
                log.info(f"trf: {g.name} is being computed by another run; waiting for {path.name}")
                waiting_logged.add(path)
        time.sleep(POLL)


def load_trf(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"chrom": str, "consensus": str})
