"""S1 engine: TRF (Benson 1999) tandem repeats with their consensus, one chunk at a time."""
import os
import tempfile

from .params import TRF_ARGS
from .util import log, run

TOOLS = ["trf"]
KEY = [TRF_ARGS]    # cache key (with the chunking) as in cenrep < 0.9, so TRF caches computed earlier stay valid
DUP_TOL = 0         # the same array seen from two chunks keeps its TRF period (within 1 bp)


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


def setup(tmpdir):
    return None


def chunk_rows(chrom: str, seq: str, offset: int, tmpdir, opt=None) -> list:
    fd, path = tempfile.mkstemp(suffix=".fa", dir=tmpdir)
    with os.fdopen(fd, "w") as f:
        f.write(f">{chrom}\n{seq}\n")
    try:
        p = run(["trf", path] + TRF_ARGS, cwd=tmpdir)
    finally:
        os.unlink(path)
    return parse_ngs(p.stdout, chrom, offset)
