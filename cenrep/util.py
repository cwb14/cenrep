"""Shared helpers: logging, subprocesses, sequences, FASTA text, intervals."""
import hashlib
import logging
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger("cenrep")


class CenrepError(RuntimeError):
    """Fatal, user-facing error (bad input, missing tool, failed command)."""


def setup_logging(verbose: bool) -> None:
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S"))
    log.handlers[:] = [h]
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    log.propagate = False


def require_tools(names) -> None:
    missing = [n for n in names if shutil.which(n) is None]
    if missing:
        raise CenrepError("required tools not on PATH: " + ", ".join(missing))


def run(cmd, **kw) -> subprocess.CompletedProcess:
    cmd = [str(c) for c in cmd]
    p = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if p.returncode != 0:
        raise CenrepError(f"command failed (exit {p.returncode}): {' '.join(cmd)}\n{p.stderr[-2000:]}")
    return p


def md5_file(path, bufsize: int = 1 << 24) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(bufsize), b""):
            h.update(b)
    return h.hexdigest()


_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]


def canonical_rotation(seq: str) -> str:
    """Lexicographically smallest rotation of a circular monomer (upper-cased).
    It starts with the smallest k-mer, so it is the spec's 'smallest 16-mer' rule with ties resolved."""
    s = seq.upper()
    return min((s[i:] + s[:i] for i in range(len(s))), default=s)


def write_fasta(path, records) -> None:
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")


def parse_fasta_text(text: str) -> list:
    out, name, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def merge_intervals(iv, gap: int = 0) -> list:
    """Merge 0-based half-open intervals separated by <= gap bp. Returns sorted [start, end] lists."""
    out = []
    for s, e in sorted(iv):
        if out and s - out[-1][1] <= gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def union_bp(iv) -> int:
    return sum(e - s for s, e in merge_intervals(iv))


def clipped_bp(iv, lo: int, hi: int) -> int:
    """Union bp of intervals clipped to [lo, hi)."""
    return union_bp([(max(s, lo), min(e, hi)) for s, e in iv if min(e, hi) > max(s, lo)])


def overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def touch(path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("")
