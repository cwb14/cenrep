# cenrep Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Annotate centromeric satellite arrays, monomers, HORs and other centromeric tandem arrays in all 40 genomes of the 10 allo10 complexes, writing `<Complex>/centromere/<Genome>_<L>_centromeric_repeat.gff` plus summaries.

**Architecture:** A Python package `cenrep` runs eight checkpointed stages per complex:

1. TRF discovery, chunked and cached by genome md5.
2. Greedy, rotation- and multimer-aware family clustering with BLAST, giving profile HMMs.
3. Genome-wide nhmmer monomer calls, chained into arrays.
4. Per-chromosome choice of the centromeric family.
5. Centromere regions, plus a minimap2 self-alignment scan for long-unit arrays.
6. HORs from cd-hit variant classes and maximal repeats found with a suffix array.
7. GFF3 and TSV output.
8. PDF figure.

A bash launcher runs the complexes one after another in a detached job.

**Tech Stack:** Python 3.10 with numpy, pandas, pysam, pyhmmer 0.12, matplotlib and pytest; TRF 4.10, BLAST+ 2.16, MAFFT 7.525, cd-hit-est 4.8.1 and minimap2. Everything comes from the existing mamba env `synLTR`.

**Spec:** `clean/cenrep/docs/2026-09-25-centromere-design.md`. The literature table is `clean/cenrep/docs/literature_cen_satellites.md`.

## Global Constraints

- **Env:** `ENV=/home/chris/bin/mambaforge/envs/synLTR`. Every command runs with `PATH=$ENV/bin:$PATH`.
- **Python:** 3.10.18, so write no 3.11+ syntax.
- **No installs** of any kind. If something is missing, stop and report exactly what is needed.
- **Threads:** at most 32 in total (`--threads 32`).
- **Long jobs:** launch with `setsid nohup ... > log 2>&1 < /dev/null &`, then verify PPID 1 with `ps -o pid,ppid,cmd`. Never use the Bash tool's `run_in_background`.
- **Memory:** never load a whole genome. The largest unit held in memory is one chromosome, fetched with pysam through the `.fai`.
- **Where writes go:**
  - `clean/<Complex>/centromere/**`
  - `clean/centromere_cache/**`
  - `clean/cenrep/**`
  - `clean/run_centromere.sh`
  - `clean/memo_centromere.md`
  - `.fai` files for `clean/Cbursapastoris/*.fa`, via `samtools faidx`
- **Output naming:** `<Complex>/centromere/<Genome>_<L>_centromeric_repeat.gff` (GFF3), `<G>_centromere_summary.tsv`, `<G>_hor_units.tsv`, `satellite_families.{tsv,fa,hmm}`, `centromere_overview.pdf`, `cenrep.{log,exit,time}`, `work/`.
- **Coordinates:** 0-based half-open internally; GFF is 1-based closed (`start+1`, `end`).
- **Logging:** the default is minimal (start, stage milestones, done, errors); `-v` adds per-step progress.
- **Errors:** fail fast with `CenrepError` for critical problems. Warn and skip for minor ones (one bad TRF line, a missing gene GFF).
- **Figures:** no titles; panel letters A–D in the upper-left; Okabe–Ito colours; vector PDF; 7 pt fonts.
- **Authorship:** never credit Claude, Anthropic or AI anywhere (code, pyproject, memo). pyproject has no `authors` field.
- **No git:** `clean/` is not a git repository, so this plan has no commit steps. Do not `git init` without asking the user.
- **Tests:** `cd /data2/chris/tetraploid/allo10/clean/cenrep && PATH=$ENV/bin:$PATH python -m pytest -q`

## Review Focus

These are the five input classes the spec implies but no task's main tests exercise. Each has a pinned test in the task that owns the code.

1. **Satellite-free genome or complex** (B genome, Coffea). Empty arrays, no primary family, every chromosome `none`. Every stage must still write valid, header-only outputs and never crash on empty DataFrames. Pinned: Task 5 `test_empty_arrays`, Task 9 `test_write_genome_without_regions`, Task 10 `test_plot_without_families`.
2. **Pathological class strings in HOR search** (Mb arrays of near-identical monomers). Must stay fast, with no O(n²) blow-up. Pinned: Task 8 `test_homogeneous_run_is_not_a_hor_and_is_fast`.
3. **GFF attribute escaping.** Flags and reasons contain `=`, `;` and `,`, which must not corrupt the attribute column. Pinned: Task 9 `test_attrs_escaping`.
4. **Changed parameters on a resumed run.** Checkpoints must not silently mix parameter sets. Pinned: Task 11 `test_resume_and_param_change`.
5. **Chromosome with a region but too few classed monomers for a HOR search** (n < 4). Must be skipped cleanly. Pinned: Task 8 `test_run_hor_tiny_region`.

---

### Task 1: Package scaffold, params, utilities, inputs, synthetic fixture

**Files:**
- Create: `clean/cenrep/pyproject.toml`, `clean/cenrep/environment.yml`
- Create: `clean/cenrep/cenrep/__init__.py`, `cenrep/__main__.py` (stub for now), `cenrep/params.py`, `cenrep/util.py`, `cenrep/genomes.py`
- Create: `clean/cenrep/tests/synth.py`, `tests/conftest.py`, `tests/test_util.py`, `tests/test_genomes.py`, `tests/test_synth.py`

**Interfaces:**
- Produces:
  - `cenrep.params.Params` (dataclass, fields below), `HELP: dict[str, str]`, `to_dict(p) -> dict`
  - `cenrep.util`:
    - `CenrepError`, `log`, `setup_logging(verbose)`
    - `require_tools(names)`, `run(cmd, **kw) -> CompletedProcess`, `md5_file(path) -> str`
    - `revcomp(s)`, `canonical_rotation(seq) -> str`
    - `write_fasta(path, records)`, `parse_fasta_text(text) -> list[(name, seq)]`
    - `merge_intervals(iv, gap=0) -> list[[s, e]]`, `union_bp(iv) -> int`, `clipped_bp(iv, lo, hi) -> int`, `overlap(a0, a1, b0, b1) -> int`
    - `touch(path)`
  - `cenrep.genomes`:
    - `Genome(name, fasta: Path, chroms: dict[str, int])`
    - `load_fai(fasta) -> dict`, `discover_genomes(complex_dir) -> list[Genome]`
    - `find_gene_gff(complex_dir, genome) -> Path|None`, `load_genes(complex_dir, genome) -> dict[str, np.ndarray]|None`, `count_in(mids, start, end) -> int`
  - `tests/synth.py`: `rseq(rng, n)`, `mutate(rng, s, rate)`, `rc(s)`, `build_complex(outdir, seed=11) -> truth dict`, constants `L=300_000`, `CEN_START=130_000`, `MON=150`
  - `tests/conftest.py`: fixtures `synth_complex` (session; returns `(cdir, truth)`) and `test_params` (session; `Params(**TEST_PARAMS)`), plus helper `have(*tools)`

- [ ] **Step 1: Create project metadata**

`clean/cenrep/pyproject.toml`:
```toml
[project]
name = "cenrep"
version = "0.1.0"
description = "Centromeric satellite, HOR and tandem-repeat annotation for allopolyploid complexes"
requires-python = ">=3.10"
dependencies = ["numpy", "pandas", "pysam", "pyhmmer>=0.12", "matplotlib"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

`clean/cenrep/environment.yml` (documentation of the env actually used; do NOT create or install it):
```yaml
# cenrep runs in the existing mamba env `synLTR`; this file only records what it needs.
name: cenrep
channels: [conda-forge, bioconda]
dependencies:
  - python=3.10
  - numpy
  - pandas
  - pysam
  - pyhmmer>=0.12
  - matplotlib
  - pytest
  - trf=4.10.0rc2
  - blast>=2.16
  - mafft=7.525
  - cd-hit=4.8.1
  - minimap2
  - samtools
```

`clean/cenrep/cenrep/__init__.py`:
```python
"""cenrep: centromeric satellite, HOR and tandem-repeat annotation."""
__version__ = "0.1.0"
```

`clean/cenrep/cenrep/__main__.py` (replaced in Task 11):
```python
raise SystemExit("cenrep CLI not built yet")
```

- [ ] **Step 2: Write `cenrep/params.py`**

```python
"""Tunable thresholds for every stage. Each field is exposed as a --flag by the CLI."""
from dataclasses import asdict, dataclass


@dataclass
class Params:
    trf_chunk: int = 5_000_000
    trf_overlap: int = 10_000
    min_period: int = 30
    max_period: int = 2000
    min_seed_span: int = 2000
    fam_ident: float = 0.90
    multi_ident: float = 0.80
    superfam_ident: float = 0.70
    min_fam_span: int = 20_000
    max_families: int = 40
    hmm_sample: int = 300
    hmm_sample_bp: int = 500_000
    min_hmm_seqs: int = 20
    nhmmer_evalue: float = 1e-5
    array_gap_factor: float = 2.0
    min_array_monomers: int = 3
    chrom_min_bp: int = 20_000
    cluster_gap: int = 500_000
    min_concentration: float = 0.5
    max_subtel_frac: float = 0.5
    subtel_min_bp: int = 1_000_000
    subtel_frac_len: float = 0.05
    max_gene_ratio: float = 0.5
    strong_prevalence: float = 0.5
    strong_cluster_bp: int = 100_000
    low_sat_bp: int = 50_000
    competing_ratio: float = 0.5
    fallback_min_bp: int = 5_000
    cen_assoc_frac: float = 0.5
    min_gap_n: int = 10
    lu_flank: int = 500_000
    lu_window: int = 1_000_000
    lu_window_overlap: int = 60_000
    lu_min_unit: int = 2_000
    lu_max_unit: int = 30_000
    lu_min_ident: float = 0.80
    hor_ident: float = 0.97
    hor_min_cov: float = 0.90
    hor_max_len: int = 200
    tandem_min_len: int = 500


HELP = {
    "trf_chunk": "TRF chunk size (bp)",
    "trf_overlap": "overlap between TRF chunks (bp)",
    "min_period": "min period of a family seed (bp)",
    "max_period": "max period of a family seed (bp)",
    "min_seed_span": "min TRF array span to seed a family (bp)",
    "fam_ident": "family identity, seeds of equal period",
    "multi_ident": "identity when one seed is a k-mer (k>=2) of the other",
    "superfam_ident": "superfamily identity (descriptive only)",
    "min_fam_span": "min union span in >=1 genome for a family to be annotated (bp)",
    "max_families": "max families annotated per complex",
    "hmm_sample": "max monomers sampled per family for MSA/HMM",
    "hmm_sample_bp": "max seed-array bp scanned per genome when sampling monomers",
    "min_hmm_seqs": "min sampled monomers for a profile HMM (else single-sequence HMM)",
    "nhmmer_evalue": "nhmmer E-value cutoff for monomer hits",
    "array_gap_factor": "max gap inside an array, in monomer lengths",
    "min_array_monomers": "min monomers per array",
    "chrom_min_bp": "min family bp on a chromosome / min main-cluster bp (bp)",
    "cluster_gap": "arrays closer than this form one cluster (bp)",
    "min_concentration": "min median share of a chromosome's family bp in its main cluster",
    "max_subtel_frac": "max share of family bp inside chromosome end zones",
    "subtel_min_bp": "end zone = max(this, subtel-frac-len x chromosome length) (bp)",
    "subtel_frac_len": "end zone as a fraction of chromosome length",
    "max_gene_ratio": "max gene density in main clusters relative to genome mean",
    "strong_prevalence": "min prevalence for 'strong' evidence",
    "strong_cluster_bp": "min median main-cluster bp for 'strong' evidence (bp)",
    "low_sat_bp": "centromere satellite bp below this -> low_satellite flag (bp)",
    "competing_ratio": "2nd cluster >= this x main cluster -> competing_cluster flag",
    "fallback_min_bp": "min cluster bp for a weak_satellite fallback call (bp)",
    "cen_assoc_frac": "min share of a non-primary family's bp inside centromeres -> cen_associated",
    "min_gap_n": "min N-run length counted as an assembly gap",
    "lu_flank": "long-unit scan flank around each centromere (bp)",
    "lu_window": "long-unit self-alignment window (bp)",
    "lu_window_overlap": "overlap between long-unit windows (bp)",
    "lu_min_unit": "min long-unit period (bp)",
    "lu_max_unit": "max long-unit period (bp)",
    "lu_min_ident": "min self-alignment identity for long-unit arrays",
    "hor_ident": "cd-hit-est identity for monomer variant classes",
    "hor_min_cov": "min model coverage for a monomer to enter the HOR search",
    "hor_max_len": "max HOR unit length in monomers (longer repeats are capped)",
    "tandem_min_len": "min length of other tandem arrays reported inside centromeres (bp)",
}


def to_dict(p: Params) -> dict:
    return asdict(p)
```

- [ ] **Step 3: Write `cenrep/util.py`**

```python
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
```

- [ ] **Step 4: Write `cenrep/genomes.py`**

```python
"""Complex-directory inputs: genomes, chromosome lengths, miniprot gene positions."""
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .util import CenrepError, log


@dataclass
class Genome:
    name: str      # file stem, e.g. Akamchatica_L
    fasta: Path
    chroms: dict   # chrom -> length, .fai order


def load_fai(fasta) -> dict:
    fai = Path(str(fasta) + ".fai")
    if not fai.exists():
        raise CenrepError(f"missing FASTA index {fai} (run: samtools faidx {fasta})")
    out = {}
    with open(fai) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or not p[1].isdigit():
                raise CenrepError(f"corrupt .fai line in {fai}: {line.strip()!r}")
            out[p[0]] = int(p[1])
    if not out:
        raise CenrepError(f"empty FASTA index {fai}")
    return out


def discover_genomes(complex_dir) -> list:
    complex_dir = Path(complex_dir)
    fas = sorted(complex_dir.glob("*.fa"))
    if len(fas) != 4:
        raise CenrepError(f"{complex_dir}: expected 4 *.fa genomes, found {len(fas)}")
    return [Genome(f.stem, f, load_fai(f)) for f in fas]


def find_gene_gff(complex_dir, genome: str):
    for gff in sorted(Path(complex_dir).glob(f"window_identity/*/{genome}_mod.gff")):
        ex = gff.parent / "run.exit"
        if ex.exists() and ex.read_text().strip() == "0":
            return gff
    return None


def load_genes(complex_dir, genome: str):
    """{chrom: sorted np.int64 array of mRNA midpoints (0-based)}, or None when no finished run exists."""
    complex_dir = Path(complex_dir)
    gff = find_gene_gff(complex_dir, genome)
    cmap = complex_dir / "window_identity" / "genomes" / "chr_map.tsv"
    if gff is None or not cmap.exists():
        log.warning(f"{genome}: no finished window_identity gene GFF; gene evidence off (no_genes)")
        return None
    back = {}
    with open(cmap) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[0] == genome:
                back[f"{p[0]}_{p[2]}"] = p[1]
    mids = {}
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 5 or p[2] != "mRNA" or p[0] not in back:
                continue
            mids.setdefault(back[p[0]], []).append((int(p[3]) - 1 + int(p[4])) // 2)
    return {c: np.sort(np.array(v, dtype=np.int64)) for c, v in mids.items()}


def count_in(mids, start: int, end: int) -> int:
    if mids is None or len(mids) == 0:
        return 0
    return int(np.searchsorted(mids, end) - np.searchsorted(mids, start))
```

- [ ] **Step 5: Write the synthetic genome builder `tests/synth.py`**

```python
"""Synthetic 4-genome complex with planted centromeric arrays, HORs, decoys, a gap and a long-unit array.

Each chromosome is L bp. Chromosomes 1-3 carry a centromere core starting at CEN_START; chromosome 4 has none.
Both ends of every chromosome carry a 90-bp subtelomeric decoy array (40 copies).
"""
from pathlib import Path

import numpy as np
import pysam

L = 300_000
CEN_START = 130_000
MON = 150


def rseq(rng, n: int) -> str:
    return "".join(rng.choice(list("ACGT"), size=n, p=[0.31, 0.19, 0.19, 0.31]))


def mutate(rng, s: str, rate: float) -> str:
    a = list(s)
    for i in np.where(rng.random(len(a)) < rate)[0]:
        a[i] = rng.choice([b for b in "ACGT" if b != a[i]])
    return "".join(a)


def rc(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def build_complex(outdir: Path, seed: int = 11) -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    shared = np.random.default_rng(seed)
    cen, dec = rseq(shared, MON), rseq(shared, 90)
    hor4 = [mutate(shared, cen, 0.08) for _ in range(4)]
    dup3 = [mutate(shared, cen, 0.08) for _ in range(3)]
    u3k = rseq(shared, 3000)
    truth = {"cen": cen, "dec": dec, "genomes": {}}
    for gi, gname in enumerate(["Pa_A", "Pb_B", "Tx_A", "Tx_B"]):
        rng = np.random.default_rng(seed + 100 + gi)

        def plain(n):
            return "".join(mutate(rng, cen, 0.08) for _ in range(n))

        def decoy():
            return "".join(mutate(rng, dec, 0.02) for _ in range(40))

        chroms, t = {}, {}
        for ci in range(1, 5):
            name = f"{gname}{ci:02d}"
            if ci == 1:
                hor = "".join(mutate(rng, v, 0.005) for _ in range(15) for v in hor4)
                parts = [("sat", plain(60)), ("hor", hor), ("sat", plain(60)), ("inv", rc(plain(60))), ("sat", plain(20))]
            elif ci == 2:
                d1 = "".join(mutate(rng, v, 0.005) for v in dup3)
                d2 = "".join(mutate(rng, v, 0.005) for v in dup3)
                parts = [("sat", plain(100)), ("dup", d1), ("sat", plain(60)), ("dup", d2), ("sat", plain(100))]
            elif ci == 3:
                lu = "".join(mutate(rng, u3k, 0.03) for _ in range(6))
                parts = [("sat", plain(120)), ("micro", "AT" * 300), ("sat", plain(120)),
                         ("gap", "N" * 500), ("lu", lu), ("sat", plain(60))]
            else:
                parts = []
            core = "".join(s for _, s in parts)
            left = decoy()
            left += rseq(rng, CEN_START - len(left))
            right = decoy()
            seq = left + core + rseq(rng, L - CEN_START - len(core) - len(right)) + right
            assert len(seq) == L
            feats, pos = {}, CEN_START
            for kind, s in parts:
                feats.setdefault(kind, []).append((pos, pos + len(s)))
                pos += len(s)
            t[name] = {"cen": (CEN_START, CEN_START + len(core)) if parts else None, **feats}
            chroms[name] = seq
        fa = outdir / f"{gname}.fa"
        with open(fa, "w") as f:
            for n, s in chroms.items():
                f.write(f">{n}\n")
                for i in range(0, len(s), 60):
                    f.write(s[i:i + 60] + "\n")
        pysam.faidx(str(fa))
        truth["genomes"][gname] = t
    return truth
```

- [ ] **Step 6: Write `tests/conftest.py`**

```python
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
```

- [ ] **Step 7: Write the failing tests**

`tests/test_util.py`:
```python
import pytest

from cenrep.util import (CenrepError, canonical_rotation, clipped_bp, merge_intervals, parse_fasta_text,
                         revcomp, run, union_bp)


def test_canonical_rotation_is_rotation_invariant():
    s = "GATTACAGATTTCC"
    rots = {canonical_rotation(s[i:] + s[:i]) for i in range(len(s))}
    assert len(rots) == 1 and len(next(iter(rots))) == len(s)


def test_merge_intervals_gap():
    assert merge_intervals([(10, 20), (0, 5), (22, 30)], gap=2) == [[0, 5], [10, 30]]
    assert merge_intervals([(0, 5), (5, 8)]) == [[0, 8]]


def test_union_and_clipped_bp():
    iv = [(0, 10), (5, 15), (20, 30)]
    assert union_bp(iv) == 25
    assert clipped_bp(iv, 8, 25) == 12


def test_revcomp_and_fasta_parse():
    assert revcomp("AACGTn") == "nACGTT"
    assert parse_fasta_text(">a x\nAC\nGT\n>b\nTT\n") == [("a", "ACGT"), ("b", "TT")]


def test_run_raises_on_failure():
    with pytest.raises(CenrepError):
        run(["false"])
```

`tests/test_genomes.py`:
```python
import pytest

from cenrep.genomes import count_in, discover_genomes, load_fai, load_genes
from cenrep.util import CenrepError


def test_load_fai_missing(tmp_path):
    fa = tmp_path / "x.fa"
    fa.write_text(">c\nACGT\n")
    with pytest.raises(CenrepError, match="missing FASTA index"):
        load_fai(fa)


def test_discover_requires_four(tmp_path):
    for n in "ab":
        (tmp_path / f"{n}.fa").write_text(">c\nA\n")
    with pytest.raises(CenrepError, match="expected 4"):
        discover_genomes(tmp_path)


def test_load_genes_maps_renamed_chroms(tmp_path):
    run_dir = tmp_path / "window_identity" / "1_x"
    run_dir.mkdir(parents=True)
    (tmp_path / "window_identity" / "genomes").mkdir()
    (tmp_path / "window_identity" / "genomes" / "chr_map.tsv").write_text(
        "Sp_A\tSp_A01\tchr1\nSp_A\tSp_A02\tchr2\nSp_B\tSp_B01\tchr1\n")
    (run_dir / "Sp_A_mod.gff").write_text(
        "Sp_A_chr1\tminiprot\tmRNA\t101\t200\t.\t+\t.\tID=a\n"
        "Sp_A_chr1\tminiprot\tCDS\t101\t200\t.\t+\t0\tParent=a\n"
        "Sp_A_chr2\tminiprot\tmRNA\t1\t10\t.\t-\t.\tID=b\n")
    assert load_genes(tmp_path, "Sp_A") is None          # run not finished (no run.exit)
    (run_dir / "run.exit").write_text("0\n")
    g = load_genes(tmp_path, "Sp_A")
    assert set(g) == {"Sp_A01", "Sp_A02"} and g["Sp_A01"].tolist() == [150]
    assert count_in(g["Sp_A01"], 100, 151) == 1 and count_in(g["Sp_A01"], 151, 300) == 0
```

`tests/test_synth.py`:
```python
import pysam


def test_synthetic_complex_layout(synth_complex):
    cdir, truth = synth_complex
    fas = sorted(cdir.glob("*.fa"))
    assert len(fas) == 4 and all((f.parent / (f.name + ".fai")).exists() for f in fas)
    t = truth["genomes"]["Pa_A"]
    assert t["Pa_A01"]["cen"] == (130_000, 130_000 + 260 * 150)
    assert t["Pa_A04"]["cen"] is None
    s, e = t["Pa_A03"]["gap"][0]
    with pysam.FastaFile(str(cdir / "Pa_A.fa")) as fa:
        assert set(fa.fetch("Pa_A03", s, e)) == {"N"}
```

- [ ] **Step 8: Run the tests**

Run: `cd /data2/chris/tetraploid/allo10/clean/cenrep && PATH=$ENV/bin:$PATH python -m pytest -q tests/test_util.py tests/test_genomes.py tests/test_synth.py`
Expected: all PASS. The modules were written in Steps 2–5, so if anything fails, fix the code, not the test.

---

### Task 2: S1 — TRF discovery (`cenrep/trf.py`)

**Files:**
- Create: `clean/cenrep/cenrep/trf.py`
- Modify: `clean/cenrep/tests/conftest.py` (add the `synth_trf` fixture)
- Test: `clean/cenrep/tests/test_trf.py`

**Interfaces:**
- Consumes: `Genome`, `Params.trf_chunk/trf_overlap`, `util.run/md5_file/log`
- Produces:
  - `TRF_ARGS`, `COLS` (list of column names), `chunks(chroms, size, overlap)`
  - `parse_ngs(text, chrom, offset) -> list[tuple]`, `dedupe(rows) -> list[tuple]`
  - `trf_genomes(genomes, cache_dir, threads, p) -> dict[genome_name, Path]` (cached `<md5>.trf.tsv.gz`), `load_trf(path) -> pd.DataFrame`
  - Row tuple order = `COLS`: `chrom, start, end, period, copies, consensus_size, pct_match, pct_indel, score, entropy, consensus`

- [ ] **Step 1: Write the failing tests** — `tests/test_trf.py`

```python
import pytest

from cenrep.genomes import discover_genomes
from cenrep.trf import chunks, dedupe, load_trf, parse_ngs, trf_genomes
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
    hit = df[(df.chrom == "Pa_A01") & df.period.between(145, 155) & (df.start < e) & (df.end > s)]
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
```

Append to `tests/conftest.py`:
```python
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest -q tests/test_trf.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'cenrep.trf'`

- [ ] **Step 3: Write `cenrep/trf.py`**

```python
"""S1: TRF tandem-repeat discovery on chunked chromosomes, cached by genome md5."""
import gzip
import os
import shutil
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import pysam

from .util import log, md5_file, run

# Benson's recommended weights (match 2, mismatch 7, indel 7, PM 80, PI 10), minscore 50,
# max period 2000 (TRF maximum), -l 6: chunks are <= 5 Mb so no array is longer.
TRF_ARGS = ["2", "7", "7", "80", "10", "50", "2000", "-l", "6", "-h", "-ngs"]
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


def trf_genomes(genomes, cache_dir, threads: int, p) -> dict:
    """Run TRF on every genome not cached yet (md5-keyed). Returns {genome name: cached .tsv.gz path}."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    out, todo = {}, {}
    for g in genomes:
        path = cache_dir / f"{md5_file(g.fasta)}.trf.tsv.gz"
        out[g.name] = path
        if path.exists():
            log.info(f"trf: {g.name} cached ({path.name})")
        else:
            todo.setdefault(path, g)          # md5-identical genomes share one run
    if not todo:
        return out
    tmp = tempfile.mkdtemp(prefix="trf_tmp_", dir=cache_dir)
    jobs = [(g.fasta, c, s, e, tmp) for g in todo.values()
            for c, s, e in chunks(g.chroms, p.trf_chunk, p.trf_overlap)]
    log.info(f"trf: {len(todo)} genome(s), {len(jobs)} chunks, {threads} workers")
    rows = {str(g.fasta): [] for g in todo.values()}
    with ProcessPoolExecutor(threads) as ex:
        for i, (fa, r) in enumerate(ex.map(_run_chunk, jobs, chunksize=1), 1):
            rows[fa].extend(r)
            log.debug(f"trf: chunk {i}/{len(jobs)}")
    for path, g in todo.items():
        kept = dedupe(rows[str(g.fasta)])
        part = Path(str(path) + ".part")
        with gzip.open(part, "wt") as f:
            f.write("\t".join(COLS) + "\n")
            for r in kept:
                f.write("\t".join(map(str, r)) + "\n")
        part.rename(path)
        log.info(f"trf: {g.name} {len(kept)} TRF rows")
    shutil.rmtree(tmp, ignore_errors=True)
    return out


def load_trf(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"chrom": str, "consensus": str})
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest -q tests/test_trf.py`
Expected: PASS (5 tests).

---

### Task 3: S2 — satellite families (`cenrep/families.py`)

**Files:**
- Create: `clean/cenrep/cenrep/families.py`
- Modify: `clean/cenrep/tests/conftest.py` (add the `synth_families` fixture)
- Test: `clean/cenrep/tests/test_families.py`

**Interfaces:**
- Consumes: TRF tables (`COLS`), `Genome`, `Params`, `util.*`
- Produces:
  - `select_seeds(trf: dict[str, DataFrame], p) -> DataFrame` (TRF cols + `genome, span, seed_id`)
  - `seed_match(slen, llen, hsps, t_same, t_multi) -> bool`, `blast_pairs(queries, subjects, workdir, threads) -> dict`
  - `cluster_seeds(seqs, order, workdir, threads, t_same, t_multi) -> dict[id, rep_id]`
  - `msa_consensus(rows, max_gap=0.5) -> str`, `single_hmm(name, seq)`, `msa_hmm(name, rows)`
  - `build_families(trf, genomes, workdir, threads, p) -> DataFrame`, which writes:
    - `workdir/families.tsv` with columns `FAM_COLS = [fam_id, rep_seed, base_seed, consensus_len, n_seeds, max_genome_span, genome_spans, hmm_mode, n_hmm_seqs, superfamily]` (superfamily = rep `fam_id`)
    - `workdir/families.fa` (names `famNNN`), `workdir/families.hmm` (HMM names `famNNN`), `workdir/seeds.tsv.gz`

- [ ] **Step 1: Write the failing tests** — `tests/test_families.py`

```python
import numpy as np
import pytest

from cenrep.families import cluster_seeds, msa_consensus, seed_match
from conftest import have
from synth import mutate, rseq


def test_seed_match_rotation_same_period():
    assert seed_match(150, 150, [(97.0, 1, 150, 40, 189)], 0.9, 0.8)


def test_seed_match_rejects_low_identity_same_period():
    assert not seed_match(150, 150, [(88.0, 1, 150, 40, 189)], 0.9, 0.8)


def test_seed_match_multimer_needs_tiling_and_integer_ratio():
    tiled = [(85.0, 1, 150, 1 + 150 * i, 150 * (i + 1)) for i in range(6)]   # 6 copies on a 900-bp dimer
    assert seed_match(150, 450, tiled, 0.9, 0.8)
    assert not seed_match(150, 450, tiled[:1], 0.9, 0.8)      # covers 1/3 of the circle
    assert not seed_match(150, 400, tiled, 0.9, 0.8)          # 400/150 is not an integer ratio


def test_msa_consensus_drops_gappy_columns():
    assert msa_consensus(["AC-GT", "ACTGT", "AC-GA"]) == "ACGT"


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_cluster_seeds_rotation_multimer_unrelated(tmp_path):
    rng = np.random.default_rng(5)
    a, b = rseq(rng, 150), rseq(rng, 150)
    seqs = {"A": a,
            "Arot": mutate(rng, a[40:] + a[:40], 0.02),
            "A3": "".join(mutate(rng, a, 0.05) for _ in range(3)),
            "B": b,
            "Asub": mutate(rng, a, 0.14)}
    assign = cluster_seeds(seqs, ["A", "Arot", "A3", "B", "Asub"], tmp_path, 2, 0.90, 0.80)
    assert assign["Arot"] == "A" and assign["A3"] == "A"
    assert assign["B"] == "B" and assign["Asub"] == "Asub"


@pytest.mark.skipif(not have("trf", "blastn", "mafft"), reason="tools missing")
def test_build_families_on_synthetic(synth_families):
    fam, workdir = synth_families
    near = lambda n: fam[(fam.consensus_len - n).abs() <= 3]
    assert len(near(150)) == 1 and near(150).hmm_mode.iloc[0] == "profile"   # CEN + its HOR multimer
    assert len(near(90)) == 1                                                  # subtelomeric decoy
    assert len(near(600)) == 0                                                 # multimer merged into CEN
    assert (workdir / "families.hmm").stat().st_size > 0
```

Append to `tests/conftest.py`:
```python
@pytest.fixture(scope="session")
def synth_families(synth_complex, synth_trf, test_params, tmp_path_factory):
    from cenrep.families import build_families
    from cenrep.genomes import discover_genomes
    wd = tmp_path_factory.mktemp("families")
    fam = build_families(synth_trf, discover_genomes(synth_complex[0]), wd, 4, test_params)
    return fam, wd
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest -q tests/test_families.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'cenrep.families'`

- [ ] **Step 3: Write `cenrep/families.py`**

```python
"""S2: satellite family library from TRF seeds (rotation-invariant, multimer-aware)."""
import gzip
import random
from pathlib import Path

import numpy as np
import pandas as pd
import pyhmmer
import pysam
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextMSA, TextSequence

from .util import canonical_rotation, log, merge_intervals, parse_fasta_text, revcomp, run, union_bp, write_fasta

ABC = Alphabet.dna()
BATCH = 200
BLAST_FMT = "6 qseqid sseqid pident length qstart qend sstart send qlen slen bitscore"
FAM_COLS = ["fam_id", "rep_seed", "base_seed", "consensus_len", "n_seeds", "max_genome_span",
            "genome_spans", "hmm_mode", "n_hmm_seqs", "superfamily"]


def select_seeds(trf: dict, p) -> pd.DataFrame:
    parts = []
    for g, df in trf.items():
        d = df[(df.period >= p.min_period) & (df.period <= p.max_period)
               & ((df.end - df.start) >= p.min_seed_span)].copy()
        d["genome"] = g
        parts.append(d)
    s = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if s.empty:
        return s
    s["span"] = s.end - s.start
    s = s.sort_values(["span", "genome", "chrom", "start"],
                      ascending=[False, True, True, True]).reset_index(drop=True)
    s["seed_id"] = [f"s{i}" for i in range(len(s))]
    return s


def seed_match(slen: int, llen: int, hsps, t_same: float, t_multi: float) -> bool:
    """Shorter seed (length slen, BLAST query) vs the longer seed's dimer (2*llen, subject).
    hsps: [(pident, qstart, qend, sstart, send), ...], BLAST 1-based coordinates."""
    ratio = llen / slen
    k = max(1, round(ratio))
    if abs(ratio - k) > 0.05 * k:
        return False
    t = 100 * (t_same if k == 1 else t_multi)
    good = [h for h in hsps if h[0] >= t]
    if not any(h[2] - h[1] + 1 >= 0.8 * slen for h in good):
        return False
    cov = np.zeros(llen, bool)
    for h in good:
        a, b = min(h[3], h[4]) - 1, max(h[3], h[4])
        cov[np.arange(a, b) % llen] = True          # circular coverage of the longer monomer
    return cov.mean() >= 0.8


def blast_pairs(queries: dict, subjects: dict, workdir: Path, threads: int) -> dict:
    """BLAST query monomers against subject dimers. Returns {(qid, sid): [hsp, ...]}."""
    workdir.mkdir(parents=True, exist_ok=True)
    q, s, db, out = workdir / "q.fa", workdir / "s.fa", workdir / "sdb", workdir / "hits.tsv"
    write_fasta(q, queries.items())
    write_fasta(s, ((k, v + v) for k, v in subjects.items()))
    run(["makeblastdb", "-in", s, "-dbtype", "nucl", "-out", db])
    run(["blastn", "-task", "blastn", "-query", q, "-db", db, "-word_size", 11, "-evalue", "1e-5",
         "-dust", "no", "-soft_masking", "false", "-max_target_seqs", 100000,
         "-num_threads", threads, "-outfmt", BLAST_FMT, "-out", out])
    pairs = {}
    with open(out) as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            pairs.setdefault((x[0], x[1]), []).append(
                (float(x[2]), int(x[4]), int(x[5]), int(x[6]), int(x[7])))
    return pairs


def _related(a, b, seqs, pairs, t_same, t_multi) -> bool:
    la, lb = len(seqs[a]), len(seqs[b])
    if la <= lb:
        return seed_match(la, lb, pairs.get((a, b), []), t_same, t_multi)
    return seed_match(lb, la, pairs.get((b, a), []), t_same, t_multi)


def cluster_seeds(seqs: dict, order: list, workdir: Path, threads: int, t_same: float, t_multi: float) -> dict:
    """Greedy centroid clustering in priority order (no chaining). Returns {id: representative id}.
    Works in batches of BATCH candidate representatives so BLAST output stays bounded."""
    assign, rest, rnd = {}, list(order), 0
    while rest:
        rnd += 1
        cand, others = rest[:BATCH], rest[BATCH:]
        sub = {i: seqs[i] for i in cand}
        pairs = blast_pairs(sub, sub, workdir / f"r{rnd}c", threads)
        reps = []
        for c in cand:
            hit = next((r for r in reps if _related(c, r, seqs, pairs, t_same, t_multi)), None)
            assign[c] = c if hit is None else hit
            if hit is None:
                reps.append(c)
        if others:
            o, r = {i: seqs[i] for i in others}, {i: seqs[i] for i in reps}
            pairs = blast_pairs(o, r, workdir / f"r{rnd}o1", threads)
            for k, v in blast_pairs(r, o, workdir / f"r{rnd}o2", threads).items():
                pairs.setdefault(k, []).extend(v)
            for x in others:
                hit = next((rp for rp in reps if _related(x, rp, seqs, pairs, t_same, t_multi)), None)
                if hit is not None:
                    assign[x] = hit
        rest = [x for x in others if x not in assign]
        log.debug(f"families: round {rnd}: {len(reps)} new representatives, {len(rest)} seeds left")
    return assign


def msa_consensus(rows: list, max_gap: float = 0.5) -> str:
    arr = np.array([list(r) for r in rows])
    out = []
    for col in arr.T:
        if (col == "-").mean() >= max_gap:
            continue
        vals, cnt = np.unique(col[col != "-"], return_counts=True)
        out.append(vals[np.argmax(cnt)])
    return "".join(out)


def single_hmm(name: str, seq: str):
    bld, bg = pyhmmer.plan7.Builder(ABC), pyhmmer.plan7.Background(ABC)
    hmm, _, _ = bld.build(TextSequence(name=name, sequence=seq).digitize(ABC), bg)
    hmm.name = name
    return hmm


def msa_hmm(name: str, rows: list):
    msa = TextMSA(name=name, sequences=[TextSequence(name=f"m{i}", sequence=r) for i, r in enumerate(rows)])
    bld, bg = pyhmmer.plan7.Builder(ABC), pyhmmer.plan7.Background(ABC)
    hmm, _, _ = bld.build_msa(msa.digitize(ABC), bg)
    hmm.name = name
    return hmm


def mafft(seqs: list, workdir: Path, threads: int) -> list:
    workdir.mkdir(parents=True, exist_ok=True)
    inp = workdir / "mafft_in.fa"
    write_fasta(inp, ((f"m{i}", s) for i, s in enumerate(seqs)))
    out = run(["mafft", "--auto", "--thread", threads, "--quiet", inp]).stdout
    return [s.upper() for _, s in parse_fasta_text(out)]


def seed_intervals(members: pd.DataFrame, p) -> list:
    """Largest merged seed intervals per genome, up to p.hmm_sample_bp per genome."""
    out = []
    for g, d in members.groupby("genome"):
        iv = [(c, s, e) for c, x in d.groupby("chrom") for s, e in merge_intervals(zip(x.start, x.end))]
        iv.sort(key=lambda t: t[1] - t[2])                  # longest first
        budget = p.hmm_sample_bp
        for c, s, e in iv:
            if budget <= 0:
                break
            e = min(e, s + budget)
            out.append((g, c, s, e))
            budget -= e - s
    return out


def sample_monomers(hmm, intervals: list, fastas: dict, p, threads: int, rng: random.Random) -> list:
    """nhmmer hmm vs the intervals; up to p.hmm_sample monomers with >= 90% model coverage,
    in consensus orientation."""
    raw = []
    for g in dict.fromkeys(i[0] for i in intervals):
        with pysam.FastaFile(str(fastas[g])) as fa:
            raw += [fa.fetch(c, s, e).upper() for gg, c, s, e in intervals if gg == g]
    if not raw:
        return []
    blk = DigitalSequenceBlock(ABC, [TextSequence(name=f"i{i}", sequence=x).digitize(ABC)
                                     for i, x in enumerate(raw)])
    out = []
    for hits in pyhmmer.hmmer.nhmmer([hmm], blk, cpus=threads, E=p.nhmmer_evalue):
        for h in hits:
            if h.evalue > p.nhmmer_evalue:
                continue
            a = h.best_domain.alignment
            if a.hmm_to - a.hmm_from + 1 < 0.9 * hmm.M:
                continue
            lo, hi = sorted((a.target_from, a.target_to))
            sub = raw[int(h.name[1:])][lo - 1:hi]
            out.append(sub if h.best_domain.strand == "+" else revcomp(sub))
    rng.shuffle(out)
    return out[:p.hmm_sample]


def build_hmm(fam_id: str, rep_seq: str, intervals, fastas, workdir: Path, threads: int, p, rng):
    """Two-pass family model. Returns (hmm, consensus, mode, n_sampled)."""
    c0 = canonical_rotation(rep_seq)
    h = single_hmm(fam_id, c0)
    s1 = sample_monomers(h, intervals, fastas, p, threads, rng)
    if len(s1) < p.min_hmm_seqs:
        log.warning(f"families: {fam_id}: only {len(s1)} monomers sampled; single-sequence HMM")
        return h, c0, "single", len(s1)
    c1 = canonical_rotation(msa_consensus(mafft(s1, workdir, threads)))
    h = single_hmm(fam_id, c1)
    s2 = sample_monomers(h, intervals, fastas, p, threads, rng)
    if len(s2) < p.min_hmm_seqs:
        log.warning(f"families: {fam_id}: only {len(s2)} monomers in pass 2; single-sequence HMM")
        return h, c1, "single", len(s2)
    h = msa_hmm(fam_id, mafft(s2, workdir, threads))
    return h, h.consensus.upper(), "profile", len(s2)


def build_families(trf: dict, genomes: list, workdir: Path, threads: int, p) -> pd.DataFrame:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    seeds = select_seeds(trf, p)
    fams, hmms = [], []
    if seeds.empty:
        log.warning("families: no TRF seeds; no satellite families")
    else:
        seqs = dict(zip(seeds.seed_id, seeds.consensus))
        assign = cluster_seeds(seqs, list(seeds.seed_id), workdir / "blast", threads, p.fam_ident, p.multi_ident)
        seeds["rep"] = seeds.seed_id.map(assign)
        seeds["fam_id"] = "."
        spans = {}
        for (rep, g), d in seeds.groupby(["rep", "genome"]):
            spans.setdefault(rep, {})[g] = sum(union_bp(zip(x.start, x.end)) for _, x in d.groupby("chrom"))
        top = sorted(((max(v.values()), r) for r, v in spans.items() if max(v.values()) >= p.min_fam_span),
                     reverse=True)[:p.max_families]
        log.info(f"families: {len(seeds)} seeds -> {len(set(assign.values()))} clusters; "
                 f"{len(top)} carried forward")
        fastas = {g.name: g.fasta for g in genomes}
        rng = random.Random(1)
        for i, (mx, rep) in enumerate(top, 1):
            fid = f"fam{i:03d}"
            mem = seeds[seeds.rep == rep]
            lens = mem.consensus.str.len()
            base = mem[lens <= 1.05 * lens.min()].iloc[0]      # largest-span shortest member = base monomer
            hmm, cons, mode, n = build_hmm(fid, base.consensus, seed_intervals(mem, p), fastas,
                                           workdir / fid, threads, p, rng)
            seeds.loc[mem.index, "fam_id"] = fid
            fams.append(dict(fam_id=fid, rep_seed=rep, base_seed=base.seed_id, consensus_len=hmm.M,
                             consensus=cons, n_seeds=len(mem), max_genome_span=mx,
                             genome_spans=",".join(f"{g}:{v}" for g, v in sorted(spans[rep].items())),
                             hmm_mode=mode, n_hmm_seqs=n))
            hmms.append(hmm)
            log.info(f"families: {fid} {hmm.M} bp, {len(mem)} seeds, max span {mx} bp, hmm={mode} (n={n})")
        if fams:
            cs = {f["fam_id"]: f["consensus"] for f in fams}
            sf = cluster_seeds(cs, list(cs), workdir / "blast_sf", threads, p.superfam_ident, p.superfam_ident)
            for f in fams:
                f["superfamily"] = sf[f["fam_id"]]
        seeds.to_csv(workdir / "seeds.tsv.gz", sep="\t", index=False)
    table = pd.DataFrame(fams, columns=FAM_COLS)
    table.to_csv(workdir / "families.tsv", sep="\t", index=False)
    write_fasta(workdir / "families.fa", ((f["fam_id"], f["consensus"]) for f in fams))
    with open(workdir / "families.hmm", "wb") as fh:
        for h in hmms:
            h.write(fh)
    return table
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest -q tests/test_families.py`
Expected: PASS (6 tests).
- If `test_build_families_on_synthetic` finds a ~600 bp family, print `seeds.tsv.gz` rows of that family and check the multimer rule (`multi_ident`, circular tiling) before touching thresholds.
- If `sseqid` in the BLAST output is not the bare seed ID, add `-parse_seqids` to `makeblastdb`.

---

### Task 4: S3 — genome-wide monomer calls and arrays (`cenrep/annotate.py`)

**Files:**
- Create: `clean/cenrep/cenrep/annotate.py`
- Modify: `clean/cenrep/tests/conftest.py` (add `synth_annot`)
- Test: `clean/cenrep/tests/test_annotate.py`

**Interfaces:**
- Consumes: HMMs (names `famNNN`), `fam_len: dict[fam_id, int]`, `Genome`, `Params`
- Produces:
  - `MON_COLS = [chrom, start, end, strand, family, bitscore, evalue, hmm_from, hmm_to, model_len, identity]`
  - `ARR_COLS = [array_id, chrom, start, end, strand, family, n_monomers, mean_identity, monomer_len]`
  - `search_chrom(hmms, chrom, seq, threads, evalue) -> list[tuple]`, `resolve_overlaps(df, max_frac=0.5) -> DataFrame`
  - `chain_arrays(mon, fam_len, p) -> (mon_with_array_id, arrays_df)`
  - `annotate_genome(genome, hmms, fam_len, threads, p, out_mon, out_arr) -> (mon, arrays)`. It writes `<G>.monomers.tsv.gz`, with columns `mon_idx` + `MON_COLS` + `array_id`, where `mon_idx == row number` and rows are sorted by (.fai order, start). It also writes `<G>.arrays.tsv` (`ARR_COLS`).

- [ ] **Step 1: Write the failing tests** — `tests/test_annotate.py`

```python
import pandas as pd
import pytest

from cenrep.annotate import MON_COLS, chain_arrays, resolve_overlaps
from conftest import have


def _mon(rows):
    return pd.DataFrame(rows, columns=MON_COLS)


def test_resolve_overlaps_keeps_best_and_adjacent():
    rows = [("c", 0, 150, "+", "f1", 100.0, 1e-20, 1, 150, 150, 0.9),
            ("c", 2, 150, "+", "f2", 80.0, 1e-15, 1, 148, 150, 0.85),   # same monomer, other family: dropped
            ("c", 145, 295, "+", "f1", 90.0, 1e-18, 1, 150, 150, 0.9)]   # 5-bp overlap only: kept
    out = resolve_overlaps(_mon(rows))
    assert out.family.tolist() == ["f1", "f1"] and out.start.tolist() == [0, 145]


def test_chain_arrays_splits_on_gap_and_strand(test_params):
    rows = [("c", i * 150, (i + 1) * 150, "+", "f", 90.0, 1e-20, 1, 150, 150, 0.9) for i in range(5)]
    rows += [("c", 1150, 1300, "+", "f", 90.0, 1e-20, 1, 150, 150, 0.9)]          # gap 400 > 300: singleton
    rows += [("c", 2000 + i * 150, 2150 + i * 150, "-", "f", 90.0, 1e-20, 1, 150, 150, 0.8) for i in range(3)]
    mon, arr = chain_arrays(_mon(rows), {"f": 150}, test_params)
    assert sorted(arr.n_monomers.tolist()) == [3, 5]
    assert set(arr.strand) == {"+", "-"} and len(mon) == 8


@pytest.mark.skipif(not have("trf", "blastn", "mafft"), reason="tools missing")
def test_annotate_synthetic_centromere(synth_annot, synth_complex, synth_families):
    mon, arr = synth_annot
    fam, _ = synth_families
    cen_fam = fam[(fam.consensus_len - 150).abs() <= 3].fam_id.iloc[0]
    s, e = synth_complex[1]["genomes"]["Pa_A"]["Pa_A01"]["cen"]
    a = arr[(arr.chrom == "Pa_A01") & (arr.family == cen_fam)]
    assert abs(a.start.min() - s) <= 150 and abs(a.end.max() - e) <= 150
    assert (a.strand == "-").any()                                    # the inverted block
    assert len(mon[(mon.chrom == "Pa_A01") & (mon.family == cen_fam)]) >= 250
    assert mon.mon_idx.tolist() == list(range(len(mon)))
```

Append to `tests/conftest.py`:
```python
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest -q tests/test_annotate.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'cenrep.annotate'`

- [ ] **Step 3: Write `cenrep/annotate.py`**

```python
"""S3: genome-wide nhmmer monomer calls (one chromosome at a time), overlap resolution, array chaining."""
import numpy as np
import pandas as pd
import pyhmmer
import pysam
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence

from .util import log

ABC = Alphabet.dna()
MON_COLS = ["chrom", "start", "end", "strand", "family", "bitscore", "evalue",
            "hmm_from", "hmm_to", "model_len", "identity"]
ARR_COLS = ["array_id", "chrom", "start", "end", "strand", "family", "n_monomers", "mean_identity", "monomer_len"]


def _identity(a) -> float:
    """Exact matches / alignment columns (HMM consensus vs target, case-insensitive)."""
    h = np.frombuffer(a.hmm_sequence.upper().encode(), np.uint8)
    t = np.frombuffer(a.target_sequence.upper().encode(), np.uint8)
    return float((h == t).mean())


def search_chrom(hmms, chrom: str, seq: str, threads: int, evalue: float) -> list:
    if not hmms:
        return []
    blk = DigitalSequenceBlock(ABC, [TextSequence(name=chrom, sequence=seq.upper()).digitize(ABC)])
    mlen = {h.name: h.M for h in hmms}
    rows = []
    for hits in pyhmmer.hmmer.nhmmer(hmms, blk, cpus=threads, E=evalue):
        fam = hits.query.name
        for h in hits:
            if h.evalue > evalue:
                continue
            d = h.best_domain
            a = d.alignment
            lo, hi = sorted((a.target_from, a.target_to))
            rows.append((chrom, lo - 1, hi, d.strand, fam, float(h.score), float(h.evalue),
                         a.hmm_from, a.hmm_to, mlen[fam], _identity(a)))
    return rows


def resolve_overlaps(df: pd.DataFrame, max_frac: float = 0.5) -> pd.DataFrame:
    """Hits overlapping by > max_frac of the shorter compete; the higher bitscore wins (greedy per component)."""
    if df.empty:
        return df
    df = df.sort_values(["start", "end"]).reset_index(drop=True)
    s, e, sc = df.start.tolist(), df.end.tolist(), df.bitscore.tolist()
    n = len(df)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def clash(i, j):
        return min(e[i], e[j]) - max(s[i], s[j]) > max_frac * min(e[i] - s[i], e[j] - s[j])

    active = []
    for i in range(n):
        active = [j for j in active if e[j] > s[i]]
        for j in active:
            if clash(i, j):
                parent[find(i)] = find(j)
        active.append(i)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    keep = np.zeros(n, bool)
    for m in groups.values():
        acc = []
        for i in sorted(m, key=lambda k: -sc[k]):
            if not any(clash(i, j) for j in acc):
                acc.append(i)
        keep[acc] = True
    return df[keep].reset_index(drop=True)


def chain_arrays(mon: pd.DataFrame, fam_len: dict, p):
    """Monomers of one family on one strand with gap <= array_gap_factor x monomer length form an array;
    arrays need >= min_array_monomers. Returns (monomers in arrays with array_id, arrays)."""
    if mon.empty:
        return mon.assign(array_id=pd.Series(dtype="int64")), pd.DataFrame(columns=ARR_COLS)
    mon = mon.sort_values(["chrom", "family", "strand", "start"]).reset_index(drop=True)
    ids = np.full(len(mon), -1, np.int64)
    rows, aid = [], 0
    ident = mon.identity.to_numpy()
    for (c, f, st), g in mon.groupby(["chrom", "family", "strand"], sort=False):
        idx = g.index.to_numpy()
        s, e = g.start.to_numpy(), g.end.to_numpy()
        maxgap = p.array_gap_factor * fam_len[f]
        brk = np.where(s[1:] - np.maximum.accumulate(e)[:-1] > maxgap)[0] + 1
        for block in np.split(np.arange(len(idx)), brk):
            if len(block) < p.min_array_monomers:
                continue
            ids[idx[block]] = aid
            rows.append((aid, c, int(s[block].min()), int(e[block].max()), st, f, len(block),
                         float(ident[idx[block]].mean()), int(fam_len[f])))
            aid += 1
    mon["array_id"] = ids
    return mon[mon.array_id >= 0].copy(), pd.DataFrame(rows, columns=ARR_COLS)


def annotate_genome(genome, hmms, fam_len: dict, threads: int, p, out_mon, out_arr):
    parts = []
    with pysam.FastaFile(str(genome.fasta)) as fa:
        for c in genome.chroms:
            d = resolve_overlaps(pd.DataFrame(search_chrom(hmms, c, fa.fetch(c), threads, p.nhmmer_evalue),
                                              columns=MON_COLS))
            parts.append(d)
            log.debug(f"annotate: {genome.name} {c}: {len(d)} monomer hits")
    mon = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=MON_COLS)
    mon, arrays = chain_arrays(mon, fam_len, p)
    order = {c: i for i, c in enumerate(genome.chroms)}
    mon = mon.assign(ci=mon.chrom.map(order)).sort_values(["ci", "start", "end"]).drop(columns="ci")
    mon = mon.reset_index(drop=True)
    mon.insert(0, "mon_idx", np.arange(len(mon), dtype=np.int64))
    arrays = arrays.assign(ci=arrays.chrom.map(order)).sort_values(["ci", "start"]).drop(columns="ci")
    mon.to_csv(out_mon, sep="\t", index=False)
    arrays.to_csv(out_arr, sep="\t", index=False)
    log.info(f"annotate: {genome.name}: {len(mon)} monomers in {len(arrays)} arrays")
    return mon, arrays.reset_index(drop=True)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest -q tests/test_annotate.py`
Expected: PASS (3 tests).

---

### Task 5: S4 — classify families (`cenrep/classify.py`)

**Files:**
- Create: `clean/cenrep/cenrep/classify.py`
- Test: `clean/cenrep/tests/test_classify.py`

**Interfaces:**
- Consumes: arrays (`ARR_COLS`), `chroms`, `genes` (`dict|None`), `Params`, `util.union_bp/clipped_bp`, `genomes.count_in`
- Produces:
  - `Cluster(start, end, bp)`, `make_clusters(iv, gap) -> list[Cluster]` (sorted bp desc), `end_zone(L, p) -> int`, `interval_gap(a, b) -> int`
  - `CLS_COLS = [family, bp, n_arrays, prevalence, concentration, subtel_frac, gene_ratio, median_main_bp, eligible, cls, strength, defines, reason]`
  - `family_metrics(arrays, chroms, genes, p) -> (DataFrame CLS_COLS, mains: dict[(fam, chrom), Cluster])`
  - `choose_centromeric(fams, mains, chroms, p) -> (fams with cls/strength/defines, per_chrom: dict[chrom, list[fam]])`
  - `assign_cen_associated(fams, arrays, regions: dict[chrom, row dict], p) -> fams`. Region rows use keys `evidence`, `cen_start`, `cen_end`; `evidence == "none"` means no region.

- [ ] **Step 1: Write the failing tests** — `tests/test_classify.py`

```python
import numpy as np
import pandas as pd

from cenrep.annotate import ARR_COLS
from cenrep.classify import assign_cen_associated, choose_centromeric, family_metrics, make_clusters
from cenrep.params import Params

P = Params()
CH = {f"c{i}": 10_000_000 for i in range(1, 5)}


def arr(fam, chrom, start, n, mlen=180):
    return (0, chrom, start, start + n * mlen, "+", fam, n, 0.9, mlen)


def test_make_clusters_bp_and_order():
    cl = make_clusters([(0, 100), (50, 200), (10_000, 10_050), (600_000, 600_300)], 500)
    assert [(c.start, c.end, c.bp) for c in cl] == [(600_000, 600_300, 300), (0, 200, 200), (10_000, 10_050, 50)]


def test_chromosome_specific_and_layered_families():
    rows = [arr("A", "c1", 5_000_000, 1000), arr("A", "c2", 5_000_000, 1000),
            arr("B", "c3", 4_000_000, 800), arr("B", "c4", 4_000_000, 800),
            arr("C", "c1", 5_300_000, 200),                                      # layered next to A on c1
            arr("D", "c1", 100_000, 300), arr("D", "c2", 9_700_000, 300),
            arr("D", "c3", 50_000, 300), arr("D", "c4", 9_800_000, 300)]         # subtelomeric
    rows += [arr("S", f"c{1 + i % 4}", 300_000 + (i // 4) * 900_000, 120) for i in range(40)]   # scattered
    fams, mains = family_metrics(pd.DataFrame(rows, columns=ARR_COLS), CH, None, P)
    fams, per_chrom = choose_centromeric(fams, mains, CH, P)
    f = fams.set_index("family")
    assert f.loc["A", "cls"] == f.loc["B", "cls"] == f.loc["C", "cls"] == "centromeric"
    assert f.loc["D", "cls"] == "other" and "subtel_frac" in f.loc["D", "reason"]
    assert f.loc["S", "cls"] == "other" and "concentration" in f.loc["S", "reason"]
    assert per_chrom["c1"] == ["A", "C"] and per_chrom["c3"] == ["B"]
    assert f.loc["A", "strength"] == "strong" and f.loc["C", "strength"] == "weak"


def test_gene_rich_cluster_is_vetoed():
    genes = {c: np.arange(0, 10_000_000, 5_000, dtype=np.int64) for c in ("c1", "c2")}
    fams, _ = family_metrics(pd.DataFrame([arr("G", "c1", 5_000_000, 500)], columns=ARR_COLS), CH, genes, P)
    assert not fams.eligible.iloc[0] and "gene_ratio" in fams.reason.iloc[0]


def test_empty_arrays():
    fams, mains = family_metrics(pd.DataFrame(columns=ARR_COLS), CH, None, P)
    fams, per_chrom = choose_centromeric(fams, mains, CH, P)
    assert fams.empty and per_chrom == {}


def test_assign_cen_associated():
    rows = [arr("A", "c1", 5_000_000, 1000), arr("X", "c1", 5_050_000, 100), arr("X", "c2", 1_500_000, 20)]
    arrays = pd.DataFrame(rows, columns=ARR_COLS)
    fams, mains = family_metrics(arrays, CH, None, P)
    fams, _ = choose_centromeric(fams, mains, CH, P)
    regions = {"c1": dict(evidence="satellite", cen_start=5_000_000, cen_end=5_180_000),
               "c2": dict(evidence="none")}
    out = assign_cen_associated(fams, arrays, regions, P).set_index("family")
    assert out.loc["X", "cls"] == "cen_associated"                            # 18 of 21.6 kb inside
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest -q tests/test_classify.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'cenrep.classify'`

- [ ] **Step 3: Write `cenrep/classify.py`**

```python
"""S4: per-family metrics and the per-chromosome centromeric-family choice."""
import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .genomes import count_in
from .util import clipped_bp, union_bp

CLS_COLS = ["family", "bp", "n_arrays", "prevalence", "concentration", "subtel_frac", "gene_ratio",
            "median_main_bp", "eligible", "cls", "strength", "defines", "reason"]


@dataclass
class Cluster:
    start: int
    end: int
    bp: int


def make_clusters(iv, gap: int) -> list:
    """Arrays closer than gap form one cluster; bp = union bp of member arrays. Sorted by bp desc, start."""
    out, members, run_end = [], [], -1
    for s, e in sorted(iv):
        if members and s - run_end > gap:
            out.append(Cluster(members[0][0], run_end, union_bp(members)))
            members, run_end = [], -1
        members.append((s, e))
        run_end = max(run_end, e)
    if members:
        out.append(Cluster(members[0][0], run_end, union_bp(members)))
    return sorted(out, key=lambda c: (-c.bp, c.start))


def end_zone(L: int, p) -> int:
    return max(p.subtel_min_bp, int(p.subtel_frac_len * L))


def end_zone_bp(iv, L: int, p) -> int:
    z = end_zone(L, p)
    if 2 * z >= L:
        return union_bp(iv)
    return clipped_bp(iv, 0, z) + clipped_bp(iv, L - z, L)


def interval_gap(a, b) -> int:
    return max(0, max(a.start, b.start) - min(a.end, b.end))


def family_metrics(arrays: pd.DataFrame, chroms: dict, genes, p):
    total = sum(chroms.values())
    n_genes = sum(len(v) for v in genes.values()) if genes else 0
    density = n_genes / total if n_genes else None
    rows, mains = [], {}
    for fam, g in arrays.groupby("family"):
        bp = subtel = gcount = gspan = 0
        shares, main_bps = [], []
        for c, d in g.groupby("chrom"):
            iv = list(zip(d.start, d.end))
            bpc = union_bp(iv)
            bp += bpc
            subtel += end_zone_bp(iv, chroms[c], p)
            m = make_clusters(iv, p.cluster_gap)[0]
            mains[(fam, c)] = m
            if bpc >= p.chrom_min_bp:
                shares.append(m.bp / bpc)
                main_bps.append(m.bp)
                if density:
                    gcount += count_in(genes.get(c), m.start, m.end)
                    gspan += m.end - m.start
        prevalence = len(shares) / len(chroms)
        conc = float(np.median(shares)) if shares else float("nan")
        subtel_frac = subtel / bp if bp else float("nan")
        gene_ratio = (gcount / gspan) / density if (density and gspan) else float("nan")
        fails = []
        if bp < p.chrom_min_bp:
            fails.append(f"bp={bp}<{p.chrom_min_bp}")
        if not conc >= p.min_concentration:
            fails.append(f"concentration={conc:.2f}<{p.min_concentration}")
        if not subtel_frac < p.max_subtel_frac:
            fails.append(f"subtel_frac={subtel_frac:.2f}>={p.max_subtel_frac}")
        if not math.isnan(gene_ratio) and gene_ratio >= p.max_gene_ratio:
            fails.append(f"gene_ratio={gene_ratio:.2f}>={p.max_gene_ratio}")
        rows.append(dict(family=fam, bp=bp, n_arrays=len(g), prevalence=prevalence, concentration=conc,
                         subtel_frac=subtel_frac, gene_ratio=gene_ratio,
                         median_main_bp=float(np.median(main_bps)) if main_bps else 0.0,
                         eligible=not fails, cls="other", strength=".", defines=".",
                         reason="; ".join(fails) if fails else "eligible"))
    return pd.DataFrame(rows, columns=CLS_COLS), mains


def choose_centromeric(fams: pd.DataFrame, mains: dict, chroms: dict, p):
    elig = set(fams.loc[fams.eligible.astype(bool), "family"])
    per_chrom = {}
    for c in chroms:
        cands = sorted(((mains[(f, c)].bp, f) for f in elig
                        if (f, c) in mains and mains[(f, c)].bp >= p.chrom_min_bp), reverse=True)
        if not cands:
            continue
        best = cands[0][1]
        bc = mains[(best, c)]
        per_chrom[c] = [best] + [f for _, f in cands[1:] if interval_gap(mains[(f, c)], bc) <= p.cluster_gap]
    defines = {}
    for c, fl in per_chrom.items():
        for f in fl:
            defines.setdefault(f, []).append(c)
    fams = fams.copy()
    for i, r in fams.iterrows():
        if r.family in defines:
            strong = r.prevalence >= p.strong_prevalence and r.median_main_bp >= p.strong_cluster_bp
            fams.loc[i, ["cls", "strength", "defines", "reason"]] = [
                "centromeric", "strong" if strong else "weak", ",".join(defines[r.family]),
                f"eligible; centromeric family on {len(defines[r.family])} chromosome(s)"]
        elif r.eligible:
            fams.loc[i, "reason"] = "eligible; never the main or co-located cluster on a chromosome"
    return fams, per_chrom


def assign_cen_associated(fams: pd.DataFrame, arrays: pd.DataFrame, regions: dict, p) -> pd.DataFrame:
    reg = {c: (r["cen_start"], r["cen_end"]) for c, r in regions.items() if r["evidence"] != "none"}
    fams = fams.copy()
    for i, r in fams.iterrows():
        if r.cls == "centromeric" or not r.bp:
            continue
        d = arrays[arrays.family == r.family]
        inside = sum(clipped_bp(list(zip(x.start, x.end)), *reg[c]) for c, x in d.groupby("chrom") if c in reg)
        if inside / r.bp >= p.cen_assoc_frac:
            fams.loc[i, ["cls", "reason"]] = ["cen_associated",
                                              f"{r.reason}; {inside / r.bp:.0%} of bp inside centromeres"]
    return fams
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest -q tests/test_classify.py`
Expected: PASS (5 tests).

---

### Task 6: S5b — long-unit tandem arrays (`cenrep/longunit.py`)

**Files:**
- Create: `clean/cenrep/cenrep/longunit.py`
- Test: `clean/cenrep/tests/test_longunit.py`

**Interfaces:**
- Consumes: `Params.lu_*`, `util.run/merge_intervals`
- Produces:
  - `LU_COLS = [chrom, start, end, unit, copies, identity]`, `windows(start, end, size, overlap)`
  - `mask_seq(seq, offset, starts, ends) -> str`, `selfaln_hits(seq, offset, tmpdir, p) -> list[(start, end, d, ident)]`, `merge_hits(hits, chrom) -> list[dict]`
  - `build_masks(frames) -> dict[chrom, (starts ndarray, ends ndarray)]`
  - `detect(fasta, targets: list[(chrom, start, end)], masks, tmpdir, threads, p) -> DataFrame LU_COLS`

- [ ] **Step 1: Write the failing tests** — `tests/test_longunit.py`

```python
import numpy as np
import pandas as pd
import pysam
import pytest

from cenrep.longunit import build_masks, detect, mask_seq, merge_hits, windows
from conftest import have
from synth import mutate, rseq


def test_windows():
    assert list(windows(0, 250, 100, 20)) == [(0, 100), (80, 180), (160, 250)]


def test_mask_seq_clips_to_window():
    starts, ends = np.array([5, 20]), np.array([8, 30])
    assert mask_seq("A" * 25, 3, starts, ends) == "AA" + "NNN" + "A" * 12 + "N" * 8


def test_build_masks_merges():
    df = pd.DataFrame({"chrom": ["c", "c", "d"], "start": [0, 5, 1], "end": [10, 20, 2]})
    m = build_masks([df])
    assert m["c"][0].tolist() == [0] and m["c"][1].tolist() == [20] and m["d"][1].tolist() == [2]


def test_merge_hits_multiples_and_distinct_arrays():
    hits = [(1000, 16000, 3000, 0.97), (1000, 16000, 6000, 0.96), (1100, 15800, 9000, 0.95),
            (15500, 30000, 4100, 0.90)]
    out = merge_hits(hits, "c")
    assert sorted((r["unit"], r["start"], r["end"]) for r in out) == [(3000, 1000, 16000), (4100, 15500, 30000)]


@pytest.mark.skipif(not have("minimap2"), reason="minimap2 missing")
def test_detect_finds_3kb_unit(tmp_path, test_params):
    rng = np.random.default_rng(9)
    u = rseq(rng, 3000)
    seq = rseq(rng, 30000) + "".join(mutate(rng, u, 0.03) for _ in range(6)) + rseq(rng, 30000)
    fa = tmp_path / "x.fa"
    fa.write_text(">c\n" + seq + "\n")
    pysam.faidx(str(fa))
    lu = detect(fa, [("c", 0, len(seq))], {}, tmp_path, 2, test_params)
    assert len(lu) == 1 and abs(lu.unit.iloc[0] - 3000) <= 150
    assert abs(lu.start.iloc[0] - 30000) <= 300 and abs(lu.end.iloc[0] - 48000) <= 300
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest -q tests/test_longunit.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'cenrep.longunit'`

- [ ] **Step 3: Write `cenrep/longunit.py`**

```python
"""S5b: tandem arrays with units > 2 kb (beyond TRF's MaxPeriod), from windowed minimap2 self-alignment.
A tandem array of unit d self-aligns with a shift of d; an alignment at least d long means >= 2 adjacent copies."""
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import pysam

from .util import log, merge_intervals, run

LU_COLS = ["chrom", "start", "end", "unit", "copies", "identity"]


def windows(start: int, end: int, size: int, overlap: int):
    s = start
    while s < end:
        e = min(end, s + size)
        yield s, e
        if e == end:
            break
        s = e - overlap


def build_masks(frames) -> dict:
    """Merged intervals to N-mask (satellite arrays, TRF rows): {chrom: (starts, ends)}."""
    iv = {}
    for df in frames:
        for c, s, e in zip(df.chrom, df.start, df.end):
            iv.setdefault(c, []).append((int(s), int(e)))
    out = {}
    for c, x in iv.items():
        m = np.array(merge_intervals(x), dtype=np.int64).reshape(-1, 2)
        out[c] = (m[:, 0], m[:, 1])
    return out


def mask_seq(seq: str, offset: int, starts, ends) -> str:
    b = bytearray(seq.encode())
    n = len(b)
    i0 = int(np.searchsorted(ends, offset, side="right"))
    i1 = int(np.searchsorted(starts, offset + n, side="left"))
    for s, e in zip(starts[i0:i1], ends[i0:i1]):
        a, z = max(int(s) - offset, 0), min(int(e) - offset, n)
        if z > a:
            b[a:z] = b"N" * (z - a)
    return b.decode()


def selfaln_hits(seq: str, offset: int, tmpdir, p) -> list:
    fd, path = tempfile.mkstemp(suffix=".fa", dir=tmpdir)
    with os.fdopen(fd, "w") as f:
        f.write(f">w\n{seq}\n")
    try:
        out = run(["minimap2", "-x", "asm20", "-c", "-DP", "--dual=no", "-t", "1", path, path]).stdout
    finally:
        os.unlink(path)
    hits = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12 or f[4] != "+":
            continue
        qs, qe, ts, te, nm, al = int(f[2]), int(f[3]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
        d = ts - qs
        if d <= 0 or not p.lu_min_unit <= d <= p.lu_max_unit or qe - qs < d or nm / al < p.lu_min_ident:
            continue
        hits.append((offset + qs, offset + te, d, nm / al))
    return hits


def merge_hits(hits, chrom: str) -> list:
    """Overlapping hits form a group; within it each hit joins the smallest base unit it is a multiple of
    (5%), and each base unit's spans are unioned into arrays."""
    out, groups, cur, cur_end = [], [], [], -1
    for h in sorted(hits):
        if cur and h[0] >= cur_end:
            groups.append(cur)
            cur, cur_end = [], -1
        cur.append(h)
        cur_end = max(cur_end, h[1])
    if cur:
        groups.append(cur)
    for g in groups:
        bases = {}
        for h in sorted(g, key=lambda x: x[2]):
            for b in bases:
                k = round(h[2] / b)
                if k >= 1 and abs(h[2] / b - k) <= 0.05 * k:
                    bases[b].append(h)
                    break
            else:
                bases[h[2]] = [h]
        for b, hs in bases.items():
            unit_hits = [h for h in hs if round(h[2] / b) == 1]
            unit = int(np.median([h[2] for h in unit_hits]))
            ident = float(np.mean([h[3] for h in unit_hits]))
            for s, e in merge_intervals((h[0], h[1]) for h in hs):
                out.append(dict(chrom=chrom, start=s, end=e, unit=unit,
                                copies=round((e - s) / unit, 2), identity=round(ident, 4)))
    return out


def detect(fasta, targets, masks: dict, tmpdir, threads: int, p) -> pd.DataFrame:
    jobs = [(c, ws, we) for c, s, e in targets for ws, we in windows(s, e, p.lu_window, p.lu_window_overlap)]
    empty = (np.zeros(0, np.int64), np.zeros(0, np.int64))

    def work(job):
        c, ws, we = job
        with pysam.FastaFile(str(fasta)) as fa:
            seq = mask_seq(fa.fetch(c, ws, we).upper(), ws, *masks.get(c, empty))
        if len(seq) - seq.count("N") < 2 * p.lu_min_unit:
            return c, []
        return c, selfaln_hits(seq, ws, tmpdir, p)

    per = {}
    with ThreadPoolExecutor(max(1, threads)) as ex:
        for c, h in ex.map(work, jobs):
            per.setdefault(c, []).extend(h)
    rows = [r for c, h in per.items() for r in merge_hits(h, c)]
    log.debug(f"longunit: {len(jobs)} windows, {len(rows)} arrays")
    return pd.DataFrame(rows, columns=LU_COLS)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest -q tests/test_longunit.py`
Expected: PASS (5 tests).

---

### Task 7: S5 — centromere region calls (`cenrep/centromere.py`)

**Files:**
- Create: `clean/cenrep/cenrep/centromere.py`
- Test: `clean/cenrep/tests/test_centromere.py`

**Interfaces:**
- Consumes: `classify.make_clusters/end_zone/end_zone_bp`, `longunit.detect/build_masks/LU_COLS`, `genomes.count_in`
- Produces:
  - `REG_COLS = [chrom, chrom_len, cen_start, cen_end, evidence, confidence, flags, main_family, main_family_bp, sat_bp, second_cluster_bp, cen_families, n_gaps, gap_bp, n_genes, gene_ratio]`
  - `satellite_call(arrays_c, primary, L, p) -> dict|None`, `longunit_call(lu_c, L, genes_c, gene_density, p) -> dict|None`
  - `finalize(chrom, L, call, fasta_handle, genes_c, chrom_gene_density, p) -> row dict`
    - Every chromosome gets a row; with no call, `evidence="none"`, `confidence="none"` and the other fields `"."`.
    - `flags` is a comma-joined list from `low_satellite`, `competing_cluster`, `near_end`, `gene_rich`, `gaps:<n>`, or `"."`.
  - `call_regions(genome, arrays, primary, genes, trf_rows, tmpdir, threads, p) -> (rows: dict[chrom, row], lu: DataFrame)`
  - `write_regions(rows, path)`, `load_regions(path) -> dict[chrom, row]`. The loader converts `cen_start`, `cen_end` and `chrom_len` to int where present.

- [ ] **Step 1: Write the failing tests** — `tests/test_centromere.py`

```python
import pandas as pd
import pysam

from cenrep.annotate import ARR_COLS
from cenrep.centromere import finalize, load_regions, longunit_call, satellite_call, write_regions
from cenrep.longunit import LU_COLS
from cenrep.params import Params


def test_satellite_call_main_cluster_and_competing():
    a = pd.DataFrame([(0, "c", 5_000_000, 5_200_000, "+", "A", 1100, 0.9, 180),
                      (1, "c", 5_300_000, 5_350_000, "+", "A", 270, 0.9, 180),
                      (2, "c", 8_000_000, 8_120_000, "+", "A", 660, 0.9, 180)], columns=ARR_COLS)
    call = satellite_call(a, {"A"}, 10_000_000, Params())
    assert (call["cen_start"], call["cen_end"]) == (5_000_000, 5_350_000)
    assert call["sat_bp"] == 250_000 and call["second_cluster_bp"] == 120_000 and call["evidence"] == "satellite"


def test_satellite_call_weak_and_absent():
    p = Params()
    small = pd.DataFrame([(0, "c", 100_000, 108_000, "+", "A", 44, 0.9, 180)], columns=ARR_COLS)
    assert satellite_call(small, {"A"}, 10_000_000, p)["evidence"] == "weak_satellite"
    assert satellite_call(small, {"B"}, 10_000_000, p) is None


def test_finalize_flags_gaps_and_near_end(tmp_path):
    seq = "A" * 1000 + "N" * 20 + "A" * 1000 + "N" * 5 + "A" * 2000
    fa = tmp_path / "x.fa"
    fa.write_text(">c\n" + seq + "\n")
    pysam.faidx(str(fa))
    p = Params(subtel_min_bp=500)
    call = dict(cen_start=100, cen_end=3000, evidence="satellite", main_family="A", main_family_bp=2000,
                sat_bp=2000, second_cluster_bp=0, cen_families="A")
    with pysam.FastaFile(str(fa)) as h:
        row = finalize("c", len(seq), call, h, None, None, p)
    assert row["n_gaps"] == 1 and row["gap_bp"] == 20
    assert row["flags"] == "low_satellite,near_end,gaps:1" and row["confidence"] == "low"


def test_none_row_when_no_call():
    row = finalize("c", 1000, None, None, None, None, Params())
    assert row["evidence"] == "none" and row["confidence"] == "none" and row["cen_start"] == "."


def test_longunit_call_requires_concentration():
    p = Params()
    lu = pd.DataFrame([("c", 4_000_000, 4_600_000, 2700, 222.2, 0.95)], columns=LU_COLS)
    call = longunit_call(lu, 10_000_000, None, None, p)
    assert call["evidence"] == "long_unit_tandem" and call["cen_start"] == 4_000_000
    scattered = pd.DataFrame([("c", s, s + 30_000, 2700, 11.1, 0.9) for s in range(1_000_000, 9_000_000, 1_000_000)],
                             columns=LU_COLS)
    assert longunit_call(scattered, 10_000_000, None, None, p) is None


def test_regions_roundtrip(tmp_path):
    rows = {"c1": finalize("c1", 1000, None, None, None, None, Params())}
    rows["c2"] = dict(rows["c1"], chrom="c2", cen_start=10, cen_end=500, evidence="satellite", confidence="high")
    write_regions(rows, tmp_path / "r.tsv")
    back = load_regions(tmp_path / "r.tsv")
    assert back["c1"]["evidence"] == "none" and back["c2"]["cen_end"] == 500 and back["c2"]["chrom_len"] == 1000
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest -q tests/test_centromere.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'cenrep.centromere'`

- [ ] **Step 3: Write `cenrep/centromere.py`**

```python
"""S5: per-chromosome centromere region call: satellite clusters first, then tandem-only fallbacks."""
import re

import numpy as np
import pandas as pd
import pysam

from . import longunit
from .classify import end_zone, end_zone_bp, make_clusters
from .genomes import count_in
from .util import log, union_bp

REG_COLS = ["chrom", "chrom_len", "cen_start", "cen_end", "evidence", "confidence", "flags", "main_family",
            "main_family_bp", "sat_bp", "second_cluster_bp", "cen_families", "n_gaps", "gap_bp", "n_genes",
            "gene_ratio"]
EMPTY = np.zeros(0, np.int64)


def satellite_call(arrays_c: pd.DataFrame, primary: set, L: int, p):
    a = arrays_c[arrays_c.family.isin(primary)]
    if a.empty:
        return None
    cl = make_clusters(list(zip(a.start, a.end)), p.cluster_gap)
    m = cl[0]
    if m.bp >= p.chrom_min_bp:
        ev = "satellite"
    elif m.bp >= p.fallback_min_bp:
        ev = "weak_satellite"
    else:
        return None
    inside = a[(a.start >= m.start) & (a.end <= m.end)]
    fam_bp = {f: union_bp(zip(g.start, g.end)) for f, g in inside.groupby("family")}
    main = max(fam_bp, key=lambda f: (fam_bp[f], f))
    return dict(cen_start=int(m.start), cen_end=int(m.end), evidence=ev, main_family=main,
                main_family_bp=fam_bp[main], sat_bp=m.bp, second_cluster_bp=cl[1].bp if len(cl) > 1 else 0,
                cen_families=",".join(sorted(fam_bp)))


def longunit_call(lu_c: pd.DataFrame, L: int, genes_c, gene_density, p):
    """Fallback 2: long-unit arrays that pass the S4 eligibility tests define a low-confidence region."""
    if lu_c.empty:
        return None
    iv = list(zip(lu_c.start, lu_c.end))
    bpc = union_bp(iv)
    cl = make_clusters(iv, p.cluster_gap)
    m = cl[0]
    if m.bp < p.chrom_min_bp or m.bp / bpc < p.min_concentration or end_zone_bp(iv, L, p) / bpc >= p.max_subtel_frac:
        return None
    if genes_c is not None and gene_density:
        if (count_in(genes_c, m.start, m.end) / (m.end - m.start)) / gene_density >= p.max_gene_ratio:
            return None
    return dict(cen_start=int(m.start), cen_end=int(m.end), evidence="long_unit_tandem", main_family="long_unit",
                main_family_bp=m.bp, sat_bp=m.bp, second_cluster_bp=cl[1].bp if len(cl) > 1 else 0,
                cen_families="long_unit")


def finalize(chrom: str, L: int, call, fa, genes_c, chrom_density, p) -> dict:
    row = {k: "." for k in REG_COLS}
    row.update(chrom=chrom, chrom_len=L, evidence="none", confidence="none")
    if call is None:
        return row
    s, e = call["cen_start"], call["cen_end"]
    runs = [m.end() - m.start() for m in re.finditer("[Nn]{%d,}" % p.min_gap_n, fa.fetch(chrom, s, e))]
    flags = []
    if call["sat_bp"] < p.low_sat_bp:
        flags.append("low_satellite")
    if call["second_cluster_bp"] > 0 and call["second_cluster_bp"] >= p.competing_ratio * call["sat_bp"]:
        flags.append("competing_cluster")
    z = end_zone(L, p)
    if s < z or e > L - z:
        flags.append("near_end")
    core = len(flags)
    n_genes, gene_ratio = ".", "."
    if genes_c is not None:
        n_genes = count_in(genes_c, s, e)
        if chrom_density:
            gene_ratio = round((n_genes / (e - s)) / chrom_density, 3)
            if gene_ratio > 1:
                flags.append("gene_rich")
    if runs:
        flags.append(f"gaps:{len(runs)}")
    conf = ["high", "medium", "low", "low"][core] if call["evidence"] == "satellite" else "low"
    row.update(call)
    row.update(confidence=conf, flags=",".join(flags) or ".", n_gaps=len(runs), gap_bp=sum(runs),
               n_genes=n_genes, gene_ratio=gene_ratio)
    return row


def call_regions(genome, arrays, primary: set, genes, trf_rows, tmpdir, threads: int, p):
    tmpdir.mkdir(parents=True, exist_ok=True)
    total = sum(genome.chroms.values())
    gd = (sum(len(v) for v in genes.values()) / total) if genes else None
    sat = {c: satellite_call(arrays[arrays.chrom == c], primary, L, p) for c, L in genome.chroms.items()}
    targets = [(c, max(0, sat[c]["cen_start"] - p.lu_flank), min(L, sat[c]["cen_end"] + p.lu_flank))
               if sat[c] else (c, 0, L) for c, L in genome.chroms.items()]
    lu = longunit.detect(genome.fasta, targets, longunit.build_masks([arrays, trf_rows]), tmpdir, threads, p)
    rows = {}
    with pysam.FastaFile(str(genome.fasta)) as fa:
        for c, L in genome.chroms.items():
            gc = None if genes is None else genes.get(c, EMPTY)
            call = sat[c] or longunit_call(lu[lu.chrom == c], L, gc, gd, p)
            rows[c] = finalize(c, L, call, fa, gc, (len(gc) / L) if gc is not None else None, p)
    log.debug(f"centromeres: {genome.name}: {len(lu)} long-unit arrays")
    return rows, lu


def write_regions(rows: dict, path) -> None:
    pd.DataFrame(list(rows.values()), columns=REG_COLS).to_csv(path, sep="\t", index=False)


def load_regions(path) -> dict:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    out = {}
    for r in df.to_dict("records"):
        r["chrom_len"] = int(r["chrom_len"])
        if r["evidence"] != "none":
            r["cen_start"], r["cen_end"] = int(r["cen_start"]), int(r["cen_end"])
        out[r["chrom"]] = r
    return out
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest -q tests/test_centromere.py`
Expected: PASS (6 tests).

---

### Task 8: S6 — HORs (`cenrep/hor.py`)

**Files:**
- Create: `clean/cenrep/cenrep/hor.py`
- Test: `clean/cenrep/tests/test_hor.py`

**Interfaces:**
- Consumes:
  - monomers table (`mon_idx` == row position, plus `MON_COLS` and `array_id`), arrays (`ARR_COLS`)
  - regions (`dict[chrom, row]`), `primary: set[fam_id]`, `Genome`, `Params.hor_*`
- Produces:
  - `suffix_array(t) -> ndarray`, `lcp_array(t, sa) -> ndarray`
  - `maximal_repeats(t, sa, lcp, max_len, min_len=2, keep=None) -> list[(len, occ ndarray)]`
  - `smallest_period(u) -> int`, `is_hor_unit(u) -> bool`
  - `assign_hors(t, reps) -> list[(unit tuple, sorted occurrence list)]`, `tandem_runs(acc, l) -> list[list]`
  - `variant_classes(seqs: dict[name, seq], workdir, threads, ident) -> dict[name, cluster int]`
  - `run_hor(genome, mon, arrays, regions, primary, workdir, threads, p) -> (mon_hor, occ, units)`. The three DataFrames have these columns:
    - `MONHOR_COLS = [mon_idx, variant_class, hor_unit, hor_occ]`
    - `OCC_COLS = [occ_id, unit, chrom, start, end, strand, array_id, n_monomers, tandem_run, tandem_copies]`
    - `UNIT_COLS = [unit, family, chrom, n_monomers, occurrences, max_tandem_copies, class_string]`

- [ ] **Step 1: Write the failing tests** — `tests/test_hor.py`

```python
import time

import numpy as np
import pandas as pd
import pytest

from cenrep.annotate import ARR_COLS, MON_COLS
from cenrep.genomes import Genome
from cenrep.hor import (assign_hors, is_hor_unit, lcp_array, maximal_repeats, run_hor, smallest_period,
                        suffix_array, tandem_runs, variant_classes)
from cenrep.params import Params
from conftest import have
from synth import mutate, rseq


def brute_sa(t):
    return sorted(range(len(t)), key=lambda i: list(t[i:]))


def brute_maximal(t, min_len=2):
    n, out = len(t), set()
    for l in range(min_len, n + 1):
        occ = {}
        for i in range(n - l + 1):
            occ.setdefault(tuple(t[i:i + l]), []).append(i)
        for ps in occ.values():
            if len(ps) < 2:
                continue
            left = {t[q - 1] if q > 0 else ("start", q) for q in ps}
            right = {t[q + l] if q + l < n else ("end", q) for q in ps}
            if len(left) > 1 and len(right) > 1:
                out.add((l, tuple(ps)))
    return out


def test_suffix_array_and_lcp_match_brute_force():
    rng = np.random.default_rng(1)
    for _ in range(60):
        t = rng.integers(0, 3, size=int(rng.integers(1, 40)))
        if rng.random() < 0.5:
            t[rng.integers(0, len(t))] = -1
        sa = suffix_array(t)
        assert sa.tolist() == brute_sa(t.tolist())
        lcp = lcp_array(t, sa)
        for i in range(1, len(t)):
            a, b = t[sa[i - 1]:].tolist(), t[sa[i]:].tolist()
            k = 0
            while k < min(len(a), len(b)) and a[k] == b[k]:
                k += 1
            assert lcp[i] == k


def test_maximal_repeats_match_brute_force():
    rng = np.random.default_rng(2)
    for _ in range(60):
        t = rng.integers(0, 3, size=int(rng.integers(2, 30)))
        sa = suffix_array(t)
        got = {(l, tuple(sorted(int(x) for x in occ)))
               for l, occ in maximal_repeats(t, sa, lcp_array(t, sa), max_len=10 ** 6)}
        assert got == brute_maximal(t.tolist())


def test_smallest_period_and_unit_rules():
    assert smallest_period((1, 2, 3, 1, 2, 3)) == 3 and smallest_period((1, 2, 1)) == 2
    assert is_hor_unit((1, 2, 3, 4)) and is_hor_unit((1, 2, 1)) and is_hor_unit((1, 2, 3, 4, 1, 2, 3))
    assert not is_hor_unit((1, 2, 1, 2)) and not is_hor_unit((5, 5, 5)) and not is_hor_unit((7,))


def _hors(t, max_len=200):
    t = np.array(t)
    sa = suffix_array(t)
    return assign_hors(t, maximal_repeats(t, sa, lcp_array(t, sa), max_len, keep=is_hor_unit))


def test_tandem_hor_prefers_full_unit_over_partial_repeat():
    units = _hors([-1] + [0, 1, 2, 3] * 14 + [0, 1, 2] + [-2])
    assert len(units) == 1
    u, acc = units[0]
    assert u == (0, 1, 2, 3) and len(acc) == 14 and tandem_runs(acc, 4) == [acc]


def test_non_tandem_duplication():
    assert _hors([-1, 9, 0, 1, 10, 11, 0, 1, 12, -2]) == [((0, 1), [2, 6])]


def test_homogeneous_run_is_not_a_hor_and_is_fast():
    t0 = time.time()
    assert _hors([-1] + [0] * 20000 + [-2]) == []
    assert time.time() - t0 < 10


@pytest.mark.skipif(not have("cd-hit-est"), reason="cd-hit-est missing")
def test_variant_classes_groups_near_identical(tmp_path):
    rng = np.random.default_rng(4)
    base = rseq(rng, 150)
    v = [mutate(rng, base, 0.08) for _ in range(3)]
    cl = variant_classes({f"m{i}": mutate(rng, v[i % 3], 0.005) for i in range(12)}, tmp_path, 2, 0.97)
    groups = {}
    for name, k in cl.items():
        groups.setdefault(k, set()).add(int(name[1:]) % 3)
    assert len(groups) == 3 and all(len(g) == 1 for g in groups.values())


@pytest.mark.skipif(not have("cd-hit-est"), reason="cd-hit-est missing")
def test_run_hor_tiny_region(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">c\n" + "ACGT" * 250 + "\n")
    import pysam
    pysam.faidx(str(fa))
    g = Genome("g", fa, {"c": 1000})
    mon = pd.DataFrame([(0, "c", 0, 150, "+", "f", 90.0, 1e-20, 1, 150, 150, 0.9, 0)],
                       columns=["mon_idx"] + MON_COLS + ["array_id"])
    arrays = pd.DataFrame([(0, "c", 0, 150, "+", "f", 1, 0.9, 150)], columns=ARR_COLS)
    regions = {"c": dict(evidence="satellite", cen_start=0, cen_end=150)}
    mh, occ, units = run_hor(g, mon, arrays, regions, {"f"}, tmp_path / "w", 1, Params())
    assert occ.empty and units.empty and len(mh) == 1
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest -q tests/test_hor.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'cenrep.hor'`

- [ ] **Step 3: Write `cenrep/hor.py`**

```python
"""S6: monomer variant classes (cd-hit-est) and higher-order repeats (HORs) via maximal repeats.

Each centromere becomes an integer string (one symbol per monomer class, consensus orientation, unique
negative separators between arrays and for unclassed monomers). A HOR unit is a maximal repeat of >= 2
monomers with >= 2 distinct classes that does not contain two full copies of a shorter period."""
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from .util import log, revcomp, run, write_fasta

MONHOR_COLS = ["mon_idx", "variant_class", "hor_unit", "hor_occ"]
OCC_COLS = ["occ_id", "unit", "chrom", "start", "end", "strand", "array_id", "n_monomers", "tandem_run",
            "tandem_copies"]
UNIT_COLS = ["unit", "family", "chrom", "n_monomers", "occurrences", "max_tandem_copies", "class_string"]


def suffix_array(t) -> np.ndarray:
    """Prefix doubling, O(n log^2 n) with numpy sorts. Shorter suffix sorts first on ties."""
    t = np.asarray(t)
    n = len(t)
    if n == 0:
        return np.zeros(0, np.int64)
    rank = np.unique(t, return_inverse=True)[1].astype(np.int64).ravel()
    k = 1
    while True:
        r2 = np.full(n, -1, np.int64)
        if k < n:
            r2[:n - k] = rank[k:]
        sa = np.lexsort((r2, rank))
        a, b = rank[sa], r2[sa]
        new = np.zeros(n, np.int64)
        new[1:] = np.cumsum((a[1:] != a[:-1]) | (b[1:] != b[:-1]))
        rank = np.empty(n, np.int64)
        rank[sa] = new
        if new[-1] == n - 1:
            return sa
        k *= 2


def lcp_array(t, sa) -> np.ndarray:
    """Kasai: lcp[i] = LCP(suffix sa[i-1], suffix sa[i]); lcp[0] = 0."""
    tl, sal = [int(x) for x in t], [int(x) for x in sa]
    n = len(tl)
    rank = [0] * n
    for i, s in enumerate(sal):
        rank[s] = i
    lcp, h = [0] * n, 0
    for i in range(n):
        r = rank[i]
        if r == 0:
            h = 0
            continue
        j = sal[r - 1]
        while i + h < n and j + h < n and tl[i + h] == tl[j + h]:
            h += 1
        lcp[r] = h
        if h:
            h -= 1
    return np.array(lcp, np.int64)


def maximal_repeats(t, sa, lcp, max_len: int, min_len: int = 2, keep=None) -> list:
    """Bottom-up lcp-interval traversal (Abouelhoda et al. 2004). A repeat is reported if its interval
    is left-diverse. lcp is capped at max_len, so the work is O(n * max_len)."""
    t = np.asarray(t)
    n = len(t)
    if n < 2:
        return []
    lcp = np.minimum(lcp, max_len)
    prev = np.where(sa > 0, t[np.maximum(sa - 1, 0)], -(10 ** 15) - np.arange(n))
    out, stack = [], [[0, 0]]
    for i in range(1, n + 1):
        cur = int(lcp[i]) if i < n else 0
        lb = i - 1
        while cur < stack[-1][0]:
            l, lb = stack.pop()
            if l >= min_len:
                if keep is None or keep(tuple(int(x) for x in t[sa[lb]:sa[lb] + l])):
                    pv = prev[lb:i]
                    if pv.min() != pv.max():
                        out.append((l, sa[lb:i]))
        if cur > stack[-1][0]:
            stack.append([cur, lb])
    return out


def smallest_period(u) -> int:
    n = len(u)
    if n == 0:
        return 0
    f, k = [0] * n, 0
    for i in range(1, n):
        while k and u[i] != u[k]:
            k = f[k - 1]
        if u[i] == u[k]:
            k += 1
        f[i] = k
    return n - f[-1]


def is_hor_unit(u) -> bool:
    return len(u) >= 2 and min(u) >= 0 and len(set(u)) >= 2 and 2 * smallest_period(u) > len(u)


def _greedy_count(occ, l) -> int:
    n, last = 0, -(10 ** 18)
    for q in occ:
        if q >= last + l:
            n, last = n + 1, q
    return n


def assign_hors(t, reps) -> list:
    """Rank units by achievable non-overlapping coverage, then length; accept occurrences whose monomers
    are unassigned; keep units with >= 2 accepted occurrences."""
    t = np.asarray(t)
    cands = []
    for l, occ in reps:
        occ = sorted(int(x) for x in occ)
        u = tuple(int(x) for x in t[occ[0]:occ[0] + l])
        if is_hor_unit(u):
            c = _greedy_count(occ, l)
            if c >= 2:
                cands.append((c * l, l, occ, u))
    cands.sort(key=lambda x: (-x[0], -x[1], x[2][0]))
    used = np.zeros(len(t), bool)
    units = []
    for _, l, occ, u in cands:
        acc = []
        for q in occ:
            if not used[q:q + l].any():
                used[q:q + l] = True
                acc.append(q)
        if len(acc) >= 2:
            units.append((u, acc))
        else:
            for q in acc:
                used[q:q + l] = False
    return units


def tandem_runs(acc, l) -> list:
    runs = [[acc[0]]]
    for q in acc[1:]:
        if q == runs[-1][-1] + l:
            runs[-1].append(q)
        else:
            runs.append([q])
    return runs


def variant_classes(seqs: dict, workdir, threads: int, ident: float) -> dict:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    fa, out = workdir / "mon.fa", workdir / "mon.cdhit"
    write_fasta(fa, seqs.items())
    run(["cd-hit-est", "-i", fa, "-o", out, "-c", ident, "-n", 10 if ident >= 0.95 else 8, "-G", 1,
         "-aS", 0.9, "-aL", 0.9, "-g", 1, "-r", 0, "-d", 0, "-M", 0, "-T", threads])
    cl, res = -1, {}
    with open(str(out) + ".clstr") as f:
        for line in f:
            if line.startswith(">Cluster"):
                cl = int(line.split()[1])
            else:
                res[line.split(">", 1)[1].split("...", 1)[0]] = cl
    return res


def run_hor(genome, mon: pd.DataFrame, arrays: pd.DataFrame, regions: dict, primary: set, workdir,
            threads: int, p):
    workdir = Path(workdir)
    assert (mon.mon_idx.to_numpy() == np.arange(len(mon))).all(), "mon_idx must equal row position"
    cov = (mon.hmm_to - mon.hmm_from + 1) / mon.model_len
    sel = mon[mon.family.isin(primary) & (cov >= p.hor_min_cov)]
    seqs = {}
    with pysam.FastaFile(str(genome.fasta)) as fa:
        for c, g in sel.groupby("chrom", sort=False):
            cs = fa.fetch(c).upper()
            for i, s, e, st in zip(g.mon_idx, g.start, g.end, g.strand):
                x = cs[s:e]
                seqs[int(i)] = x if st == "+" else revcomp(x)
    klass = {}
    for f in sorted(primary):
        ids = sel.mon_idx[sel.family == f].astype(int).tolist()
        if ids:
            for name, k in variant_classes({f"m{i}": seqs[i] for i in ids}, workdir / f"cdhit_{f}",
                                           threads, p.hor_ident).items():
                klass[int(name[1:])] = (f, k)
    code = {key: n for n, key in enumerate(sorted(set(klass.values())))}
    inv = {v: k for k, v in code.items()}
    m_start, m_end = mon.start.to_numpy(), mon.end.to_numpy()
    by_array = {a: g.mon_idx.astype(int).tolist() for a, g in mon.groupby("array_id", sort=False)}
    hor_unit, hor_occ, occ_rows, unit_rows = {}, {}, [], []
    for c, r in regions.items():
        if r["evidence"] == "none":
            continue
        a = arrays[(arrays.chrom == c) & arrays.family.isin(primary)
                   & (arrays.start >= r["cen_start"]) & (arrays.end <= r["cen_end"])].sort_values("start")
        t, pos2mon, pos2arr, sep = [], [], [], -1
        for aid, st in zip(a.array_id, a.strand):
            idxs = by_array.get(aid, [])
            for i in (idxs[::-1] if st == "-" else idxs):
                if i in klass:
                    t.append(code[klass[i]])
                else:
                    t.append(sep)
                    sep -= 1
                pos2mon.append(i)
                pos2arr.append((aid, st))
            t.append(sep)
            sep -= 1
            pos2mon.append(-1)
            pos2arr.append(None)
        if sum(x >= 0 for x in t) < 4:
            continue
        t = np.array(t, np.int64)
        sa = suffix_array(t)
        units = assign_hors(t, maximal_repeats(t, sa, lcp_array(t, sa), p.hor_max_len, keep=is_hor_unit))
        units.sort(key=lambda x: (-len(x[0]) * len(x[1]), x[1][0]))
        n_occ = 0
        for k, (u, acc) in enumerate(units, 1):
            uid, l = f"HOR_{c}_{k}", len(u)
            runs = tandem_runs(acc, l)
            for ri, rr in enumerate(runs, 1):
                for q in rr:
                    n_occ += 1
                    oid = f"hor_{c}_{n_occ}"
                    mids = pos2mon[q:q + l]
                    aid, st = pos2arr[q]
                    occ_rows.append((oid, uid, c, int(m_start[mids].min()), int(m_end[mids].max()), st, aid, l,
                                     f"{uid}.r{ri}", len(rr)))
                    for mi in mids:
                        hor_unit[mi], hor_occ[mi] = uid, oid
            unit_rows.append((uid, inv[u[0]][0], c, l, len(acc), max(len(rr) for rr in runs),
                              "-".join(str(inv[x][1]) for x in u)))
        log.debug(f"hor: {genome.name} {c}: {len(t)} symbols, {len(units)} HOR units")
    mh = pd.DataFrame([(i, k[1], hor_unit.get(i, "."), hor_occ.get(i, ".")) for i, k in sorted(klass.items())],
                      columns=MONHOR_COLS)
    return mh, pd.DataFrame(occ_rows, columns=OCC_COLS), pd.DataFrame(unit_rows, columns=UNIT_COLS)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest -q tests/test_hor.py`
Expected: PASS (8 tests).

---

### Task 9: S7 — GFF3 and tables (`cenrep/output.py`)

**Files:**
- Create: `clean/cenrep/cenrep/output.py`
- Create: `clean/cenrep/tests/gffcheck.py` (GFF3 validator used by tests)
- Test: `clean/cenrep/tests/test_output.py`

**Interfaces:**
- Consumes: all stage tables; a genome data dict `d` with these keys:
  - `genome`, `arrays`, `monomers`, `fams` (CLS_COLS), `regions`, `longunit`
  - `mon_hor`, `hor_occ`, `hor_units`, `trf`, `genes`
- Produces:
  - `attrs(dict) -> str`, `trf_loci(rows) -> list[dict]`
  - `final_names(fam, arrays_by_genome, complex_name) -> dict[fam_id, name]`
  - `SUMMARY_COLS`, `write_genome(d, names, outdir, p)`, which writes `<G>_centromeric_repeat.gff`, `<G>_centromere_summary.tsv` and `<G>_hor_units.tsv`
  - `FAMTAB_COLS`, `write_families(fam, data, names, hmm_path, fa_path, outdir)`, which writes `satellite_families.{tsv,fa,hmm}`

- [ ] **Step 1: Write the GFF validator `tests/gffcheck.py`**

```python
def check_gff(path, chroms: dict) -> list:
    """Assert GFF3 validity: 9 columns, coordinates in range, unique IDs, Parents defined earlier and
    containing their children. Returns [(chrom, start, end, type, attrs)]."""
    ids, feats = {}, []
    with open(path) as f:
        assert f.readline().startswith("##gff-version 3")
        for line in f:
            if line.startswith("#"):
                continue
            x = line.rstrip("\n").split("\t")
            assert len(x) == 9, line
            c, s, e = x[0], int(x[3]), int(x[4])
            assert c in chroms and 1 <= s <= e <= chroms[c], line
            a = dict(kv.split("=", 1) for kv in x[8].split(";"))
            assert a["ID"] not in ids, a["ID"]
            if "Parent" in a:
                pc, ps, pe = ids[a["Parent"]]
                assert pc == c and ps <= s and e <= pe, line
            ids[a["ID"]] = (c, s, e)
            feats.append((c, s, e, x[2], a))
    return feats
```

- [ ] **Step 2: Write the failing tests** — `tests/test_output.py`

```python
import numpy as np
import pandas as pd
import pysam

from cenrep.annotate import ARR_COLS, MON_COLS
from cenrep.centromere import finalize
from cenrep.classify import CLS_COLS
from cenrep.genomes import Genome
from cenrep.hor import MONHOR_COLS, OCC_COLS, UNIT_COLS
from cenrep.longunit import LU_COLS
from cenrep.output import attrs, final_names, trf_loci, write_genome
from cenrep.params import Params
from cenrep.trf import COLS as TRF_COLS
from gffcheck import check_gff


def test_attrs_escaping():
    assert attrs({"ID": "a", "flags": "gaps:1,x=y;z", "none": None}) == "ID=a;flags=gaps:1,x%3Dy%3Bz"


def test_trf_loci_merge_and_top_row():
    rows = pd.DataFrame([("c", 0, 100, 2, 50.0, 2, 90, 0, 50, 1.0, "AT"),
                         ("c", 50, 300, 4, 60.0, 4, 90, 0, 80, 1.5, "ATGC"),
                         ("c", 400, 500, 10, 10.0, 10, 90, 0, 40, 1.9, "ACGTACGTAA")], columns=TRF_COLS)
    loci = trf_loci(rows)
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
    trf = pd.DataFrame([("Sp_A01", 5000, 5600, 150, 4.0, 150, 90, 1, 500, 1.9, "A" * 150),
                        ("Sp_A01", 6500, 7100, 2, 300.0, 2, 95, 0, 900, 1.0, "AT")], columns=TRF_COLS)
    mh = pd.DataFrame([(i, i % 2, "HOR_Sp_A01_1" if i in (0, 1) else ".", "hor_Sp_A01_1" if i in (0, 1) else ".")
                       for i in range(7)], columns=MONHOR_COLS)
    occ = pd.DataFrame([("hor_Sp_A01_1", "HOR_Sp_A01_1", "Sp_A01", 5000, 5300, "+", 0, 2, "HOR_Sp_A01_1.r1", 1)],
                       columns=OCC_COLS)
    units = pd.DataFrame([("HOR_Sp_A01_1", "fam001", "Sp_A01", 2, 2, 1, "0-1")], columns=UNIT_COLS)
    return dict(genome=g, arrays=arrays, monomers=mon, fams=_fams([("fam001", "centromeric"), ("fam002", "cen_associated")]),
                regions=rows, longunit=pd.DataFrame(columns=LU_COLS), mon_hor=mh, hor_occ=occ, hor_units=units,
                trf=trf, genes=None)


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
    assert tr[4]["period"] == "2" and tr[4]["Parent"] == "cen_Sp_A01"
    s = pd.read_csv(tmp_path / "Sp_A_centromere_summary.tsv", sep="\t", dtype=str)
    assert s.set_index("chrom").loc["Sp_A02", "confidence"] == "none"
    assert s.set_index("chrom").loc["Sp_A01", "cen_start"] == "5001"


def test_write_genome_without_regions(tmp_path):
    d = _data(tmp_path, with_region=False)
    d["fams"] = _fams([("fam001", "other"), ("fam002", "other")])
    write_genome(d, {"fam001": "Sp_sat01_150", "fam002": "Sp_sat02_100"}, tmp_path, Params())
    assert check_gff(tmp_path / "Sp_A_centromeric_repeat.gff", d["genome"].chroms) == []
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `python -m pytest -q tests/test_output.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'cenrep.output'`

- [ ] **Step 4: Write `cenrep/output.py`**

```python
"""S7: final family names, per-genome GFF3 and summary tables, complex-level family files."""
import json
from dataclasses import asdict
from datetime import date
from pathlib import Path

import pandas as pd
import pyhmmer

from . import __version__
from .util import clipped_bp, overlap, parse_fasta_text, union_bp, write_fasta

SUMMARY_COLS = ["chrom", "chrom_len", "cen_start", "cen_end", "cen_len", "evidence", "confidence", "flags",
                "main_family", "main_family_bp", "cen_sat_bp", "sat_fraction", "n_arrays", "largest_array_bp",
                "n_monomers", "mean_identity", "hor_monomers", "hor_fraction", "top_tandem_hor", "n_gaps",
                "gap_bp", "n_genes", "gene_ratio", "other_tandem_bp"]
FAMTAB_COLS = ["family", "fam_id", "superfamily", "consensus_len", "hmm_mode", "genome", "class", "strength",
               "defines", "bp", "n_arrays", "prevalence", "concentration", "subtel_frac", "gene_ratio",
               "median_main_bp", "eligible", "reason"]
LEVEL = {"centromere": 0, "satellite_DNA": 1, "tandem_repeat": 1, "higher_order_repeat": 2, "satellite_monomer": 3}
_ESC = {"%": "%25", ";": "%3B", "=": "%3D", "&": "%26", "\t": "%09", "\n": "%0A", "\r": "%0D"}


def _esc(v) -> str:
    return "".join(_ESC.get(ch, ch) for ch in str(v))


def attrs(d: dict) -> str:
    """GFF3 column 9; commas are kept (multi-value lists), reserved characters percent-encoded."""
    return ";".join(f"{k}={_esc(v)}" for k, v in d.items() if v is not None)


def trf_loci(rows: pd.DataFrame) -> list:
    """Merge overlapping TRF rows; a locus takes period/copies/pct_match/consensus of its top-scoring row."""
    loci, cur = [], None
    for r in rows.sort_values("start").itertuples(index=False):
        if cur is not None and r.start < cur["end"]:
            cur["end"] = max(cur["end"], int(r.end))
            if r.score > cur["score"]:
                cur.update(period=int(r.period), copies=r.copies, pct_match=int(r.pct_match),
                           consensus=r.consensus, score=r.score)
        else:
            if cur is not None:
                loci.append(cur)
            cur = dict(start=int(r.start), end=int(r.end), period=int(r.period), copies=r.copies,
                       pct_match=int(r.pct_match), consensus=r.consensus, score=r.score)
    if cur is not None:
        loci.append(cur)
    return loci


def final_names(fam: pd.DataFrame, arrays_by_genome: dict, complex_name: str) -> dict:
    bp = {f: 0 for f in fam.fam_id}
    for a in arrays_by_genome.values():
        for (f, _), d in a.groupby(["family", "chrom"]):
            if f in bp:
                bp[f] += union_bp(zip(d.start, d.end))
    order = sorted(fam.itertuples(index=False), key=lambda r: (-bp[r.fam_id], r.fam_id))
    return {r.fam_id: f"{complex_name}_sat{i:02d}_{r.consensus_len}" for i, r in enumerate(order, 1)}


def write_genome(d: dict, names: dict, outdir, p) -> None:
    outdir = Path(outdir)
    g, regions, fams = d["genome"], d["regions"], d["fams"]
    cls = dict(zip(fams.family, fams.cls))
    primary = {f for f, c in cls.items() if c == "centromeric"}
    reg = {c: r for c, r in regions.items() if r["evidence"] != "none"}
    order = {c: i for i, c in enumerate(g.chroms)}
    nm = lambda f: names.get(f, f)

    mon = d["monomers"].merge(d["mon_hor"], on="mon_idx", how="left")
    mon["variant_class"] = ["." if pd.isna(v) else str(int(v)) for v in mon.variant_class]
    mon["hor_unit"] = mon.hor_unit.fillna(".")
    mon["classed"], mon["inhor"] = mon.variant_class != ".", mon.hor_unit != "."
    per_arr = mon.groupby("array_id").agg(classed=("classed", "sum"), inhor=("inhor", "sum"))

    arrays = d["arrays"].copy()
    ov, inside = [], []
    for r in arrays.itertuples(index=False):
        rg = reg.get(r.chrom)
        ov.append(overlap(r.start, r.end, rg["cen_start"], rg["cen_end"]) if rg else 0)
        inside.append(bool(rg) and r.start >= rg["cen_start"] and r.end <= rg["cen_end"])
    arrays["ov"], arrays["inside"] = ov, inside
    inc = arrays[arrays.family.isin(primary) | (arrays.ov > 0)].copy()
    inc = inc.assign(ci=inc.chrom.map(order)).sort_values(["ci", "start", "end"])
    inc["gid"] = [f"arr_{c}_{k}" for c, k in zip(inc.chrom, inc.groupby("chrom").cumcount() + 1)]
    gid = dict(zip(inc.array_id, inc.gid))

    stats, tandem = {}, []
    for c, rg in reg.items():
        s0, e0 = rg["cen_start"], rg["cen_end"]
        a = inc[(inc.chrom == c) & inc.inside & inc.family.isin(primary)]
        m = mon[mon.array_id.isin(a.array_id)]
        classed, inhor = int(m.classed.sum()), int(m.inhor.sum())
        units = d["hor_units"][d["hor_units"].chrom == c]
        top = "."
        if len(units) and units.max_tandem_copies.max() > 1:
            u = units.sort_values(["max_tandem_copies", "n_monomers"], ascending=False).iloc[0]
            top = f"{u.unit}:{u.n_monomers}x{u.max_tandem_copies}"
        stats[c] = dict(cen_sat_bp=union_bp(zip(a.start, a.end)), n_arrays=len(a),
                        largest_array_bp=int((a.end - a.start).max()) if len(a) else 0, n_monomers=len(m),
                        mean_identity=round(float(m.identity.mean()), 3) if len(m) else ".",
                        hor_monomers=inhor, hor_fraction=round(inhor / classed, 3) if classed else ".",
                        top_tandem_hor=top)
        sat_iv = list(zip(inc.start[inc.chrom == c], inc.end[inc.chrom == c]))
        cand = []
        t = d["trf"]
        t = t[(t.chrom == c) & (t.start < e0) & (t.end > s0)]
        # drop satellite-derived TRF rows BEFORE merging, so a satellite row that bleeds a few bp into an
        # adjacent microsatellite cannot swallow it into one locus
        t = t.loc[[clipped_bp(sat_iv, s, e) < 0.5 * (e - s) for s, e in zip(t.start, t.end)]]
        for L in trf_loci(t):
            cand.append((L["start"], L["end"], dict(period=L["period"], copies=L["copies"],
                                                    pct_match=L["pct_match"], unit_source="trf",
                                                    consensus=L["consensus"])))
        lu = d["longunit"]
        for r in lu[lu.chrom == c].itertuples(index=False):
            cand.append((int(r.start), int(r.end), dict(period=int(r.unit), copies=r.copies, identity=r.identity,
                                                         unit_source="selfaln")))
        k = 0
        for s, e, a_ in sorted(cand, key=lambda x: (x[0], x[1])):
            if e - s < p.tandem_min_len or s < s0 or e > e0 or clipped_bp(sat_iv, s, e) >= 0.5 * (e - s):
                continue
            k += 1
            if len(a_.get("consensus", "")) > 100:
                a_["consensus_len"] = len(a_["consensus"])
                a_["consensus"] = a_["consensus"][:100] + "..."
            tandem.append((c, s, e, dict(ID=f"tr_{c}_{k}", Parent=f"cen_{c}", **a_)))

    feats = []
    for c, rg in reg.items():
        st = stats[c]
        feats.append((c, rg["cen_start"], rg["cen_end"], "centromere", ".", ".", dict(
            ID=f"cen_{c}", evidence=rg["evidence"], confidence=rg["confidence"], flags=rg["flags"],
            main_family=nm(rg["main_family"]), satellite_bp=st["cen_sat_bp"],
            sat_fraction=round(st["cen_sat_bp"] / (rg["cen_end"] - rg["cen_start"]), 3),
            n_arrays=st["n_arrays"], hor_fraction=st["hor_fraction"], n_genes=rg["n_genes"], gaps=rg["n_gaps"])))
    for r in inc.itertuples(index=False):
        hf = "."
        if r.array_id in per_arr.index and per_arr.loc[r.array_id, "classed"]:
            hf = round(per_arr.loc[r.array_id, "inhor"] / per_arr.loc[r.array_id, "classed"], 3)
        a_ = dict(ID=r.gid)
        if r.inside:
            a_["Parent"] = f"cen_{r.chrom}"
        a_.update(family=nm(r.family), family_class=cls.get(r.family, "other"), monomer_len=r.monomer_len,
                  n_monomers=r.n_monomers, mean_identity=round(r.mean_identity, 3), hor_fraction=hf,
                  in_centromere="true" if r.ov >= 0.5 * (r.end - r.start) else "false")
        feats.append((r.chrom, r.start, r.end, "satellite_DNA", ".", r.strand, a_))
    m = mon[mon.array_id.isin(inc.array_id)]
    m = m.assign(ci=m.chrom.map(order)).sort_values(["ci", "start"])
    m["k"] = m.groupby("chrom").cumcount() + 1
    for r in m.itertuples(index=False):
        feats.append((r.chrom, r.start, r.end, "satellite_monomer", f"{r.bitscore:.1f}", r.strand, dict(
            ID=f"mon_{r.chrom}_{r.k}", Parent=gid[r.array_id], family=nm(r.family), identity=round(r.identity, 3),
            model_cov=round((r.hmm_to - r.hmm_from + 1) / r.model_len, 3), variant_class=r.variant_class,
            hor=r.hor_unit)))
    uocc = dict(zip(d["hor_units"].unit, d["hor_units"].occurrences))
    for r in d["hor_occ"].itertuples(index=False):
        if r.array_id in gid:
            feats.append((r.chrom, r.start, r.end, "higher_order_repeat", ".", r.strand, dict(
                ID=r.occ_id, Parent=gid[r.array_id], unit=r.unit, unit_monomers=r.n_monomers,
                unit_occurrences=uocc.get(r.unit, "."), tandem_copies=r.tandem_copies)))
    for c, s, e, a_ in tandem:
        feats.append((c, s, e, "tandem_repeat", ".", ".", a_))
    feats.sort(key=lambda f: (order[f[0]], f[1], LEVEL[f[3]], -f[2]))

    with open(outdir / f"{g.name}_centromeric_repeat.gff", "w") as f:
        f.write("##gff-version 3\n")
        for c, L in g.chroms.items():
            f.write(f"##sequence-region {c} 1 {L}\n")
        f.write(f"# cenrep {__version__} {date.today()} genome={g.name} "
                f"params={json.dumps(asdict(p), separators=(',', ':'))}\n")
        for c, s, e, t, sc, st, a_ in feats:
            f.write(f"{c}\tcenrep\t{t}\t{s + 1}\t{e}\t{sc}\t{st}\t.\t{attrs(a_)}\n")

    rows = []
    for c, L in g.chroms.items():
        rg = regions.get(c)
        if rg is None or rg["evidence"] == "none":
            rows.append({**{k: "." for k in SUMMARY_COLS}, "chrom": c, "chrom_len": L,
                         "evidence": "none", "confidence": "none"})
            continue
        st, s0, e0 = stats[c], rg["cen_start"], rg["cen_end"]
        oth = inc[(inc.chrom == c) & inc.inside & ~inc.family.isin(primary)]
        other_iv = [(s, e) for cc, s, e, _ in tandem if cc == c] + list(zip(oth.start, oth.end))
        rows.append(dict(chrom=c, chrom_len=L, cen_start=s0 + 1, cen_end=e0, cen_len=e0 - s0,
                         evidence=rg["evidence"], confidence=rg["confidence"], flags=rg["flags"],
                         main_family=nm(rg["main_family"]), main_family_bp=rg["main_family_bp"],
                         cen_sat_bp=st["cen_sat_bp"], sat_fraction=round(st["cen_sat_bp"] / (e0 - s0), 3),
                         n_arrays=st["n_arrays"], largest_array_bp=st["largest_array_bp"],
                         n_monomers=st["n_monomers"], mean_identity=st["mean_identity"],
                         hor_monomers=st["hor_monomers"], hor_fraction=st["hor_fraction"],
                         top_tandem_hor=st["top_tandem_hor"], n_gaps=rg["n_gaps"], gap_bp=rg["gap_bp"],
                         n_genes=rg["n_genes"], gene_ratio=rg["gene_ratio"], other_tandem_bp=union_bp(other_iv)))
    pd.DataFrame(rows, columns=SUMMARY_COLS).to_csv(outdir / f"{g.name}_centromere_summary.tsv", sep="\t", index=False)
    hu = d["hor_units"].copy()
    hu["family"] = hu.family.map(nm)
    hu.to_csv(outdir / f"{g.name}_hor_units.tsv", sep="\t", index=False)


def write_families(fam: pd.DataFrame, data: dict, names: dict, hmm_path, fa_path, outdir) -> None:
    outdir = Path(outdir)
    rows = []
    for gname, d in data.items():
        f = d["fams"].set_index("family")
        for r in fam.itertuples(index=False):
            row = dict(family=names[r.fam_id], fam_id=r.fam_id, superfamily=names.get(r.superfamily, r.superfamily),
                       consensus_len=r.consensus_len, hmm_mode=r.hmm_mode, genome=gname)
            if r.fam_id in f.index:
                x = f.loc[r.fam_id]
                row.update({"class": x.cls, "strength": x.strength, "defines": x.defines, "bp": x.bp,
                            "n_arrays": x.n_arrays, "prevalence": x.prevalence, "concentration": x.concentration,
                            "subtel_frac": x.subtel_frac, "gene_ratio": x.gene_ratio,
                            "median_main_bp": x.median_main_bp, "eligible": x.eligible, "reason": x.reason})
            else:
                row.update({"class": "absent", "bp": 0, "n_arrays": 0, "reason": "no arrays in this genome"})
            rows.append(row)
    pd.DataFrame(rows, columns=FAMTAB_COLS).to_csv(outdir / "satellite_families.tsv", sep="\t", index=False,
                                                   float_format="%.3f", na_rep=".")
    cons = dict(parse_fasta_text(Path(fa_path).read_text())) if Path(fa_path).exists() else {}
    write_fasta(outdir / "satellite_families.fa",
                ((f"{names[i]} fam_id={i} hmm={m}", cons.get(i, "")) for i, m in zip(fam.fam_id, fam.hmm_mode)))
    hmms = []
    if len(fam):
        with pyhmmer.plan7.HMMFile(hmm_path) as hf:
            hmms = list(hf)
    with open(outdir / "satellite_families.hmm", "wb") as fh:
        for h in hmms:
            h.name = names[h.name]
            h.write(fh)
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest -q tests/test_output.py`
Expected: PASS (5 tests).

---

### Task 10: S8 — overview figure (`cenrep/plot.py`)

**Files:**
- Create: `clean/cenrep/cenrep/plot.py`
- Test: `clean/cenrep/tests/test_plot.py`

**Interfaces:**
- Consumes: `panels`, a list of 4 dicts with keys:
  - `genome`, `chroms`
  - `arrays` (ARR_COLS; `family` holds final names; extra column `cls`)
  - `regions` (`dict[chrom, row]`), `genes` (`dict|None`)
- Produces: `plot_overview(path, panels)` writes a vector PDF; `bin_bp(arrays, L, bin_size) -> ndarray`

- [ ] **Step 0: Load the `dataviz` skill** (Skill tool) before writing the plot code. Keep this task's CLAUDE.md figure rules: no titles, panel letters A–D in the upper-left, colourblind-safe, 7 pt, PDF. Adjust colours or marks only if the skill's validator flags them.

- [ ] **Step 1: Write the failing tests** — `tests/test_plot.py`

```python
import numpy as np
import pandas as pd

from cenrep.annotate import ARR_COLS
from cenrep.plot import bin_bp, plot_overview


def _panel(name, with_fam=True):
    chroms = {f"{name}01": 3_000_000, f"{name}02": 2_000_000}
    rows = [(0, f"{name}01", 1_000_000, 1_250_000, "+", "Cx_sat01_178", 1400, 0.9, 178)] if with_fam else []
    a = pd.DataFrame(rows, columns=ARR_COLS).assign(cls="centromeric")
    reg = {f"{name}01": dict(evidence="satellite", cen_start=1_000_000, cen_end=1_250_000),
           f"{name}02": dict(evidence="none")} if with_fam else {c: dict(evidence="none") for c in chroms}
    genes = {f"{name}01": np.arange(0, 3_000_000, 20_000, dtype=np.int64)}
    return dict(genome=name, chroms=chroms, arrays=a, regions=reg, genes=genes)


def test_bin_bp_splits_across_bins():
    a = pd.DataFrame({"start": [50_000], "end": [250_000]})
    assert bin_bp(a, 300_000, 100_000).tolist() == [50_000, 100_000, 50_000]


def test_plot_writes_pdf(tmp_path):
    out = tmp_path / "o.pdf"
    plot_overview(out, [_panel(n) for n in ("Aa_A", "Ab_B", "Ac_A", "Ad_B")])
    assert out.stat().st_size > 5_000


def test_plot_without_families(tmp_path):
    out = tmp_path / "o.pdf"
    plot_overview(out, [_panel(n, with_fam=False) for n in ("Aa_A", "Ab_B", "Ac_A", "Ad_B")])
    assert out.exists()
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest -q tests/test_plot.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'cenrep.plot'`

- [ ] **Step 3: Write `cenrep/plot.py`**

```python
"""S8: centromere overview: one panel per genome (A-D), one track per chromosome.
Track = satellite bp per 100-kb bin (primary families coloured, cen-associated grey), scaled gene density line,
and a box for the centromere call (dashed for fallback calls). No titles; panel letters upper-left."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"]
ASSOC = "#999999"
BIN = 100_000
H = 0.8


def bin_bp(arrays, L: int, bin_size: int = BIN) -> np.ndarray:
    out = np.zeros(max(1, int(np.ceil(L / bin_size))))
    for s, e in zip(arrays.start, arrays.end):
        for b in range(int(s) // bin_size, (int(e) - 1) // bin_size + 1):
            out[b] += min(e, (b + 1) * bin_size) - max(s, b * bin_size)
    return out


def plot_overview(path, panels: list) -> None:
    plt.rcParams.update({"font.size": 7, "pdf.fonttype": 42, "ps.fonttype": 42})
    fams = sorted({f for pn in panels for f in pn["arrays"].loc[pn["arrays"].cls == "centromeric", "family"]})
    color = {f: OKABE_ITO[i % len(OKABE_ITO)] for i, f in enumerate(fams)}
    nmax = max(len(pn["chroms"]) for pn in panels)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 1.2 + 0.32 * nmax * 2), constrained_layout=True)
    for ax, pn, letter in zip(axes.flat, panels, "ABCD"):
        chroms = list(pn["chroms"].items())
        n = len(chroms)
        for row, (c, L) in enumerate(chroms):
            y = n - 1 - row
            ax.add_patch(Rectangle((0, y), L / 1e6, H, fc="#EEEEEE", ec="none", zorder=0))
            a = pn["arrays"][pn["arrays"].chrom == c]
            x = (np.arange(int(np.ceil(L / BIN))) + 0.5) * BIN / 1e6
            bottom = np.zeros(len(x))
            layers = [(a[(a.family == f) & (a.cls == "centromeric")], color[f]) for f in fams]
            layers.append((a[a.cls == "cen_associated"], ASSOC))
            for sel, col in layers:
                if sel.empty:
                    continue
                h = bin_bp(sel, L)[:len(x)] / BIN * H
                ax.bar(x, h, width=BIN / 1e6, bottom=y + bottom, color=col, linewidth=0, zorder=2)
                bottom += h
            genes = pn["genes"]
            if genes is not None and len(genes.get(c, [])):
                cnt = np.bincount(genes[c] // BIN, minlength=len(x))[:len(x)]
                ax.plot(x, y + cnt / max(cnt.max(), 1) * H, color="black", lw=0.4, zorder=3)
            r = pn["regions"].get(c)
            if r is not None and r["evidence"] != "none":
                ax.add_patch(Rectangle((r["cen_start"] / 1e6, y - 0.06), (r["cen_end"] - r["cen_start"]) / 1e6,
                                       H + 0.12, fill=False, ec="black", lw=0.7,
                                       ls="-" if r["evidence"] == "satellite" else "--", zorder=4))
        ax.set_yticks([n - 1 - i + H / 2 for i in range(n)])
        ax.set_yticklabels([c for c, _ in chroms])
        ax.set_ylim(-0.3, n + 0.1)
        ax.set_xlim(0, max(L for _, L in chroms) / 1e6 * 1.01)
        ax.set_xlabel("Position (Mb)")
        ax.set_ylabel(pn["genome"])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.text(-0.2, 1.02, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")
    handles = [Patch(color=color[f], label=f) for f in fams]
    handles += [Patch(color=ASSOC, label="cen-associated satellite"),
                Line2D([], [], color="black", lw=0.6, label="gene density (scaled)"),
                Patch(fill=False, ec="black", label="centromere call (dashed: fallback)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    fig.savefig(path, dpi=300)
    plt.close(fig)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest -q tests/test_plot.py`
Expected: PASS (3 tests).

---

### Task 11: Pipeline orchestration, CLI, end-to-end test

**Files:**
- Create: `clean/cenrep/cenrep/pipeline.py`, `clean/cenrep/cenrep/cli.py`
- Modify: `clean/cenrep/cenrep/__main__.py` (replace the stub)
- Test: `clean/cenrep/tests/test_e2e.py`

**Interfaces:**
- Consumes: every stage module above
- Produces:
  - `pipeline.STAGES`, `pipeline.TOOLS`, `Pipeline(complex_dir, threads, trf_cache, p, force=None).run()`
  - `cli.main(argv=None) -> int`; `python -m cenrep run ...`

- [ ] **Step 1: Write the failing end-to-end tests** — `tests/test_e2e.py`

```python
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


def test_resume_and_param_change(synth_run):
    from cenrep.pipeline import Pipeline
    cdir, _, root = synth_run
    Pipeline(cdir, 2, root / "cache", Params(**TEST_PARAMS)).run()          # all checkpoints: no-op
    changed = Params(**dict(TEST_PARAMS, cluster_gap=60_000))
    with pytest.raises(CenrepError, match="parameters differ"):
        Pipeline(cdir, 2, root / "cache", changed)
    Pipeline(cdir, 2, root / "cache", changed, force="plot")                 # explicit force is accepted


def test_cli_help():
    r = subprocess.run([sys.executable, "-m", "cenrep", "run", "--help"], capture_output=True, text=True)
    assert r.returncode == 0 and "--cluster-gap" in r.stdout and "--complex-dir" in r.stdout
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python -m pytest -q tests/test_e2e.py`
Expected: FAIL with `ModuleNotFoundError: No module named 'cenrep.pipeline'`

- [ ] **Step 3: Write `cenrep/pipeline.py`**

```python
"""Stage orchestration. Checkpoints: <complex>/centromere/work/<stage>.done; a rerun resumes by default."""
import json
import time
from pathlib import Path

import pandas as pd
import pyhmmer

from . import __version__, annotate, centromere, classify, families, hor, output, plot, trf
from .genomes import discover_genomes, load_genes
from .params import Params, to_dict
from .util import CenrepError, log, require_tools, touch

STAGES = ["trf", "families", "annotate", "classify", "centromeres", "hor", "output", "plot"]
TOOLS = ["trf", "makeblastdb", "blastn", "mafft", "cd-hit-est", "minimap2"]


class Pipeline:
    def __init__(self, complex_dir, threads: int, trf_cache, p: Params, force: str = None):
        self.cdir = Path(complex_dir).resolve()
        self.name = self.cdir.name
        self.out = self.cdir / "centromere"
        self.work = self.out / "work"
        self.threads, self.cache, self.p = threads, Path(trf_cache).resolve(), p
        self.genomes = discover_genomes(self.cdir)
        self.work.mkdir(parents=True, exist_ok=True)
        pj, cur = self.work / "params.json", to_dict(p)
        done_any = any((self.work / f"{s}.done").exists() for s in STAGES)
        if not force and done_any and pj.exists() and json.loads(pj.read_text()) != cur:
            raise CenrepError(f"parameters differ from the checkpointed run in {self.work}; rerun with "
                              f"--force <earliest stage affected by the change>")
        if force:
            for s in STAGES[STAGES.index(force):]:
                (self.work / f"{s}.done").unlink(missing_ok=True)
        pj.write_text(json.dumps(cur, indent=1))
        self._genes = {}

    # ---------------------------------------------------------------- helpers
    def _p(self, stage: str, name: str) -> Path:
        return self.work / stage / name

    def genes(self, g):
        if g.name not in self._genes:
            self._genes[g.name] = load_genes(self.cdir, g.name)
        return self._genes[g.name]

    def trf_table(self, g):
        return trf.load_trf(self._p("trf", f"{g.name}.trf.tsv.gz"))

    def fam_table(self):
        return pd.read_csv(self._p("families", "families.tsv"), sep="\t", dtype={"superfamily": str})

    def monomers(self, g):
        return pd.read_csv(self._p("annotate", f"{g.name}.monomers.tsv.gz"), sep="\t", dtype={"chrom": str})

    def arrays(self, g):
        return pd.read_csv(self._p("annotate", f"{g.name}.arrays.tsv"), sep="\t", dtype={"chrom": str})

    def fams(self, g, stage: str = "centromeres"):
        return pd.read_csv(self._p(stage, f"{g.name}.families.tsv"), sep="\t")

    def regions(self, g):
        return centromere.load_regions(self._p("centromeres", f"{g.name}.regions.tsv"))

    # ---------------------------------------------------------------- driver
    def run(self) -> None:
        require_tools(TOOLS)
        t0 = time.time()
        log.info(f"cenrep {__version__}: {self.name} [{', '.join(g.name for g in self.genomes)}], "
                 f"{self.threads} threads")
        for s in STAGES:
            if (self.work / f"{s}.done").exists():
                log.info(f"{s}: checkpoint found, skipping")
                continue
            t = time.time()
            log.info(f"{s}: start")
            (self.work / s).mkdir(exist_ok=True)
            getattr(self, f"stage_{s}")()
            touch(self.work / f"{s}.done")
            log.info(f"{s}: done ({time.time() - t:.0f} s)")
        log.info(f"all stages done ({(time.time() - t0) / 60:.1f} min); outputs in {self.out}")

    # ---------------------------------------------------------------- stages
    def stage_trf(self):
        paths = trf.trf_genomes(self.genomes, self.cache, self.threads, self.p)
        for g in self.genomes:
            link = self._p("trf", f"{g.name}.trf.tsv.gz")
            link.unlink(missing_ok=True)
            link.symlink_to(paths[g.name])

    def stage_families(self):
        families.build_families({g.name: self.trf_table(g) for g in self.genomes}, self.genomes,
                                self.work / "families", self.threads, self.p)

    def stage_annotate(self):
        fam = self.fam_table()
        hmms = []
        if len(fam):
            with pyhmmer.plan7.HMMFile(self._p("families", "families.hmm")) as f:
                hmms = list(f)
        fam_len = dict(zip(fam.fam_id, fam.consensus_len))
        for g in self.genomes:
            annotate.annotate_genome(g, hmms, fam_len, self.threads, self.p,
                                     self._p("annotate", f"{g.name}.monomers.tsv.gz"),
                                     self._p("annotate", f"{g.name}.arrays.tsv"))

    def stage_classify(self):
        for g in self.genomes:
            fams, mains = classify.family_metrics(self.arrays(g), g.chroms, self.genes(g), self.p)
            fams, per_chrom = classify.choose_centromeric(fams, mains, g.chroms, self.p)
            fams.to_csv(self._p("classify", f"{g.name}.families.tsv"), sep="\t", index=False)
            prim = fams[fams.cls == "centromeric"]
            log.info(f"classify: {g.name}: centromeric = {', '.join(prim.family) or 'none'}; "
                     f"{len(per_chrom)}/{len(g.chroms)} chromosomes have one")

    def stage_centromeres(self):
        for g in self.genomes:
            arrays = self.arrays(g)
            fams = self.fams(g, "classify")
            primary = set(fams.family[fams.cls == "centromeric"])
            rows, lu = centromere.call_regions(g, arrays, primary, self.genes(g), self.trf_table(g),
                                               self.work / "centromeres" / f"tmp_{g.name}", self.threads, self.p)
            fams = classify.assign_cen_associated(fams, arrays, rows, self.p)
            centromere.write_regions(rows, self._p("centromeres", f"{g.name}.regions.tsv"))
            lu.to_csv(self._p("centromeres", f"{g.name}.longunit.tsv"), sep="\t", index=False)
            fams.to_csv(self._p("centromeres", f"{g.name}.families.tsv"), sep="\t", index=False)
            conf = pd.Series([r["confidence"] for r in rows.values()]).value_counts().to_dict()
            log.info(f"centromeres: {g.name}: {conf}; {len(lu)} long-unit arrays")

    def stage_hor(self):
        for g in self.genomes:
            fams = self.fams(g)
            primary = set(fams.family[fams.cls == "centromeric"])
            mh, occ, units = hor.run_hor(g, self.monomers(g), self.arrays(g), self.regions(g), primary,
                                         self.work / "hor" / f"tmp_{g.name}", self.threads, self.p)
            mh.to_csv(self._p("hor", f"{g.name}.monomer_hor.tsv.gz"), sep="\t", index=False)
            occ.to_csv(self._p("hor", f"{g.name}.hor_occ.tsv"), sep="\t", index=False)
            units.to_csv(self._p("hor", f"{g.name}.hor_units.tsv"), sep="\t", index=False)
            log.info(f"hor: {g.name}: {len(units)} HOR units, {len(occ)} occurrences, "
                     f"{(mh.hor_unit != '.').sum()}/{len(mh)} classed monomers in HORs")

    def genome_data(self, g) -> dict:
        rd = lambda n: pd.read_csv(self._p("hor", n), sep="\t", dtype={"chrom": str})
        return dict(genome=g, arrays=self.arrays(g), monomers=self.monomers(g), fams=self.fams(g),
                    regions=self.regions(g),
                    longunit=pd.read_csv(self._p("centromeres", f"{g.name}.longunit.tsv"), sep="\t",
                                         dtype={"chrom": str}),
                    mon_hor=rd(f"{g.name}.monomer_hor.tsv.gz"), hor_occ=rd(f"{g.name}.hor_occ.tsv"),
                    hor_units=rd(f"{g.name}.hor_units.tsv"), trf=self.trf_table(g), genes=self.genes(g))

    def stage_output(self):
        fam = self.fam_table()
        data = {g.name: self.genome_data(g) for g in self.genomes}
        names = output.final_names(fam, {k: d["arrays"] for k, d in data.items()}, self.name)
        pd.DataFrame(sorted(names.items()), columns=["fam_id", "name"]).to_csv(
            self._p("output", "names.tsv"), sep="\t", index=False)
        for g in self.genomes:
            output.write_genome(data[g.name], names, self.out, self.p)
        output.write_families(fam, data, names, self._p("families", "families.hmm"),
                              self._p("families", "families.fa"), self.out)

    def stage_plot(self):
        names = dict(pd.read_csv(self._p("output", "names.tsv"), sep="\t").itertuples(index=False))
        panels = []
        for g in self.genomes:
            fams = self.fams(g)
            a = self.arrays(g)
            a = a.assign(cls=a.family.map(dict(zip(fams.family, fams.cls))).fillna("other"),
                         family=a.family.map(names))
            panels.append(dict(genome=g.name, chroms=g.chroms, arrays=a, regions=self.regions(g),
                               genes=self.genes(g)))
        plot.plot_overview(self.out / "centromere_overview.pdf", panels)
```

- [ ] **Step 4: Write `cenrep/cli.py` and replace `cenrep/__main__.py`**

`cenrep/cli.py`:
```python
"""Command line: cenrep run --complex-dir DIR --threads N --trf-cache DIR [--force STAGE] [-v] [thresholds]"""
import argparse
from dataclasses import fields
from pathlib import Path

from . import __version__
from .params import HELP, Params
from .util import CenrepError, log, setup_logging


def main(argv=None) -> int:
    from .pipeline import STAGES, Pipeline
    ap = argparse.ArgumentParser(prog="cenrep", description=(
        "Centromeric satellite, HOR and tandem-repeat annotation for one complex dir holding 4 genomes. "
        "Outputs go to <complex-dir>/centromere/."))
    ap.add_argument("--version", action="version", version=f"cenrep {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run all stages; resumes from checkpoints in centromere/work/")
    r.add_argument("--complex-dir", required=True, type=Path, help="dir with exactly 4 *.fa (+ .fai)")
    r.add_argument("--threads", type=int, default=8, help="max threads/processes [8]")
    r.add_argument("--trf-cache", required=True, type=Path, help="shared TRF cache dir (keyed by genome md5)")
    r.add_argument("--force", choices=STAGES, help="recompute this stage and all later ones")
    r.add_argument("-v", "--verbose", action="store_true", help="per-step progress and sanity checks")
    g = r.add_argument_group("thresholds")
    for f in fields(Params):
        g.add_argument("--" + f.name.replace("_", "-"), type=type(f.default), default=f.default,
                       help=f"{HELP[f.name]} [{f.default}]")
    a = ap.parse_args(argv)
    setup_logging(a.verbose)
    p = Params(**{f.name: getattr(a, f.name) for f in fields(Params)})
    try:
        Pipeline(a.complex_dir, a.threads, a.trf_cache, p, force=a.force).run()
    except CenrepError as e:
        log.error(str(e))
        return 1
    return 0
```

`cenrep/__main__.py`:
```python
from .cli import main

raise SystemExit(main())
```

- [ ] **Step 5: Run the end-to-end tests to verify they pass**

Run: `python -m pytest -q tests/test_e2e.py`
Expected: PASS (7 tests), in about 1–3 min. If a boundary or HOR assertion fails:
- Inspect `centromere/work/` of the tmp run. The path is printed by pytest `--basetemp`; rerun with `--basetemp=/tmp/claude-1007/.../scratchpad/e2e` to keep it.
- Fix the code, not the test's tolerances. A tolerance change is only acceptable if the planted truth was mis-specified; in that case fix `synth.py` and say why.

- [ ] **Step 6: Run the full suite**

Run: `python -m pytest -q`
Expected: every test passes, 0 failures, 0 errors.

---

### Task 12: Launcher `clean/run_centromere.sh`

**Files:**
- Create: `clean/run_centromere.sh`

**Interfaces:**
- Consumes: `python -m cenrep run` (Task 11)
- Produces: a detached job that writes `clean/centromere_run.log`, `clean/<Complex>/centromere/cenrep.{log,exit,time}`, and the TRF cache `clean/centromere_cache/trf/`

- [ ] **Step 1: Write the launcher**

```bash
#!/usr/bin/env bash
# Centromeric repeat annotation (cenrep) for the allo10 complexes.
# Runs complexes one after another at 32 threads, inside ONE detached job (setsid nohup; survives logout).
#   input : clean/<Complex>/*.fa (4 per complex); gene GFFs from clean/<Complex>/window_identity (optional)
#   output: clean/<Complex>/centromere/  (GFFs, summaries, figure, cenrep.log/.exit/.time)
#   cache : clean/centromere_cache/trf/  (TRF per genome md5; shared Brassica diploids run once)
# Usage : bash run_centromere.sh [Complex ...]     (default: every complex dir with .done)
#         DRY=1 bash run_centromere.sh [Complex ...] (validate + list, launch nothing)
# Status: cat centromere_run.log ; cat */centromere/cenrep.exit ; tail -n3 */centromere/cenrep.log
set -euo pipefail

THREADS=32
ENV=/home/chris/bin/mambaforge/envs/synLTR
CLEAN=/data2/chris/tetraploid/allo10/clean
SRC=$CLEAN/cenrep

die() { echo "ERROR: $*" >&2; exit 1; }
[[ -x $ENV/bin/python ]] || die "missing env python: $ENV/bin/python"
for t in trf makeblastdb blastn mafft cd-hit-est minimap2 samtools; do
  [[ -x $ENV/bin/$t ]] || die "missing $t in $ENV"
done
[[ -f $SRC/cenrep/cli.py ]] || die "missing cenrep source in $SRC"

cd "$CLEAN"
if (( $# )); then complexes=("$@"); else
  complexes=(); for d in */.done; do complexes+=("${d%/.done}"); done
fi
(( ${#complexes[@]} )) || die "no complexes found"
for c in "${complexes[@]}"; do
  shopt -s nullglob; fa=("$c"/*.fa); shopt -u nullglob
  (( ${#fa[@]} == 4 )) || die "$c: expected 4 *.fa, found ${#fa[@]}"
done

if [[ -n ${DRY:-} ]]; then printf 'would run: %s\n' "${complexes[@]}"; exit 0; fi

if [[ -z ${CENREP_DETACHED:-} ]]; then
  CENREP_DETACHED=1 setsid nohup bash "$0" "${complexes[@]}" > "$CLEAN/centromere_run.log" 2>&1 < /dev/null &
  echo "launched pid $! (log: $CLEAN/centromere_run.log)"
  exit 0
fi

export PATH=$ENV/bin:$PATH PYTHONPATH=$SRC
for c in "${complexes[@]}"; do
  out=$CLEAN/$c/centromere
  mkdir -p "$out"; rm -f "$out/cenrep.exit"
  for f in "$CLEAN/$c"/*.fa; do [[ -f $f.fai ]] || samtools faidx "$f"; done
  echo "$(date '+%F %T') start $c"
  set +e
  /usr/bin/time -v -o "$out/cenrep.time" python -m cenrep run --complex-dir "$CLEAN/$c" \
      --threads "$THREADS" --trf-cache "$CLEAN/centromere_cache/trf" > "$out/cenrep.log" 2>&1
  rc=$?
  set -e
  echo "$rc" > "$out/cenrep.exit"
  echo "$(date '+%F %T') done $c exit=$rc"
done
echo "$(date '+%F %T') all done"
```

- [ ] **Step 2: Syntax check and dry run**

Run: `bash -n /data2/chris/tetraploid/allo10/clean/run_centromere.sh && DRY=1 bash /data2/chris/tetraploid/allo10/clean/run_centromere.sh`
Expected: `would run:` followed by all 10 complex names, no errors.

Run: `DRY=1 bash /data2/chris/tetraploid/allo10/clean/run_centromere.sh NoSuch`
Expected: `ERROR: NoSuch: expected 4 *.fa, found 0` and a non-zero exit.

---

### Task 13: Real data — pilot on Asuecica, validate, then all complexes

**Files:**
- Create: `clean/memo_centromere.md`
- Modify (memory): `/home/chris/.claude/projects/-data2-chris-tetraploid-allo10-clean/memory/` (new project memory + `MEMORY.md` pointer)

- [ ] **Step 1: Launch the pilot (Asuecica, detached) and verify it is detached**

Run: `bash /data2/chris/tetraploid/allo10/clean/run_centromere.sh Asuecica`
Then: `ps -o pid,ppid,etime,cmd -p $(pgrep -f "run_centromere.sh Asuecica" | head -1)`
Expected: PPID 1. Poll `tail -n5 clean/Asuecica/centromere/cenrep.log` every ~10 min (Monitor with an until-loop, not sleep). Expect TRF ~20–40 min for the 4 × ~140 Mb genomes.

- [ ] **Step 2: Validate the pilot against the spec (§7.1) and the literature**

Check each of the following, and write the numbers into the memo:
1. `Asuecica/centromere/cenrep.exit` is `0`.
2. `satellite_families.tsv`:
   - Athaliana_T and Asuecica_T have a centromeric family of 176–180 bp. Its consensus must be ≥ 80% identical to the CEN178 template in the spec; check with `blastn -subject`.
   - Aarenosa_A and Asuecica_A have a ~178 bp pAa-type family. This is a separate family, because CEN180 vs pAa is 58–80%.
3. Athaliana_T:
   - the centromeric-family monomer count in the GFF (`grep -c satellite_monomer`, restricted to that family) is within ±10% of 66,131;
   - there are 5 centromere rows with `confidence` high or medium, each 1–6 Mb.
4. The AthCEN159 (~159 bp) family is present as `cen_associated`, or as a second `centromeric` family on some chromosomes. Record which.
5. `centromere_overview.pdf` shows one satellite block per chromosome, sitting in a gene-density dip.

If any check fails, STOP and diagnose with superpowers:systematic-debugging before launching the rest.

- [ ] **Step 3: Launch the remaining 9 complexes (one detached job)**

Run: `bash /data2/chris/tetraploid/allo10/clean/run_centromere.sh Akamchatica Bcarinata Bhybridum Bjuncea Bnapus Carabica Cbursapastoris Msookensis Ominuta`
Verify PPID 1 as in Step 1. Poll `clean/centromere_run.log`.

- [ ] **Step 4: Per-complex validation against the literature table**

For every complex:
- Check that `cenrep.exit` is 0, and read `satellite_families.tsv`, the four summaries and the PDF.
- Compare against the expected table in spec §7.2:

  | complex | what to check |
  |---|---|
  | Akamchatica | chromosome-specific AlyCEN168/179 split |
  | Brassica A/C | CentBr1/2 ~176 bp; C03/C04/C08/C09 CentBr2-only; A03/A05 CentBr2 |
  | Brassica B | expected `none`/`low`, or a weak pBN35-like call (~496 bp) |
  | Bhybridum | CentBd 156 vs CentBs 157 as separate families; Bd5 low/weak |
  | Ominuta | CentO 155/165 (B) vs CentO-C 126 (C); O. officinalis CEN2 none; CEN7 CentO |
  | Msookensis | Cent728 on all chromosomes |
  | Carabica | no short satellite; long-unit arrays reported |

- Check subgenome vs progenitor consistency (spec §7.3).
- Check the shared diploids across complexes (spec §7.5): identical TRF cache files. Centromere spans agree within 1 kb on the same chromosome between the two complexes.
- Record every disagreement and its explanation in the memo. If a threshold must change, say which and why, rerun only with `--force <stage>`, and record it.

- [ ] **Step 5: Write `clean/memo_centromere.md`** (terse key:value + fenced commands, per CLAUDE.md)

It must contain:
- date, purpose, env and tool versions, inputs and their provenance;
- the exact launch commands and the per-complex `cenrep run` command;
- the outputs layout, a per-complex results table (families, lengths, classes, centromeres by confidence, runtime from `cenrep.time`), and the literature comparison;
- gotchas (expected negatives, shared diploids, parameters changed if any);
- the pytest command and its result;
- `claude --resume "<this session id>"`.

- [ ] **Step 6: Save a project memory**

Write `memory/cenrep-centromere-annotation.md` (type `project`) with:
- where the outputs live;
- the key decisions: tandem-only, satellite-free → `none`, family identity 90%, and the per-chromosome family choice;
- the expected negatives.

Add a one-line pointer to `MEMORY.md`. Link `[[window-identity-allo10-settings]]`.
