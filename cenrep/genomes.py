"""Complex-directory inputs: genomes, chromosome lengths, gene positions."""
import gzip
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
    """One genome per *.fa, or bgzip *.fa.gz with .fai + .gzi (read through pysam); name = file name without
    .fa[.gz]. Genomes are ordered by file name."""
    complex_dir = Path(complex_dir)
    fas = sorted(list(complex_dir.glob("*.fa")) + list(complex_dir.glob("*.fa.gz")), key=lambda f: f.name)
    if not fas:
        raise CenrepError(f"{complex_dir}: no *.fa genomes (plain, or bgzip *.fa.gz)")
    seen = {}
    for f in fas:
        name = f.name[:-len(".fa.gz")] if f.name.endswith(".fa.gz") else f.stem
        if name in seen:
            raise CenrepError(f"{complex_dir}: genome {name} is given twice ({seen[name].name}, {f.name})")
        if f.name.endswith(".gz") and not Path(str(f) + ".gzi").exists():
            raise CenrepError(f"{f}: a compressed genome must be bgzip with .fai + .gzi indexes "
                              f"(missing {f.name}.gzi; run: samtools faidx {f})")
        seen[name] = f
    return [Genome(n, f, load_fai(f)) for n, f in seen.items()]


def find_gene_gff(complex_dir, genome: str):
    for gff in sorted(Path(complex_dir).glob(f"window_identity/*/{genome}_mod.gff")):
        ex = gff.parent / "run.exit"
        if ex.exists() and ex.read_text().strip() == "0":
            return gff
    return None


def native_gene_gff(complex_dir, genome: str):
    for ext in (".gff3", ".gff3.gz"):
        f = Path(complex_dir) / "genes" / f"{genome}{ext}"
        if f.exists():
            return f
    return None


def load_native_genes(gff, chroms=None) -> tuple:
    """One midpoint per gene locus (isoforms counted once); loci with a locus_type other than coding skipped.
    Returns ({chrom: sorted midpoints}, number of gene lines on sequences not in chroms)."""
    mids, foreign = {}, 0
    with (gzip.open(gff, "rt") if str(gff).endswith(".gz") else open(gff)) as f:
        for n, line in enumerate(f, 1):
            if line.startswith("#"):
                continue
            p = line.rstrip("\r\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            if chroms is not None and p[0] not in chroms:
                foreign += 1
                continue
            attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            if attrs.get("locus_type", "coding") != "coding":
                continue
            try:
                mids.setdefault(p[0], []).append((int(p[3]) - 1 + int(p[4])) // 2)
            except ValueError:
                raise CenrepError(f"{gff}:{n}: bad gene coordinates {p[3]!r}..{p[4]!r}") from None
    return {c: np.sort(np.array(v, dtype=np.int64)) for c, v in mids.items()}, foreign


def load_genes(complex_dir, genome: str, chroms=None):
    """{chrom: sorted np.int64 array of gene midpoints (0-based)}, or None without gene evidence.
    Source: <complex>/genes/<genome>.gff3[.gz] (gene lines, native seqids; only chroms when given), else the
    miniprot mRNA of a finished window_identity run."""
    complex_dir = Path(complex_dir)
    native = native_gene_gff(complex_dir, genome)
    if native is not None:
        mids, foreign = load_native_genes(native, chroms)
        if not mids:
            log.warning(f"{genome}: no gene in {native} lies on a sequence of the FASTA "
                        f"({foreign} on other sequence names); gene evidence off (no_genes)")
            return None
        if foreign:
            log.info(f"{genome}: {foreign} genes in {native.name} on sequences not in the FASTA ignored")
        return mids
    gff = find_gene_gff(complex_dir, genome)
    cmap = complex_dir / "window_identity" / "genomes" / "chr_map.tsv"
    if gff is None or not cmap.exists():
        log.warning(f"{genome}: no genes/{genome}.gff3 or finished window_identity gene GFF; "
                    f"gene evidence off (no_genes)")
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
