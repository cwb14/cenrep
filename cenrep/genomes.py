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


FASTA_EXT = (".fa", ".fasta", ".fna")


def _fasta_name(f: Path):
    """Genome name = file name without its FASTA extension (and .gz), or None if f is not a FASTA."""
    n = f.name[:-3] if f.name.endswith(".gz") else f.name
    ext = next((e for e in FASTA_EXT if n.endswith(e)), None)
    return n[:-len(ext)] if ext else None


def discover_genomes(complex_dir) -> list:
    """One genome per FASTA (*.fa, *.fasta, *.fna; plain with .fai, or bgzip .gz with .fai + .gzi; read through
    pysam); name = file name without the extension. Genomes are ordered by file name."""
    complex_dir = Path(complex_dir)
    fas = sorted((f for f in complex_dir.iterdir() if f.is_file() and _fasta_name(f)), key=lambda f: f.name) \
        if complex_dir.is_dir() else []
    if not fas:
        raise CenrepError(f"{complex_dir}: no *.fa genomes (*.fa/*.fasta/*.fna, plain or bgzip .gz)")
    seen = {}
    for f in fas:
        name = _fasta_name(f)
        if name in seen:
            raise CenrepError(f"{complex_dir}: genome {name} is given twice ({seen[name].name}, {f.name})")
        if f.name.endswith(".gz") and not Path(str(f) + ".gzi").exists():
            raise CenrepError(f"{f}: a compressed genome must be bgzip with .fai + .gzi indexes "
                              f"(missing {f.name}.gzi; run: samtools faidx {f})")
        seen[name] = f
    return [Genome(n, f, load_fai(f)) for n, f in seen.items()]


def gene_gff(complex_dir, genome: str):
    for ext in (".gff3", ".gff3.gz"):
        f = Path(complex_dir) / "genes" / f"{genome}{ext}"
        if f.exists():
            return f
    return None


def load_gene_gff(gff, chroms=None) -> tuple:
    """One midpoint per gene locus from `gene` lines (isoforms counted once). Non-coding loci are skipped: a
    locus_type other than coding, or a gene_biotype/biotype other than protein_coding (NCBI/Ensembl `gene` lines
    include tRNA, ncRNA and pseudogene loci). Returns ({chrom: sorted midpoints}, gene lines on sequences not in
    chroms)."""
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
            if attrs.get("gene_biotype", attrs.get("biotype", "protein_coding")) != "protein_coding":
                continue
            try:
                mids.setdefault(p[0], []).append((int(p[3]) - 1 + int(p[4])) // 2)
            except ValueError:
                raise CenrepError(f"{gff}:{n}: bad gene coordinates {p[3]!r}..{p[4]!r}") from None
    return {c: np.sort(np.array(v, dtype=np.int64)) for c, v in mids.items()}, foreign


def load_genes(complex_dir, genome: str, chroms=None):
    """{chrom: sorted np.int64 array of gene midpoints (0-based)}, or None without gene evidence.
    Source: <complex>/genes/<genome>.gff3[.gz] (gene lines, native seqids; only chroms when given)."""
    gff = gene_gff(complex_dir, genome)
    if gff is None:
        log.warning(f"{genome}: no genes/{genome}.gff3[.gz]; gene evidence off (no_genes)")
        return None
    mids, foreign = load_gene_gff(gff, chroms)
    if not mids:
        log.warning(f"{genome}: no gene in {gff} lies on a sequence of the FASTA "
                    f"({foreign} on other sequence names); gene evidence off (no_genes)")
        return None
    if foreign:
        log.info(f"{genome}: {foreign} genes in {gff.name} on sequences not in the FASTA ignored")
    return mids


def count_in(mids, start: int, end: int) -> int:
    if mids is None or len(mids) == 0:
        return 0
    return int(np.searchsorted(mids, end) - np.searchsorted(mids, start))
