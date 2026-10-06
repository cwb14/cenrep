"""Command line: cenrep run GENOME_DIR [-o OUTDIR] [-t THREADS] [--trf-cache DIR] [--classify-from T=S] [--force STAGE]
[-v]. Every threshold is an internal constant (params.py); none is an option."""
import argparse
from pathlib import Path

from . import __version__
from .params import Params
from .util import CenrepError, log, setup_logging


def main(argv=None) -> int:
    from .pipeline import STAGES, Pipeline
    ap = argparse.ArgumentParser(prog="cenrep", description=(
        "Centromeric satellite, HOR and tandem-repeat annotation. Finds satellite families de novo, shared across "
        "all genomes in GENOME_DIR, and calls one centromere per chromosome with its evidence and confidence."))
    ap.add_argument("--version", action="version", version=f"cenrep {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run all stages; resumes from checkpoints in OUTDIR/work/",
                       description="Run all stages; a rerun resumes from checkpoints in OUTDIR/work/.")
    r.add_argument("genome_dir", nargs="?", type=Path, metavar="GENOME_DIR",
                   help="dir with one FASTA per genome (*.fa/*.fasta/*.fna, plain or bgzip, indexed with "
                        "samtools faidx); optional GENOME_DIR/genes/<genome>.gff3[.gz] (gene loci) and "
                        "GENOME_DIR/tes/<genome>.gff3[.gz] (LTR retrotransposons with clade=)")
    r.add_argument("--complex-dir", dest="complex_dir", type=Path, help=argparse.SUPPRESS)   # alias of GENOME_DIR
    r.add_argument("-o", "--outdir", type=Path, help="output dir [GENOME_DIR/centromere]")
    r.add_argument("-t", "--threads", type=int, default=8, help="threads [8]")
    r.add_argument("--trf-cache", type=Path, metavar="DIR", help="TRF result cache, shareable between runs [OUTDIR/work/trf_cache]")
    r.add_argument("--classify-from", default="", metavar="T=S[+S2],...",
                   help="genome T (e.g. a contig assembly) takes the centromeric families of genome(s) S")
    r.add_argument("--force", choices=STAGES, help="recompute this stage and all later ones")
    r.add_argument("-v", "--verbose", action="store_true", help="per-step progress and sanity checks")
    a = ap.parse_args(argv)
    gdir = a.genome_dir or a.complex_dir
    if gdir is None or (a.genome_dir and a.complex_dir and a.genome_dir != a.complex_dir):
        ap.error("give one GENOME_DIR")
    setup_logging(a.verbose)
    try:
        Pipeline(gdir, a.threads, a.trf_cache, Params(classify_from=a.classify_from), force=a.force,
                 out_dir=a.outdir).run()
    except CenrepError as e:
        log.error(str(e))
        return 1
    return 0
