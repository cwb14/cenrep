"""Command line: cenrep run --complex-dir DIR --threads N --trf-cache DIR [--classify-from T=S] [--force STAGE] [-v]
[thresholds]"""
import argparse
from dataclasses import fields
from pathlib import Path

from . import __version__
from .params import HELP, Params
from .util import CenrepError, log, setup_logging


def main(argv=None) -> int:
    from .pipeline import STAGES, Pipeline
    ap = argparse.ArgumentParser(prog="cenrep", description=(
        "Centromeric satellite, HOR and tandem-repeat annotation for one complex dir of genomes (one *.fa or bgzip "
        "*.fa.gz each, "
        "chromosome-level; contig assemblies via --classify-from). Satellite families are shared across all "
        "genomes. Outputs go to <complex-dir>/centromere/."))
    ap.add_argument("--version", action="version", version=f"cenrep {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run all stages; resumes from checkpoints in centromere/work/")
    r.add_argument("--complex-dir", required=True, type=Path,
                   help="dir with one *.fa (+ .fai) or bgzip *.fa.gz (+ .fai, .gzi) per genome; optional "
                        "genes/<genome>.gff3[.gz]")
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
