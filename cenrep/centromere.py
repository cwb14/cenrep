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


def finalize(chrom: str, L: int, call, fa, genes_c, chrom_density, p, borrowed: bool = False) -> dict:
    """borrowed (contig assemblies): near_end is replaced by contig_end (call within contig_end_bp of a contig
    end: the array may continue in another contig), which does not lower confidence."""
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
    if not borrowed and (s < z or e > L - z):
        flags.append("near_end")
    if e - s >= p.diffuse_min_len and call["sat_bp"] < p.diffuse_max_frac * (e - s):
        flags.append("diffuse")                        # broad, satellite-sparse block (e.g. Coffea heterochromatin)
    core = len(flags)
    if borrowed and (s < p.contig_end_bp or e > L - p.contig_end_bp):
        flags.append("contig_end")
    n_genes, gene_ratio = ".", "."
    if genes_c is not None:
        n_genes = count_in(genes_c, s, e)
        if chrom_density:
            gene_ratio = round((n_genes / (e - s)) / chrom_density, 3)
            if gene_ratio > 1:
                flags.append("gene_rich")
    if runs:
        flags.append(f"gaps:{len(runs)}")
    conf = ["high", "medium", "low", "low", "low"][core] if call["evidence"] == "satellite" else "low"
    row.update(call)
    row.update(confidence=conf, flags=",".join(flags) or ".", n_gaps=len(runs), gap_bp=sum(runs),
               n_genes=n_genes, gene_ratio=gene_ratio)
    return row


def call_regions(genome, arrays, primary: set, genes, trf_rows, tmpdir, threads: int, p, borrowed: bool = False):
    """borrowed (contig assemblies with families from --classify-from): only full satellite calls count. A contig
    without the centromeric satellite is no evidence of a satellite-free centromere, so there is no weak-satellite
    or long-unit fallback, and the long-unit scan covers call flanks only."""
    tmpdir.mkdir(parents=True, exist_ok=True)
    total = sum(genome.chroms.values())
    gd = (sum(len(v) for v in genes.values()) / total) if genes else None
    by_chrom = dict(tuple(arrays.groupby("chrom"))) if len(arrays) else {}
    sat = {c: satellite_call(by_chrom[c], primary, L, p) if c in by_chrom else None for c, L in genome.chroms.items()}
    if borrowed:
        sat = {c: x if x and x["evidence"] == "satellite" else None for c, x in sat.items()}
    targets = [(c, max(0, sat[c]["cen_start"] - p.lu_flank), min(L, sat[c]["cen_end"] + p.lu_flank))
               if sat[c] else (c, 0, L) for c, L in genome.chroms.items() if sat[c] or not borrowed]
    lu = longunit.detect(genome.fasta, targets, longunit.build_masks([arrays, trf_rows]), tmpdir, threads, p)
    rows = {}
    lu_by = dict(tuple(lu.groupby("chrom"))) if len(lu) else {}
    with pysam.FastaFile(str(genome.fasta)) as fa:
        for c, L in genome.chroms.items():
            gc = None if genes is None else genes.get(c, EMPTY)
            call = sat[c]
            if call is None and not borrowed:
                call = longunit_call(lu_by.get(c, lu.iloc[:0]), L, gc, gd, p)
            rows[c] = finalize(c, L, call, fa, gc, (len(gc) / L) if gc is not None else None, p, borrowed)
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
