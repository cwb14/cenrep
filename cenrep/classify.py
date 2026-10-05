"""S4: per-family metrics and the per-chromosome centromeric-family choice."""
import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .genomes import count_in
from .params import RELATED_LEN_TOL
from .util import clipped_bp, log, union_bp

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


def related(f: str, g: str, fam_info) -> bool:
    """Families that may share a centromere as layers: same superfamily or period within 15%
    (CentBr1/2, CentO 155/165, CEN180/AthCEN159). Without family info every pair counts as related."""
    if fam_info is None:
        return True
    a, b = fam_info[f], fam_info[g]
    return a["superfamily"] == b["superfamily"] or abs(a["length"] / b["length"] - 1) <= RELATED_LEN_TOL


def family_metrics(arrays: pd.DataFrame, chroms: dict, genes, p, fam_info=None):
    """fam_info: {fam: {superfamily, length, rdna}} from the families stage (optional)."""
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
        rdna = fam_info[fam]["rdna"] if fam_info and fam in fam_info else "."
        if rdna != ".":
            fails.append(f"rDNA ({rdna})")
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


def choose_centromeric(fams: pd.DataFrame, mains: dict, chroms: dict, p, fam_info=None):
    elig = set(fams.loc[fams.eligible.astype(bool), "family"])
    by_chrom = {c: sorted(((mains[(f, c)].bp, f) for f in elig
                           if (f, c) in mains and mains[(f, c)].bp >= p.chrom_min_bp), reverse=True) for c in chroms}
    # genome-wide dominant family: strong AND its kind is the largest cluster on >= strong_prevalence of
    # chromosomes; presence alone let a widespread pericentromeric satellite (B. napus A06, 228 bp) take precedence
    largest = [cands[0][1] for cands in by_chrom.values() if cands]
    dominant = {r.family for r in fams.itertuples(index=False) if r.family in elig
                and r.prevalence >= p.strong_prevalence and r.median_main_bp >= p.strong_cluster_bp
                and sum(related(w, r.family, fam_info) for w in largest) >= p.strong_prevalence * len(chroms)}
    per_chrom = {}
    for c, cands in by_chrom.items():
        if not cands:
            continue
        # a family related to a genome-wide dominant family wins where it has a qualifying cluster;
        # otherwise a larger pericentromeric satellite (e.g. TR238 next to CentBr) would define the centromere
        pref = [bf for bf in cands if any(related(bf[1], s, fam_info) for s in dominant)]
        best = (pref or cands)[0][1]
        bc = mains[(best, c)]
        per_chrom[c] = [best] + [f for _, f in cands if f != best and interval_gap(mains[(f, c)], bc) <= p.cluster_gap
                                 and related(f, best, fam_info)]
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


def apply_rdna(fams: pd.DataFrame, rdna: dict, p) -> pd.DataFrame:
    """Families whose largest arrays sit among 45S rRNA genes (rdna_neighbourhood) are rDNA spacer repeats (NORs):
    never centromeric."""
    fams = fams.copy()
    for i, r in fams.iterrows():
        n, k = rdna.get(r.family, (0, 0))
        if n and k / n >= p.rdna_array_frac:
            fams.loc[i, ["eligible", "reason"]] = [False, f"rDNA-associated: 45S rRNA gene within {p.rdna_flank} bp "
                                                          f"of {k}/{n} largest arrays (rDNA spacer repeat)"]
    return fams


def borrow_centromeric(fams: pd.DataFrame, mains: dict, sources: dict, p):
    """Contig assemblies: the per-chromosome tests (prevalence, end zones, concentration) are meaningless on
    contigs, so the target takes the union of its source genomes' centromeric families (strongest strength);
    its own metrics (eligible, reason) stay in the table for reference only. sources: {genome: classify table}.
    defines = target sequences where the family has a main cluster >= chrom_min_bp."""
    rank = {"strong": 2, "weak": 1}
    cen = {}
    for t in sources.values():
        for r in t[t.cls == "centromeric"].itertuples(index=False):
            if rank.get(r.strength, 0) >= rank.get(cen.get(r.family), 0):
                cen[r.family] = r.strength
    src = "+".join(sources)
    missing = sorted(set(cen) - set(fams.family))
    if missing:
        log.info(f"classify: centromeric in {src} but without arrays in the target: {', '.join(missing)}")
    per_chrom = {}
    fams = fams.copy()
    for i, r in fams.iterrows():
        if r.family in cen:
            chroms = [c for (f, c), m in mains.items() if f == r.family and m.bp >= p.chrom_min_bp]
            for c in chroms:
                per_chrom.setdefault(c, []).append(r.family)
            fams.loc[i, ["cls", "strength", "defines", "reason"]] = [
                "centromeric", cen[r.family], ",".join(chroms) or ".",
                f"borrowed from {src} ({cen[r.family]} there); own contig metrics: {r.reason}"]
        else:
            fams.loc[i, "reason"] = f"not centromeric in {src}; own contig metrics: {r.reason}"
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
