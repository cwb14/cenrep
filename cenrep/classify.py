"""S4: per-family metrics and the per-chromosome centromeric-family choice."""
import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .genomes import count_in
from .te import agrees
from .params import RELATED_LEN_TOL
from .util import clipped_bp, log, union_bp

CLS_COLS = ["family", "bp", "n_arrays", "prevalence", "concentration", "subtel_frac", "gene_ratio",
            "median_main_bp", "eligible", "cls", "strength", "defines", "dominant", "reason"]


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
    """The terminal subtel_frac_len of a chromosome at each end. No absolute floor: a 1-Mb floor made every chromosome
    <= 2 Mb all end zone (Cryptococcus) and covered a CENH3 centromere 1.0 Mb from a chromosome end (Fragaria x
    ananassa Chr1D, Type II satellite; Jin 2025)."""
    return int(p.subtel_frac_len * L)


def terminal_bp(iv, L: int, p, mids=None) -> int:
    """bp of iv in the subtelomeres: the stretch between each chromosome end and its outermost gene (the textbook
    subtelomere). A centromere has a gene-bearing arm on both sides even next to an end (Fragaria x ananassa Chr1D:
    CENH3-bound Type II array 1.0 Mb from the end, 5 genes beyond it); subtelomeric repeats have none beyond them,
    whether at one end or both (rice TrsC, Solanum 182-bp blocks). Without >= 2 gene loci on the sequence: the
    terminal subtel_frac_len at each end."""
    if mids is not None and len(mids) >= 2 and mids[0] < mids[-1]:
        lo, hi = int(mids[0]), int(mids[-1])
    else:
        z = end_zone(L, p)
        lo, hi = z, L - z
    return clipped_bp(iv, 0, lo) + clipped_bp(iv, hi, L)


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
    """fam_info: {fam: {superfamily, length, rdna}} from the families stage (optional). A family with most of its bp in
    subtelomeres (terminal_bp) is subtelomeric."""
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
            subtel += terminal_bp(iv, chroms[c], p, genes.get(c) if genes else None)
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
                         eligible=not fails, cls="other", strength=".", defines=".", dominant=False,
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
    fams["dominant"] = fams.family.isin(dominant)
    for i, r in fams.iterrows():
        if r.family in defines:
            strong = r.prevalence >= p.strong_prevalence and r.median_main_bp >= p.strong_cluster_bp
            fams.loc[i, ["cls", "strength", "defines", "reason"]] = [
                "centromeric", "strong" if strong else "weak", ",".join(defines[r.family]),
                f"eligible; centromeric family on {len(defines[r.family])} chromosome(s)"]
        elif r.eligible:
            fams.loc[i, "reason"] = "eligible; never the main or co-located cluster on a chromosome"
    return fams, per_chrom


def apply_rdna(fams: pd.DataFrame, rdna: dict, p, elsewhere: dict = None) -> pd.DataFrame:
    """Families whose largest arrays sit among 45S rRNA genes (rdna_neighbourhood) are rDNA spacer repeats (NORs):
    never centromeric. elsewhere: {fam: (genome, k, n)} tagged so in another genome of the complex. rRNA genes are
    often missing from assemblies, so their presence next to the arrays in one genome is the evidence, and their
    absence in another is not (Lemna: one family flanked by 45S genes at 5/9 and 6/10 arrays in two genomes, 1/6 in a
    third, where it was called centromeric)."""
    fams = fams.copy()
    for i, r in fams.iterrows():
        n, k = rdna.get(r.family, (0, 0))
        if n and k / n >= p.rdna_array_frac:
            fams.loc[i, ["eligible", "reason"]] = [False, f"rDNA-associated: 45S rRNA gene within {p.rdna_flank} bp "
                                                          f"of {k}/{n} largest arrays (rDNA spacer repeat)"]
        elif elsewhere and r.family in elsewhere and r.eligible:
            g, k2, n2 = elsewhere[r.family]
            fams.loc[i, ["eligible", "reason"]] = [False, f"rDNA-associated in {g}: 45S rRNA gene within "
                                                          f"{p.rdna_flank} bp of {k2}/{n2} largest arrays there"]
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


def te_rejected(fams: pd.DataFrame, mains: dict, cores, fam_info=None) -> dict:
    """Centromeric family kinds (related families: subfamilies, CentBr1/2, CEN180/AthCEN159) whose main clusters
    disagree with the retrotransposon core on the majority of the chromosomes they define, for EVERY concentrated
    clade (cores: one {chrom: core} per candidate clade, or a single dict). A chromosome agrees if any member of the
    kind agrees there; judged per kind, so a subfamily on a few chromosomes is not rejected alone; judged against
    every clade, because which concentrated clade is centromeric is not known. Returns {family: (disagreeing,
    tested)} (counts for the most supportive clade) for every member of each rejected kind."""
    cands = [cores] if isinstance(cores, dict) else list(cores)
    cen = [r for r in fams[fams.cls == "centromeric"].itertuples(index=False)]
    kinds = []                                           # connected components of `related`
    for r in cen:
        hit = [k for k in kinds if any(related(r.family, x.family, fam_info) for x in k)]
        merged = [r] + [x for k in hit for x in k]
        kinds = [k for k in kinds if k not in hit] + [merged]
    out = {}
    for k in kinds:
        tallies = []
        for cc in cands:
            votes = {}
            for r in k:
                for c in str(r.defines).split(","):
                    m, core = mains.get((r.family, c)), cc.get(c)
                    if core and m:
                        votes[c] = votes.get(c, False) or agrees(m.start, m.end, core)
            if votes:
                tallies.append((sum(not v for v in votes.values()), len(votes)))
        if tallies and all(bad * 2 > n for bad, n in tallies):
            best = min(tallies, key=lambda t: t[0] / t[1])
            out.update({r.family: best for r in k})
    return out


def choose_with_te(fams: pd.DataFrame, mains: dict, chroms: dict, p, fam_info, te):
    """choose_centromeric, then drop families away from the retrotransposon centromere and choose again, until no
    centromeric family disagrees: a family that was outcompeted by a rejected one (e.g. the CentT550-like satellite
    beside interstitial blocks in Aegilops) gets its chance. te: (clade, {chrom: significant core}, decisive) or
    None; only decisive evidence (CRM, te.marker) rejects. Returns (fams, per_chrom, {rejected family: (disagreeing, tested)})."""
    rejected = {}
    while True:
        out, per_chrom = choose_centromeric(fams, mains, chroms, p, fam_info)
        rej = te_rejected(out, mains, [te[1]], fam_info) if te and te[2] else {}
        if not rej:
            return out, per_chrom, rejected
        rejected.update(rej)
        fams = fams.copy()
        for i, r in fams.iterrows():
            if r.family in rej:
                k, n = rej[r.family]
                fams.loc[i, ["eligible", "reason"]] = [
                    False, f"main clusters away from the {te[0]} retrotransposon centromere on {k}/{n} chromosomes"]
