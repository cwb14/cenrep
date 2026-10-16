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
                "gap_bp", "n_genes", "gene_ratio", "other_tandem_bp", "te_clade", "te_core_start", "te_core_end"]
FAMTAB_COLS = ["family", "fam_id", "superfamily", "consensus_len", "hmm_mode", "rdna", "genome", "class", "role", "strength",
               "defines", "bp", "n_arrays", "prevalence", "concentration", "subtel_frac", "gene_ratio",
               "median_main_bp", "eligible", "dominant", "in_ltr", "ltr_clade", "reason"]
LEVEL = {"centromere": 0, "satellite_DNA": 1, "tandem_repeat": 1, "higher_order_repeat": 2, "satellite_monomer": 3}
_ESC = {"%": "%25", ";": "%3B", "=": "%3D", "&": "%26", "\t": "%09", "\n": "%0A", "\r": "%0D"}


def _esc(v) -> str:
    return "".join(_ESC.get(ch, ch) for ch in str(v))


def attrs(d: dict) -> str:
    """GFF3 column 9; commas are kept (multi-value lists), reserved characters percent-encoded."""
    return ";".join(f"{k}={_esc(v)}" for k, v in d.items() if v is not None)


def tandem_loci(rows: pd.DataFrame) -> list:
    """Merge overlapping seed-table rows; a locus takes period/copies/pct_match/consensus of its top-scoring row."""
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
    by = lambda df: {c: x for c, x in df.groupby("chrom", sort=False)} if len(df) else {}
    inc_by, mon_by, units_by, tr_by, lu_by = (by(x) for x in (inc, mon, d["hor_units"], d["seeds"], d["longunit"]))
    for c, rg in reg.items():                            # per-sequence tables: ~3,000 contigs stay cheap
        s0, e0 = rg["cen_start"], rg["cen_end"]
        ic = inc_by.get(c, inc.iloc[:0])
        a = ic[ic.inside & ic.family.isin(primary)]
        mc = mon_by.get(c, mon.iloc[:0])
        m = mc[mc.array_id.isin(a.array_id)]
        classed, inhor = int(m.classed.sum()), int(m.inhor.sum())
        units = units_by.get(c, d["hor_units"].iloc[:0])
        top = "."
        if len(units) and units.max_tandem_copies.max() > 1:
            u = units.sort_values(["max_tandem_copies", "n_monomers"], ascending=False).iloc[0]
            top = f"{u.unit}:{u.n_monomers}x{u.max_tandem_copies}"
        stats[c] = dict(cen_sat_bp=union_bp(zip(a.start, a.end)), n_arrays=len(a),
                        largest_array_bp=int((a.end - a.start).max()) if len(a) else 0, n_monomers=len(m),
                        mean_identity=round(float(m.identity.mean()), 3) if len(m) else ".",
                        hor_monomers=inhor, hor_fraction=round(inhor / classed, 3) if classed else ".",
                        top_tandem_hor=top)
        sat_iv = list(zip(ic.start, ic.end))
        cand = []
        t = tr_by.get(c, d["seeds"].iloc[:0])
        t = t[(t.start < e0) & (t.end > s0)]
        # drop satellite-derived seed rows BEFORE merging, so a satellite row that bleeds a few bp into an
        # adjacent microsatellite cannot swallow it into one locus
        t = t.loc[[clipped_bp(sat_iv, s, e) < 0.5 * (e - s) for s, e in zip(t.start, t.end)]]
        for L in tandem_loci(t):
            cand.append((L["start"], L["end"], dict(period=L["period"], copies=L["copies"],
                                                    pct_match=L["pct_match"], unit_source=p.seed_engine,
                                                    consensus=L["consensus"])))
        for r in lu_by.get(c, d["longunit"].iloc[:0]).itertuples(index=False):
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
            n_arrays=st["n_arrays"], hor_fraction=st["hor_fraction"], n_genes=rg["n_genes"], gaps=rg["n_gaps"],
            te_clade=rg.get("te_clade", "."),
            te_core=f"{int(rg['te_core_start']) + 1}-{rg['te_core_end']}" if rg.get("te_core_start", ".") != "." else ".")))
    for r in inc.itertuples(index=False):
        hf = "."
        if r.array_id in per_arr.index and per_arr.loc[r.array_id, "classed"]:
            hf = round(per_arr.loc[r.array_id, "inhor"] / per_arr.loc[r.array_id, "classed"], 3)
        a_ = dict(ID=r.gid)
        if r.inside:
            a_["Parent"] = f"cen_{r.chrom}"
        rdna = (d.get("fam_info") or {}).get(r.family, {}).get("rdna", ".")
        if rdna != ".":
            a_["rdna"] = rdna
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

    write_arrays(g, arrays, fams, d.get("fam_info"), nm, order, outdir)

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
                         n_genes=rg["n_genes"], gene_ratio=rg["gene_ratio"], other_tandem_bp=union_bp(other_iv),
                         te_clade=rg.get("te_clade", "."),
                         te_core_start=int(rg["te_core_start"]) + 1 if rg.get("te_core_start", ".") != "." else ".",
                         te_core_end=rg.get("te_core_end", ".")))
    pd.DataFrame(rows, columns=SUMMARY_COLS).to_csv(outdir / f"{g.name}_centromere_summary.tsv", sep="\t", index=False)
    hu = d["hor_units"].copy()
    hu["family"] = hu.family.map(nm)
    hu.to_csv(outdir / f"{g.name}_hor_units.tsv", sep="\t", index=False)


def write_arrays(g, arrays: pd.DataFrame, fams: pd.DataFrame, fam_info, nm, order: dict, outdir) -> None:
    """<genome>_satellite_arrays.gff: every array of every family cenrep keeps, on every sequence, with the family's
    role in this genome (classify.assign_roles) and its rDNA tag. arrays: with column ov (bp inside the centromere
    call)."""
    role = dict(zip(fams.family, fams.role)) if "role" in fams else {}
    tag = dict(zip(fams.family, fams.rdna_tag)) if "rdna_tag" in fams else {}
    a = arrays.assign(ci=arrays.chrom.map(order)).sort_values(["ci", "start", "end"])
    k = a.groupby("chrom").cumcount() + 1
    with open(Path(outdir) / f"{g.name}_satellite_arrays.gff", "w") as f:
        f.write("##gff-version 3\n")
        for c, L in g.chroms.items():
            f.write(f"##sequence-region {c} 1 {L}\n")
        f.write(f"# cenrep {__version__} {date.today()} genome={g.name}: every array of every satellite family; "
                f"role = what the family is in this genome (README, Outputs)\n")
        for r, n in zip(a.itertuples(index=False), k):
            rd = tag.get(r.family, ".")
            if not isinstance(rd, str) or rd in ("", "."):
                rd = (fam_info or {}).get(r.family, {}).get("rdna", ".")
            a_ = dict(ID=f"sat_{r.chrom}_{n}", family=nm(r.family), role=role.get(r.family, "."),
                      rdna=rd if rd != "." else None, monomer_len=r.monomer_len, n_monomers=r.n_monomers,
                      mean_identity=round(r.mean_identity, 3),
                      in_centromere="true" if r.ov >= 0.5 * (r.end - r.start) else "false")
            f.write(f"{r.chrom}\tcenrep\tsatellite_DNA\t{r.start + 1}\t{r.end}\t.\t{r.strand}\t.\t{attrs(a_)}\n")


def write_families(fam: pd.DataFrame, data: dict, names: dict, hmm_path, fa_path, outdir) -> None:
    outdir = Path(outdir)
    rows = []
    for gname, d in data.items():
        f = d["fams"].set_index("family")
        for r in fam.itertuples(index=False):
            row = dict(family=names[r.fam_id], fam_id=r.fam_id, superfamily=names.get(r.superfamily, r.superfamily),
                       consensus_len=r.consensus_len, hmm_mode=r.hmm_mode, rdna=getattr(r, "rdna", "."), genome=gname)
            if r.fam_id in f.index:
                x = f.loc[r.fam_id]
                row.update({"class": x.cls, "role": getattr(x, "role", "."), "strength": x.strength,
                            "defines": x.defines, "bp": x.bp,
                            "n_arrays": x.n_arrays, "prevalence": x.prevalence, "concentration": x.concentration,
                            "subtel_frac": x.subtel_frac, "gene_ratio": x.gene_ratio,
                            "median_main_bp": x.median_main_bp, "eligible": x.eligible,
                            "dominant": getattr(x, "dominant", False), "in_ltr": getattr(x, "in_ltr", float("nan")),
                            "ltr_clade": getattr(x, "ltr_clade", "."), "reason": x.reason})
            else:
                row.update({"class": "absent", "role": "absent", "bp": 0, "n_arrays": 0,
                            "reason": "no arrays in this genome"})
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
