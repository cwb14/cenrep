"""Stage orchestration. Checkpoints: <complex>/centromere/work/<stage>.done; a rerun resumes by default."""
import json
import time
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pyhmmer

from . import __version__, annotate, centromere, classify, families, hor, output, plot, te, trf
from .genomes import discover_genomes, load_genes
from .params import MAX_CHROMOSOMES, Params, format_borrow, parse_borrow, to_dict
from .util import CenrepError, log, require_tools

STAGES = ["trf", "families", "annotate", "classify", "centromeres", "hor", "output", "plot"]
# results version per stage: bump a stage when its code can change outputs; a checkpoint's fingerprint carries the
# highest version among it and earlier stages, so a bump invalidates that stage and later ones only.
# 0.3.0 = C1/C2/I1 review fixes; 0.3.2 = classify: strong-family precedence per chromosome;
# 0.3.3 = annotate: shortest-unit rule (drop_multimer_hits); 0.3.4 = annotate: coverage-based overlap resolution;
# 0.3.5 = families: shared monomer frame for related families + near-identical family merge (0.97);
#         centromeres: diffuse flag; plot: stacked bars capped at one bin;
# 0.3.6 = families: trim run-through (terminally redundant) consensus + chain groups; annotate: chain per group.
# 0.3.7 = classify: precedence only for dominant families (largest cluster on >= strong_prevalence of chromosomes).
# 0.4.0 = any number of genomes; annotate: evalue_ref_len; classify: classify_from (borrowed families for contig
#         assemblies) + genes/<genome>.gff3 gene loci; centromeres: borrowed genomes; plot: N panels.
# 0.4.1 = plot: gene density line grey, under the satellite bars.
# 0.4.2 = classify: 45S rDNA spacer repeats (arrays flanked by 18S/5.8S/25S genes) are never centromeric.
# 0.4.3 = inputs: bgzip *.fa.gz genomes (results unchanged: no stage version bump).
# 0.4.4 = families: models founded on a seed whose consensus spans its whole TRF period (base_unit).
# 0.6.0 = classify: dominant column; centromeres: retrotransposon evidence (tes/), honesty cap without it.
STAGE_VERSION = {"trf": "0.3.0", "families": "0.4.4", "annotate": "0.4.0", "classify": "0.6.0",
                 "centromeres": "0.6.0", "hor": "0.3.0", "output": "0.3.0", "plot": "0.4.1"}
# earliest stage each parameter affects; a stage's checkpoint fingerprint covers its own and all earlier stages'
STAGE_PARAMS = {
    "trf": ["trf_chunk", "trf_overlap"],
    "families": ["min_period", "max_period", "min_seed_span", "fam_ident", "multi_ident", "superfam_ident",
                 "absorb_ident", "dedup_ident", "chain_ident", "min_fam_span", "max_families", "hmm_sample", "hmm_sample_bp", "min_hmm_seqs",
                 "nhmmer_evalue"],
    "annotate": ["array_gap_factor", "min_array_monomers", "evalue_ref_len"],
    "classify": ["chrom_min_bp", "cluster_gap", "min_concentration", "max_subtel_frac", "subtel_min_bp",
                 "subtel_frac_len", "max_gene_ratio", "strong_prevalence", "strong_cluster_bp", "classify_from",
                 "rdna_flank", "rdna_arrays", "rdna_array_frac"],
    "centromeres": ["low_sat_bp", "competing_ratio", "fallback_min_bp", "cen_assoc_frac", "min_gap_n",
                    "diffuse_min_len", "diffuse_max_frac", "lu_flank",
                    "lu_window", "lu_window_overlap", "lu_min_unit", "lu_max_unit", "lu_min_ident",
                    "contig_end_bp"],
    "hor": ["hor_ident", "hor_min_cov", "hor_max_len"],
    "output": ["tandem_min_len"],
    "plot": [],
}
TOOLS = ["trf", "makeblastdb", "blastn", "mafft", "cd-hit-est", "minimap2", "cmsearch"]


class Pipeline:
    def __init__(self, complex_dir, threads: int, trf_cache, p: Params, force: str = None, out_dir=None):
        self.cdir = Path(complex_dir).resolve()
        self.name = self.cdir.name
        self.out = Path(out_dir).resolve() if out_dir else self.cdir / "centromere"
        self.work = self.out / "work"
        self.threads = threads
        self.cache = Path(trf_cache).resolve() if trf_cache else self.work / "trf_cache"
        genomes = discover_genomes(self.cdir)
        self.borrow = self._check_borrow(parse_borrow(p.classify_from), {g.name for g in genomes})
        self.p = p = replace(p, classify_from=format_borrow(self.borrow))
        # sources are classified before the genomes that borrow from them
        self.genomes = [g for g in genomes if g.name not in self.borrow] + \
                       [g for g in genomes if g.name in self.borrow]
        for g in self.genomes:
            if g.name not in self.borrow and len(g.chroms) > MAX_CHROMOSOMES:
                log.warning(f"{g.name}: {len(g.chroms)} sequences; classify treats each as a chromosome "
                            f"(prevalence = share of sequences). Use a chromosome-level FASTA, or --classify-from "
                            f"for a contig assembly")
        self.work.mkdir(parents=True, exist_ok=True)
        self.force = force
        stale = next((st for st in STAGES if self._done(st) and self._stored(st) != self.fingerprint(st)), None)
        if stale and (force is None or STAGES.index(force) > STAGES.index(stale)):
            raise CenrepError(f"code version or parameters changed since stage '{stale}' was checkpointed in "
                              f"{self.work}; rerun with --force {stale} (or an earlier stage)")
        if force:
            for st in STAGES[STAGES.index(force):]:
                (self.work / f"{st}.done").unlink(missing_ok=True)
        (self.work / "params.json").write_text(json.dumps(to_dict(p), indent=1))
        self._genes, self._te = {}, {}

    @staticmethod
    def _check_borrow(borrow: dict, names: set) -> dict:
        for tgt, srcs in borrow.items():
            for n in [tgt] + srcs:
                if n not in names:
                    raise CenrepError(f"--classify-from: unknown genome {n!r} (genomes: {', '.join(sorted(names))})")
            if tgt in srcs:
                raise CenrepError(f"--classify-from: {tgt} borrows from itself")
            for s in srcs:
                if s in borrow:
                    raise CenrepError(f"--classify-from: source {s} of {tgt} is itself borrowed")
        return borrow

    def fingerprint(self, stage: str) -> dict:
        cur = to_dict(self.p)
        keys = [k for st in STAGES[:STAGES.index(stage) + 1] for k in STAGE_PARAMS[st]]
        upto = STAGES[:STAGES.index(stage) + 1]
        version = max((STAGE_VERSION[st] for st in upto), key=lambda v: tuple(int(x) for x in v.split(".")))
        return {"version": version, "params": {k: cur[k] for k in keys}}

    def _done(self, stage: str) -> bool:
        return (self.work / f"{stage}.done").exists()

    def _stored(self, stage: str) -> dict:
        try:
            return json.loads((self.work / f"{stage}.done").read_text() or "{}")
        except json.JSONDecodeError:
            return {}

    # ---------------------------------------------------------------- helpers
    def _p(self, stage: str, name: str) -> Path:
        return self.work / stage / name

    def genes(self, g):
        if g.name not in self._genes:
            self._genes[g.name] = load_genes(self.cdir, g.name, g.chroms)
        return self._genes[g.name]

    def trf_table(self, g):
        return trf.load_trf(self._p("trf", f"{g.name}.trf.tsv.gz"))

    def fam_table(self):
        return pd.read_csv(self._p("families", "families.tsv"), sep="\t",
                           dtype={"superfamily": str, "rdna": str, "chain_group": str},
                           keep_default_na=False)

    def fam_info(self) -> dict:
        f = self.fam_table()
        return {r.fam_id: dict(superfamily=r.superfamily, length=int(r.consensus_len), rdna=r.rdna)
                for r in f.itertuples(index=False)}

    def monomers(self, g):
        return pd.read_csv(self._p("annotate", f"{g.name}.monomers.tsv.gz"), sep="\t", dtype={"chrom": str})

    def arrays(self, g):
        return pd.read_csv(self._p("annotate", f"{g.name}.arrays.tsv"), sep="\t", dtype={"chrom": str})

    def fams(self, g, stage: str = "centromeres"):
        name = g if isinstance(g, str) else g.name
        return pd.read_csv(self._p(stage, f"{name}.families.tsv"), sep="\t")

    def regions(self, g):
        return centromere.load_regions(self._p("centromeres", f"{g.name}.regions.tsv"))

    # ---------------------------------------------------------------- driver
    def run(self) -> None:
        require_tools(TOOLS)
        if not families.RRNA_CM.exists():
            raise CenrepError(f"rRNA models missing: {families.RRNA_CM} (part of the cenrep package; reinstall it)")
        t0 = time.time()
        log.info(f"cenrep {__version__}: {self.name} [{', '.join(g.name for g in self.genomes)}], "
                 f"{self.threads} threads")
        rerun = False                                   # once a stage runs, every later stage runs too
        for s in STAGES:
            if not rerun and self._done(s):
                log.info(f"{s}: checkpoint found, skipping")
                continue
            rerun = True
            t = time.time()
            log.info(f"{s}: start")
            (self.work / s).mkdir(exist_ok=True)
            getattr(self, f"stage_{s}")()
            (self.work / f"{s}.done").write_text(json.dumps(self.fingerprint(s)))
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
        groups = dict(zip(fam.fam_id, fam.chain_group)) if "chain_group" in fam else None
        for g in self.genomes:
            annotate.annotate_genome(g, hmms, fam_len, self.threads, self.p,
                                     self._p("annotate", f"{g.name}.monomers.tsv.gz"),
                                     self._p("annotate", f"{g.name}.arrays.tsv"), groups)

    def stage_classify(self):
        for g in self.genomes:
            info = self.fam_info()
            fams, mains = classify.family_metrics(self.arrays(g), g.chroms, self.genes(g), self.p, info)
            if g.name in self.borrow:
                src = {s: self.fams(s, "classify") for s in self.borrow[g.name]}
                fams, per_chrom = classify.borrow_centromeric(fams, mains, src, self.p)
                how = f" (borrowed from {'+'.join(src)})"
            else:
                elig = list(fams.family[fams.eligible.astype(bool)])
                rd = families.rdna_neighbourhood(g.fasta, g.chroms, self.arrays(g), elig,
                                                 self._p("classify", f"rdna_{g.name}"), self.threads,
                                                 self.p.rdna_flank, self.p.rdna_arrays)
                fams = classify.apply_rdna(fams, rd, self.p)
                tagged = [f"{f} ({k}/{n})" for f, (n, k) in rd.items() if n and k / n >= self.p.rdna_array_frac]
                if tagged:
                    log.info(f"classify: {g.name}: rDNA spacer repeats (not centromeric): {', '.join(tagged)}")
                tev = self.te_evidence(g)
                fams, per_chrom, rejected = classify.choose_with_te(fams, mains, g.chroms, self.p, info, tev)
                if rejected:
                    log.info(f"classify: {g.name}: away from the {tev[0]} retrotransposon centromere (not "
                             f"centromeric): {', '.join(f'{f} ({k}/{n})' for f, (k, n) in rejected.items())}")
                how = ""
            fams.to_csv(self._p("classify", f"{g.name}.families.tsv"), sep="\t", index=False)
            prim = fams[fams.cls == "centromeric"]
            log.info(f"classify: {g.name}: centromeric{how} = {', '.join(prim.family) or 'none'}; "
                     f"{len(per_chrom)}/{len(g.chroms)} sequences have one")

    def te_evidence(self, g):
        """(clade, {chrom: significant core}, decisive) of the centromere-marking LTR clade (te.marker), or None
        (cached)."""
        if g.name not in self._te:
            self._te[g.name] = self._te_evidence(g)
        return self._te[g.name]

    def _te_evidence(self, g):
        gff = te.te_gff(self.cdir, g.name)
        if gff is None:
            log.info(f"{g.name}: no tes/{g.name}.gff3[.gz]; satellite evidence only")
            return None
        cands = te.candidate_clades(te.load_te_gff(gff, g.chroms), g.chroms)
        if not cands:
            log.info(f"{g.name}: no LTR clade concentrated at one locus on most chromosomes; "
                     f"satellite evidence only")
            return None
        m = te.marker(cands)
        log.info(f"{g.name}: LTR clades concentrated at one locus per sequence: "
                 f"{', '.join(f'{c} ({len(x)}/{len(g.chroms)})' for c, x in cands)}; marker {m[0]}"
                 f"{' (decisive)' if m[2] else ' (flags only)'}")
        return m

    def stage_centromeres(self):
        for g in self.genomes:
            arrays = self.arrays(g)
            fams = self.fams(g, "classify")
            primary = set(fams.family[fams.cls == "centromeric"])
            borrowed = g.name in self.borrow
            tev = None if borrowed else self.te_evidence(g)
            dominant = fams["dominant"].astype(str).eq("True").any() if "dominant" in fams else True
            capped = tev is None and not borrowed and bool(primary) and not dominant

            def call_with(prim):
                return centromere.call_regions(g, arrays, prim, self.genes(g), self.trf_table(g),
                                               self.work / "centromeres" / f"tmp_{g.name}", self.threads, self.p,
                                               borrowed=borrowed, te=tev, capped=capped)
            rows, lu = call_with(primary)
            if capped:
                log.info(f"centromeres: {g.name}: no dominant satellite and no retrotransposon evidence; "
                         f"calls capped at low confidence")
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
                    hor_units=rd(f"{g.name}.hor_units.tsv"), trf=self.trf_table(g), genes=self.genes(g),
                    fam_info=self.fam_info())

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
                               genes=self.genes(g), borrowed=g.name in self.borrow))
        plot.plot_overview(self.out / "centromere_overview.pdf", panels)
