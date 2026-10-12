"""cenrep's internal constants: every analysis threshold in one place, each with a one-line reason.

None of them is a command-line option: a user should not have to tune 51 numbers, and a value tuned to one dataset
does not travel. `Params` holds the stage-fingerprinted values; checkpoints (work/<stage>.done) are keyed by its
field names, so a field is renamed only together with a stage-version bump. Values that express one concept share
one named constant. 'prior' marks an absolute size with no data-derived rule yet. Plot layout and lock polling
constants live with their modules."""
from dataclasses import asdict, dataclass

from .util import CenrepError

# ---- shared values: one concept, one number
TRF_MAX_PERIOD = 2000          # TRF's maximum period; longer units are left to the long-unit detector
FAMILY_IDENT = 0.90            # same satellite family: TRF consensuses are near-homogeneous, 90% splits variants >10% apart
SAME_VARIANT_IDENT = 0.97      # same monomer variant: Naish et al. 2021 HOR rule, <= 5 differences per 178 bp
RELATED_IDENT = 0.70           # related families (superfamily, multimer folding, shared monomer frame)
FULL_MONOMER_COV = 0.90        # a hit covering >= 90% of the model is a whole monomer
S_MIN = 20_000                 # prior: smallest satellite cluster treated as centromere-sized (bp); other sizes are rungs
LEN_TOL = 0.05                 # length tolerance for "k copies of a unit" and "same-length" comparisons
RELATED_LEN_TOL = 0.15         # length tolerance for related families (subfamilies with indels)
FULL_PERIOD_TOL = 0.01         # a TRF consensus within 1% (>= 1 bp) of its period is one whole repeat unit
MAX_CHROMOSOMES = 50           # more sequences than this: probably not a chromosome-level assembly (warning, plot)
LU_MAX_UNIT = 30_000           # longest long-unit period searched (no principled upper bound; 2 copies fit a window)

# ---- engineering and search constants (cost or tool settings; not tuned on data)
# TRF: Benson's recommended weights (match 2, mismatch 7, indel 7, PM 80, PI 10), minscore 50, max period,
# -l 6 (Mb): chunks are <= 5 Mb, so no array is longer
TRF_ARGS = ["2", "7", "7", "80", "10", "50", str(TRF_MAX_PERIOD), "-l", "6", "-h", "-ngs"]
BLAST_BATCH = 200              # seed representatives per BLAST round (bounds output size)
ANNOT_CHUNK = 2_000_000        # nhmmer chunk (bp): keeps all threads busy; E-values rescaled to chromosome length
MULTIMER_RATIO = 1.9           # a model >= 1.9 x longer than another is a multimer of it (k >= 2, minus indels)
MULTIMER_MIN_COV = 0.5         # drop a multimer hit when shorter-unit hits cover the majority of it
OVERLAP_MAX_COV = 0.2          # accept a hit only if <= 20% of it is already covered by better hits
OVERHANG_MIN_BP = 20           # terminal redundancy shorter than this is a chance match, not a run-through
MSA_MAX_GAP = 0.5              # consensus keeps alignment columns occupied in the majority of sequences
NHMMER_REPORT_E = 10.0         # hits kept up to HMMER's default reporting E 10 per search (pyhmmer returns weaker ones
                               # too); significance is then decided per hit or per run of copies (annotate.run_significant)
RDNA_TAG_EVALUE = 1e-5         # cmsearch E for an rRNA gene inside a family consensus dimer
RDNA_FLANK_EVALUE = 1e-10      # cmsearch --hmmonly E for 45S genes in array flanks (stricter: long flank queries)

# ---- retrotransposon evidence (docs/2026-10-05-retrotransposon-evidence.md)
TE_BINS = 500                  # density profile resolution: 1/500 of the chromosome
TE_MIN_BIN = 1_000             # ... but bins of at least 1 kb
TE_SMOOTH = 5                  # smoothing window, in bins
TE_MIN_ELEMENTS = 5            # fewer elements of a clade on a chromosome: no peak estimate (and one nested
                               # cluster of 2-3 elements cannot pass alone); mosses carry ~10 per chromosome
TE_DERIVED_FRAC = 0.5          # a call's main family with most array bp inside LTR retrotransposons: flag te_derived
TE_ALPHA = 1e-3                # core significance: binomial tail x (L / core length) < 0.001 (conventional level)
TE_CENTROMERIC_CLADE = "CRM"   # CR chromoviruses (CRM maize, CRR rice, CRB Brassica, CRS sugarcane) are
                               # centromere-specific in plants; other clades are often pericentromeric


@dataclass
class Params:
    trf_chunk: int = 5_000_000
    trf_overlap: int = 10_000
    min_period: int = 30
    max_period: int = TRF_MAX_PERIOD
    min_seed_span: int = TRF_MAX_PERIOD
    fam_ident: float = FAMILY_IDENT
    multi_ident: float = 0.80
    superfam_ident: float = RELATED_IDENT
    absorb_ident: float = RELATED_IDENT
    dedup_ident: float = SAME_VARIANT_IDENT
    chain_ident: float = FAMILY_IDENT
    min_fam_span: int = S_MIN
    max_families: int = 40
    hmm_sample: int = 300
    hmm_sample_bp: int = 500_000
    min_hmm_seqs: int = 20
    nhmmer_evalue: float = 1e-5
    array_gap_factor: float = 2.0
    min_array_monomers: int = 3
    evalue_ref_len: int = 0
    chrom_min_bp: int = S_MIN
    cluster_gap: int = 500_000
    min_concentration: float = 0.5
    max_subtel_frac: float = 0.5
    subtel_frac_len: float = 0.05
    max_gene_ratio: float = 0.5
    strong_prevalence: float = 0.5
    strong_cluster_bp: int = 5 * S_MIN
    rdna_flank: int = 15_000
    rdna_arrays: int = 10
    rdna_array_frac: float = 0.5
    classify_from: str = ""
    low_sat_bp: int = 5 * S_MIN // 2
    competing_ratio: float = 0.5
    fallback_min_bp: int = S_MIN // 4
    cen_assoc_frac: float = 0.5
    min_gap_n: int = 10
    diffuse_min_len: int = 5_000_000
    diffuse_max_frac: float = 0.2
    lu_flank: int = 500_000
    lu_window: int = 1_000_000
    lu_window_overlap: int = 2 * LU_MAX_UNIT
    lu_min_unit: int = TRF_MAX_PERIOD
    lu_max_unit: int = LU_MAX_UNIT
    lu_min_ident: float = 0.80
    contig_end_bp: int = 10_000
    hor_ident: float = SAME_VARIANT_IDENT
    hor_min_cov: float = FULL_MONOMER_COV
    hor_max_len: int = 200
    tandem_min_len: int = 500


# one line per field: what it does and why it has this value
REASON = {
    "trf_chunk": "TRF chunk (bp); speed only, bounded by TRF -l 6 (6 Mb)",
    "trf_overlap": "overlap between TRF chunks (bp); >= 2 x TRF_MAX_PERIOD so no array is lost at a chunk edge",
    "min_period": "shortest seed period (bp): satellites, not microsatellites (misses 5-12-bp satellites, e.g. Drosophila)",
    "max_period": "longest seed period (bp) = TRF's maximum",
    "min_seed_span": "shortest TRF array that seeds a family (bp) = one longest unit; bounds BLAST work only",
    "fam_ident": "identity for one family (seeds of equal period); see FAMILY_IDENT",
    "multi_ident": "identity when one seed is a k-mer of another: a HOR-period consensus averages divergent monomers",
    "superfam_ident": "related-family identity: superfamily table, shared monomer frame (so the HMMs), co-centromeric "
                      "and dominance tests",
    "absorb_ident": "identity for folding a k-mer family into its monomer family = RELATED_IDENT",
    "dedup_ident": "refined families this identical (same length) are one variant; see SAME_VARIANT_IDENT",
    "chain_ident": "subfamilies this identical are chained into shared arrays = FAMILY_IDENT",
    "min_fam_span": "prior: union seed span in >= 1 genome for a family to be annotated (bp) = S_MIN",
    "max_families": "families annotated per complex; cost ceiling (annotate time grows with models x genome size)",
    "hmm_sample": "monomers sampled per family for MSA/HMM (MAFFT cost)",
    "hmm_sample_bp": "seed-array bp scanned per genome when sampling monomers (cost)",
    "min_hmm_seqs": "below ~20 sampled monomers, profile emissions are noise: single-sequence HMM",
    "nhmmer_evalue": "nhmmer E-value for monomer hits (conventional)",
    "array_gap_factor": "max gap inside an array, in monomer lengths (scales with the unit)",
    "min_array_monomers": "a tandem array has >= 3 copies",
    "evalue_ref_len": "score sequences shorter than this as if this long (bp; 0 = off), so contigs get a chromosome's "
                      "nhmmer cutoff",
    "chrom_min_bp": "prior: min family bp on a chromosome / min main cluster (bp) = S_MIN",
    "cluster_gap": "prior: arrays closer than this form one cluster; also the co-centromeric distance (bp)",
    "min_concentration": "the main cluster holds the majority of a chromosome's family bp (median over chromosomes)",
    "max_subtel_frac": "majority share: of a family's bp beyond the outermost genes (terminal) or near ends, and of genomes in the subtelomeric vote",
    "subtel_frac_len": "end zone (near_end flag; near-end grade; subtelomere when a sequence has < 2 gene loci), fraction of length",
    "max_gene_ratio": "centromeres are gene-poor: main clusters at >= half the genome's gene density are not centromeric",
    "strong_prevalence": "'strong' and dominant families are centromeric on the majority of chromosomes",
    "strong_cluster_bp": "prior: min median main-cluster bp for 'strong' evidence = 5 x S_MIN",
    "rdna_flank": "window each side of a tested array searched for 45S rRNA genes (bp); fits plant 45S units",
    "rdna_arrays": "largest arrays per eligible family tested for flanking 45S genes (cost)",
    "rdna_array_frac": "the majority of tested arrays flanked by 45S genes -> rDNA spacer repeat, not centromeric",
    "classify_from": "TARGET=SRC[+SRC2][,...]: TARGET (e.g. a contig assembly) takes the centromeric families of its "
                     "SRC genomes (the --classify-from option)",
    "low_sat_bp": "prior: centromere satellite bp below this -> low_satellite flag = S_MIN x 2.5",
    "competing_ratio": "a second cluster >= half the main one -> competing_cluster flag",
    "fallback_min_bp": "prior: min cluster bp for a weak_satellite fallback call = S_MIN / 4",
    "cen_assoc_frac": "a non-primary family with the majority of its bp inside centromeres is cen_associated",
    "min_gap_n": "assemblers write gaps as runs of >= 10 N; shorter runs are ambiguity codes",
    "diffuse_min_len": "prior: min centromere length for the diffuse flag (bp)",
    "diffuse_max_frac": "diffuse flag when the satellite fraction of a long region is below this",
    "lu_flank": "long-unit scan flank around each centromere (bp)",
    "lu_window": "long-unit self-alignment window (bp); cost, >= 4 x lu_max_unit",
    "lu_window_overlap": "overlap between long-unit windows = 2 x LU_MAX_UNIT, so no unit pair is split",
    "lu_min_unit": "shortest long-unit period = TRF_MAX_PERIOD (hand-off from TRF)",
    "lu_max_unit": "longest long-unit period = LU_MAX_UNIT",
    "lu_min_ident": "min self-alignment identity for long-unit arrays: the divergence minimap2 asm20 aligns through",
    "contig_end_bp": "borrowed genomes: a call this close to a contig end is flagged contig_end (bp)",
    "hor_ident": "cd-hit-est identity for monomer variant classes = SAME_VARIANT_IDENT",
    "hor_min_cov": "only whole monomers enter the HOR search = FULL_MONOMER_COV",
    "hor_max_len": "longest HOR unit searched, in monomers (runtime bound)",
    "tandem_min_len": "shortest other tandem array reported inside centromeres (bp; reporting filter)",
}


def to_dict(p: Params) -> dict:
    return asdict(p)


def parse_borrow(spec: str) -> dict:
    """'T=S1+S2,U=S1' -> {T: [S1, S2], U: [S1]} (sources sorted); '' -> {}."""
    out = {}
    for item in filter(None, (x.strip() for x in spec.split(","))):
        tgt, _, src = (x.strip() for x in item.partition("="))
        srcs = sorted({x.strip() for x in src.split("+")})
        if not tgt or "=" not in item or not all(srcs):
            raise CenrepError(f"--classify-from: {item!r} is not TARGET=SOURCE[+SOURCE...]")
        if tgt in out:
            raise CenrepError(f"--classify-from: {tgt} is given more than once")
        out[tgt] = srcs
    return out


def format_borrow(borrow: dict) -> str:
    return ",".join(f"{t}={'+'.join(s)}" for t, s in sorted(borrow.items()))
