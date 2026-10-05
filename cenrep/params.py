"""Tunable thresholds for every stage. Each field is exposed as a --flag by the CLI."""
from dataclasses import asdict, dataclass

from .util import CenrepError


@dataclass
class Params:
    trf_chunk: int = 5_000_000
    trf_overlap: int = 10_000
    min_period: int = 30
    max_period: int = 2000
    min_seed_span: int = 2000
    fam_ident: float = 0.90
    multi_ident: float = 0.80
    superfam_ident: float = 0.70
    absorb_ident: float = 0.70
    dedup_ident: float = 0.97
    chain_ident: float = 0.90
    min_fam_span: int = 20_000
    max_families: int = 40
    hmm_sample: int = 300
    hmm_sample_bp: int = 500_000
    min_hmm_seqs: int = 20
    nhmmer_evalue: float = 1e-5
    array_gap_factor: float = 2.0
    min_array_monomers: int = 3
    evalue_ref_len: int = 0
    chrom_min_bp: int = 20_000
    cluster_gap: int = 500_000
    min_concentration: float = 0.5
    max_subtel_frac: float = 0.5
    subtel_min_bp: int = 1_000_000
    subtel_frac_len: float = 0.05
    max_gene_ratio: float = 0.5
    strong_prevalence: float = 0.5
    strong_cluster_bp: int = 100_000
    rdna_flank: int = 15_000
    rdna_arrays: int = 10
    rdna_array_frac: float = 0.5
    classify_from: str = ""
    low_sat_bp: int = 50_000
    competing_ratio: float = 0.5
    fallback_min_bp: int = 5_000
    cen_assoc_frac: float = 0.5
    min_gap_n: int = 10
    diffuse_min_len: int = 5_000_000
    diffuse_max_frac: float = 0.2
    lu_flank: int = 500_000
    lu_window: int = 1_000_000
    lu_window_overlap: int = 60_000
    lu_min_unit: int = 2_000
    lu_max_unit: int = 30_000
    lu_min_ident: float = 0.80
    contig_end_bp: int = 10_000
    hor_ident: float = 0.97
    hor_min_cov: float = 0.90
    hor_max_len: int = 200
    tandem_min_len: int = 500


HELP = {
    "trf_chunk": "TRF chunk size (bp)",
    "trf_overlap": "overlap between TRF chunks (bp)",
    "min_period": "min period of a family seed (bp)",
    "max_period": "max period of a family seed (bp)",
    "min_seed_span": "min TRF array span to seed a family (bp)",
    "fam_ident": "family identity, seeds of equal period",
    "multi_ident": "identity when one seed is a k-mer (k>=2) of the other",
    "superfam_ident": "superfamily identity (descriptive only)",
    "absorb_ident": "identity for folding a k-mer (k>=2) family into its monomer family",
    "dedup_ident": "refined families at least this identical (same length) are merged",
    "chain_ident": "subfamilies at least this identical are chained into shared arrays",
    "min_fam_span": "min union span in >=1 genome for a family to be annotated (bp)",
    "max_families": "max families annotated per complex",
    "hmm_sample": "max monomers sampled per family for MSA/HMM",
    "hmm_sample_bp": "max seed-array bp scanned per genome when sampling monomers",
    "min_hmm_seqs": "min sampled monomers for a profile HMM (else single-sequence HMM)",
    "nhmmer_evalue": "nhmmer E-value cutoff for monomer hits",
    "array_gap_factor": "max gap inside an array, in monomer lengths",
    "min_array_monomers": "min monomers per array",
    "evalue_ref_len": "score sequences shorter than this as if this long, so contigs get a chromosome's "
                      "nhmmer cutoff (0 = off) (bp)",
    "chrom_min_bp": "min family bp on a chromosome / min main-cluster bp (bp)",
    "cluster_gap": "arrays closer than this form one cluster (bp)",
    "min_concentration": "min median share of a chromosome's family bp in its main cluster",
    "max_subtel_frac": "max share of family bp inside chromosome end zones",
    "subtel_min_bp": "end zone = max(this, subtel-frac-len x chromosome length) (bp)",
    "subtel_frac_len": "end zone as a fraction of chromosome length",
    "max_gene_ratio": "max gene density in main clusters relative to genome mean",
    "strong_prevalence": "min prevalence for 'strong' evidence",
    "strong_cluster_bp": "min median main-cluster bp for 'strong' evidence (bp)",
    "rdna_flank": "window on each side of a tested array searched for 45S rRNA genes (bp)",
    "rdna_arrays": "largest arrays per eligible family tested for flanking 45S rRNA genes",
    "rdna_array_frac": "min share of tested arrays with a flanking 45S rRNA gene -> rDNA spacer repeat, not centromeric",
    "classify_from": "TARGET=SRC[+SRC2][,...]: TARGET (e.g. a contig assembly) takes the centromeric families "
                     "of its SRC genomes instead of its own per-chromosome tests",
    "low_sat_bp": "centromere satellite bp below this -> low_satellite flag (bp)",
    "competing_ratio": "2nd cluster >= this x main cluster -> competing_cluster flag",
    "fallback_min_bp": "min cluster bp for a weak_satellite fallback call (bp)",
    "cen_assoc_frac": "min share of a non-primary family's bp inside centromeres -> cen_associated",
    "min_gap_n": "min N-run length counted as an assembly gap",
    "diffuse_min_len": "min centromere length for the diffuse flag (bp)",
    "diffuse_max_frac": "diffuse flag when satellite fraction of a region >= diffuse-min-len is below this",
    "lu_flank": "long-unit scan flank around each centromere (bp)",
    "lu_window": "long-unit self-alignment window (bp)",
    "lu_window_overlap": "overlap between long-unit windows (bp)",
    "lu_min_unit": "min long-unit period (bp)",
    "lu_max_unit": "max long-unit period (bp)",
    "lu_min_ident": "min self-alignment identity for long-unit arrays",
    "contig_end_bp": "borrowed genomes: a call this close to a contig end is flagged contig_end (bp)",
    "hor_ident": "cd-hit-est identity for monomer variant classes",
    "hor_min_cov": "min model coverage for a monomer to enter the HOR search",
    "hor_max_len": "max HOR unit length in monomers (longer repeats are capped)",
    "tandem_min_len": "min length of other tandem arrays reported inside centromeres (bp)",
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
