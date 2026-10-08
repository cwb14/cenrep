# cenrep

> Draft (2026-10-06): validation numbers will be updated after the final release of the current batch.

Finds the centromeric satellite DNA of a set of related genomes de novo, and calls one centromere per chromosome
with the evidence behind it and a confidence level. You need no known monomer, no repeat library and no species
knowledge. Optional evidence from gene annotation and from LTR retrotransposons makes the calls more reliable.

Built for allopolyploid complexes (a polyploid and its progenitors, analysed together so they share one family
namespace), but it works on any set of chromosome-level genomes, or on a single one.

## What it does

1. **Tandem repeats.** TRF scans every chromosome in chunks. Results are cached by file checksum, so a genome shared
   by several runs is scanned once.
2. **Satellite families.** Tandem arrays of ≥ 2 kb are clustered into families shared by all genomes:
   - 90% identity within a family; multimers are folded into their monomer.
   - Each family gets a profile HMM built from sampled copies, founded on one whole monomer.
   - rDNA-derived repeats are tagged with Rfam models; they are never centromeric.
3. **Annotation.** nhmmer finds every monomer of every family; neighbouring monomers form arrays.
   - A copy of a short unit (< ~60 bp) can never be significant on its own, but a run of whole copies is: the run's
     E-value is the product of its copies' E-values times the chance that each next copy falls within the array gap.
   - Run-only copies never override a hit that is significant on its own.
4. **Classification, per genome.** A family is a centromere candidate when, on most chromosomes:
   - its main cluster holds most of its sequence;
   - it is not mainly near chromosome ends;
   - it is not gene-rich;
   - it is not next to 45S rDNA genes.
   On each chromosome the largest candidate cluster is chosen, preferring the genome's dominant satellite.
5. **Centromere calls.** One region per chromosome, with evidence (`satellite`, `weak_satellite`,
   `long_unit_tandem`, `retrotransposon` or `none`), a confidence level and flags.
6. **HORs.** Monomers are sorted into variant classes (cd-hit-est, 97%), and higher-order repeats are found as
   repeated runs of classes.
7. **Outputs:** GFF3, summary tables and an overview figure.

### Retrotransposon evidence (optional, recommended for plants)

Many centromeres carry no satellite: wheat and oat (CRM elements), banana, potato (about half its centromeres),
perennial soybeans, *Cryptococcus*. A satellite-only tool then picks an interstitial satellite block. If you supply
LTR retrotransposons with their clade labels, cenrep finds the clade that piles up at one locus per chromosome:
- **Finding the marker clade.** The most concentrated clade wins (a significant density core on most chromosomes:
  binomial test, Bonferroni-corrected, p < 0.001). CRM is preferred whenever it qualifies, because CRM (CR
  chromoviruses) is the centromere-specific lineage in plants.
- **Rejecting satellites.** With CRM as the marker, satellite families whose clusters lie away from the CRM core on
  most chromosomes are rejected.
- **Placing satellite-free centromeres.** A chromosome with no satellite call takes the core as its centromere
  (`evidence=retrotransposon`).
- **Moving distant satellite calls.** With CRM as the marker, a satellite call far from the core is replaced by
  the core.
- **CRM marks the centromere neighbourhood, not always the CENH3 domain itself.** In tobacco, CENH3 binds Tekay and
  avoids CRM, yet CRM clusters at the centromeres: calls are within 5 Mb of independent centromere positions
  on 44/48 chromosomes. In Solanaceae with CRM outside centromeres (eggplant), a CRM-decisive call would be wrong.
- **Other clades only flag.** When the marker is another clade, it flags agreement (`te_concordant`) or
  disagreement (`te_discordant`) but changes nothing.

Without retrotransposon input, a genome with no dominant satellite has every call capped at `low` confidence
(flag `no_dominant_satellite`): satellite evidence alone cannot tell a centromeric array from an interstitial one.

## Install

```bash
mamba env create -f environment.yml      # python, pyhmmer, pysam, pandas, matplotlib, TRF, BLAST+, MAFFT,
mamba activate cenrep                    # cd-hit, Infernal, minimap2, samtools
pip install .                            # provides the `cenrep` command
pytest -q                                # optional: about 150 tests on small synthetic genomes, ~3 min
```

## Run

```bash
samtools faidx genome1.fa                # every FASTA needs an index (bgzip .fa.gz also needs .gzi: samtools faidx does both)
cenrep run my_complex/ -t 16             # outputs in my_complex/centromere/
```

`my_complex/` holds:

| file | required | content |
|---|---|---|
| `<genome>.fa` / `.fasta` / `.fna` (`.gz` bgzip) + `.fai` | yes | one FASTA per genome, chromosome-level (≤ 50 sequences) |
| `genes/<genome>.gff3[.gz]` | no | `gene` lines (protein-coding only, if gene types are given); used to veto gene-rich candidates |
| `tes/<genome>.gff3[.gz]` | no | LTR retrotransposons (type containing `LTR_retrotransposon`) with a `clade=` attribute, e.g. from LTRquest, or EDTA + TEsorter |

Options: `-o` output dir, `-t` threads, `--trf-cache` a cache shared between runs, `--classify-from T=S` (a contig
assembly T takes its centromeric families from chromosome-level genome S), `--force STAGE` (recompute from a stage),
`-v`. A rerun resumes from the last finished stage. Every threshold is an internal constant, with its reason, in
`cenrep/params.py`.

## Outputs (`<outdir>/`)

| file | content |
|---|---|
| `<genome>_centromeric_repeat.gff` | GFF3: `centromere` (one per called chromosome), `satellite_DNA` arrays, `satellite_monomer`, `higher_order_repeat`, other `tandem_repeat`s inside the centromere |
| `<genome>_centromere_summary.tsv` | one row per chromosome: call coordinates, evidence, confidence, flags, main family, satellite bp, HOR fraction, gaps, genes, retrotransposon clade and core |
| `<genome>_hor_units.tsv` | HOR units: monomer composition, copies, tandem runs |
| `satellite_families.tsv/.fa/.hmm` | every family: class per genome (centromeric / cen_associated / other / absent), reason, consensus, HMM |
| `centromere_overview.pdf` | one panel per genome, one track per chromosome |

**Confidence.**
- `high`: a satellite call with no warning flags.
- `medium`: one warning flag, or a retrotransposon call backed by CRM, or a call that agrees with CRM.
- `low`: weak or long-unit evidence, several flags, or no dominant satellite and no retrotransposon evidence.

Flags name each reason: `low_satellite`, `competing_cluster`, `near_end`, `diffuse`, `gene_rich`, `gaps:N`,
`te_concordant`, `te_discordant`, `satellite_away_from_core`, `no_dominant_satellite`.

## Validation

- **Synthetic genomes** with planted satellites, HORs, a decoy, a gap, a long-unit array and a satellite-free
  chromosome are run by the test suite.
- **Literature known answers.** Reproduced:
  - *A. thaliana* CEN180 on Col-CEN: 67,241 monomers, against 66,131 in Naish et al. 2021 (v0.3.7).
  - pAa (*A. arenosa*), CentBr (*Brassica* A), CentO (*Oryza* B), CentBd/CentBs (*Brachypodium*) and Cent728
    (*Mimulus*) on every chromosome.
  - Correct negatives: the satellite-free *Brassica* B genome and *Oryza* C CEN2.
- **CENH3 ChIP-seq**, the direct test, on *Aegilops geniculata* (Yunusbaev et al. 2026, same assembly):
  - The CRM core overlaps the CENH3 peak on 14/14 chromosomes.
  - Calls overlap on 10/14 (v0.7.0 rules replayed on the v0.6.2 run; full rerun pending); the other 4 sit on
    satellite arrays 2–11 Mb away, inside wide CRM cores.

## Limitations

- Units < 30 bp (*Drosophila*) and > 30 kb are outside the search. Units of 30–60 bp are found as runs of copies.
  Arrays of a composite unit (two different short units in tandem, e.g. 48 + 39 bp) are described by the composite.
- Non-tandem centromeres are reported as `none` (no retrotransposon input) or by their retrotransposon core.
- Without retrotransposon input, a dominant but non-centromeric satellite is called centromeric (oat: the 335–360-bp
  family on every chromosome). The `no_dominant_satellite` cap catches only genomes where no family dominates; give
  `tes/` input for grasses and other genomes with retrotransposon-based centromeres.
  Holocentric chromosomes have no diagnostic yet.
- At most 40 families per run, the largest. Satellite-rich genomes (oat, wheat relatives) drop hundreds of small
  families; the centromeric ones tested so far were kept.
- Solanaceae: short (150–190-bp) satellites are often pericentromeric, while centromeres are kb-unit arrays or
  satellite-free. CENH3 ChIP on *S. verrucosum* (same assembly) puts calls within 1 Mb on 7/11 chromosomes but
  2–10 Mb off on 4, some at high confidence.
- Under-assembled centromeres give short calls. Lengths reflect the assembly.
- In large grass genomes, calls on small satellite arrays inside wide CRM cores can be offset from CENH3 by a few Mb.
- Chromosome ends: a family is subtelomeric when most of it lies beyond the outermost genes (needs `genes/`; else
  the terminal 5% of each sequence). An acrocentric centromere next to an end is kept but flagged `near_end`
  (*Fragaria* Chr1D, CENH3-bound, 1.0 Mb from the end). Gene models reaching into a subtelomere can let a small
  subtelomeric family through (flagged `near_end`).
