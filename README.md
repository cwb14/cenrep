# cenrep

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
pytest -q                                # optional: about 160 tests on small synthetic genomes, ~3 min
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
| `genes/<genome>.gff3[.gz]` | no | `gene` lines (protein-coding only, if gene types are given); used to veto gene-rich candidates and to define subtelomeres (beyond the outermost gene) |
| `tes/<genome>.gff3[.gz]` | no | LTR retrotransposons (type containing `LTR_retrotransposon`) with a `clade=` attribute, e.g. from LTRquest, or EDTA + TEsorter; used for the marker clade and for `in_ltr` |

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
| `<genome>_satellite_arrays.gff` | GFF3: every array of every family, on every sequence (e.g. ~21,000 in *Poa annua* A): `family` (as in `satellite_families.tsv`), `role`, `rdna` where tagged, monomer length and count, identity, `in_centromere` |
| `satellite_families.tsv/.fa/.hmm` | every family: class, role and (with `tes/`) `in_ltr` / `ltr_clade` per genome, reason, consensus, HMM |
| `centromere_overview.pdf` | one panel per genome, one track per chromosome |

**Family roles.** The tests that decide what is centromeric also say what the other families are. Each family gets
one role per genome (first match wins):
- `rDNA`: 5S rRNA (Rfam match on the consensus; `rdna=5S_rRNA`) or 45S spacer repeat (largest arrays among 45S
  rRNA genes, in this or another genome of the run; `rdna=45S_spacer`).
- `centromeric`, `cen_associated` (most bp inside centromere calls, never the main family).
- `subtelomeric`: most bp beyond the outermost genes or in the terminal 5% of a sequence, and either so in at least
  half of the genomes or in blocks at the ends of most of this genome's chromosomes. An interior genome is never
  overruled. A family terminal in one genome on one chromosome is not enough: truncated assembly ends can make a
  centromeric array look terminal.
- `interstitial`: concentrated at one locus per chromosome, but not the centromere. Knobs and other heterochromatic
  blocks land here; cenrep cannot tell them apart.
- `dispersed`: spread along the chromosomes. `minor`: < 20 kb in this genome, too little to judge.
- `unassessed`: contig assemblies (`--classify-from`), where per-chromosome tests do not apply.
Not reported: telomeres and other units < 30 bp (outside the search), 45S/5S genes themselves, and tandem repeats
outside the ≤ 40 families kept per run.

**Confidence.**
- `high`: a satellite call with no warning flags.
- `medium`: one warning flag, or a retrotransposon call backed by CRM, or a call that agrees with CRM.
- `low`: weak or long-unit evidence, several flags, or no dominant satellite and no retrotransposon evidence.

Flags name each reason: `low_satellite`, `competing_cluster`, `near_end`, `diffuse`, `gene_rich`, `gaps:N`,
`te_concordant`, `te_discordant`, `satellite_away_from_core`, `no_dominant_satellite`, `te_derived` (the main family
has most of its array bp inside clade-classified LTR retrotransposons; needs `tes/`).

## Validation

Results of v0.8.9 on 185 genomes from 38 polyploid complexes (plants, fish, fungi): 93% of 2,331 chromosomes
called (1,771 satellite, 301 retrotransposon, 97 long-unit tandem; 162 `none`).

- **Synthetic genomes** with planted satellites, HORs, a decoy, a gap, a long-unit array and a satellite-free
  chromosome are run by the test suite (164 tests).
- **CENH3 data**, the direct test:
  - *Aegilops geniculata* (CENH3 ChIP-seq peaks, Yunusbaev et al. 2026, same assembly): calls overlap the CENH3
    domain on 10/14 chromosomes. The other 4 sit on satellite arrays 2–11 Mb away, inside wide CRM cores.
  - *Solanum verrucosum* (Zhang et al. 2014 CENH3 ChIP reads, mapped to the same assembly): 11/12 chromosomes
    called, 7/12 within 1 Mb of the CENH3 domain. Its short satellites are mostly pericentromeric (see
    Limitations).
  - *Fragaria × ananassa* (CENH3 ChIP, Jin et al. 2025): 28/28 centromeres on the published satellite type,
    including the Type II (Chr1D) and Type III (Chr1C) exceptions.
  - *Brachypodium* (Chen et al. 2024, same assemblies): CentBd/CentBs found and kept apart by subgenome; monomer
    counts within 2% of the published ones (*B. hybridum* D 18,100 vs 17,838; S 87,200 vs 86,896).
  - *Nicotiana tabacum*: 44/48 chromosomes placed by CENH3 ChIP-clone clusters or by CENH3/Hi-C centromere windows
    transferred from a near-T2T *N. sylvestris* assembly.
- **Literature known answers**, recovered on every or nearly every chromosome: CEN180 (*A. thaliana*), pAa
  (*A. arenosa*), CentBr (*Brassica* A/C), CL1 (radish), CentO and CentO-C (*Oryza*), Cent728 (*Mimulus*),
  the 166-bp *Prunus* satellite, the *Carassius* centromere probe.
  - *A. thaliana* CEN180 on Col-CEN: 67,241 monomers vs 66,131 in Naish et al. 2021 (v0.3.7).
  - Correct negatives: the satellite-free *Brassica* B genomes, oat A/D (CRM centromeres), mosses (*Sphagnum*,
    Copia centromeres), *Cryptococcus* (`none`).
- **Held-out test** (5 complexes; pass/fail criteria written before they were run, outputs unseen while the code
  changed):
  - *Saccharum* (SCEN, positive): pass, 36/36 chromosomes, all high confidence.
  - *Magallana* (no uniform centromere sequence): pass.
  - *Glycine* and *Musa* (CENH3/FISH: no centromeric satellite): fail. A tandem repeat inside the retrotransposon
    centromeres is reported as the satellite. The location agrees with the CRM core; the calls are flagged
    `te_derived`.

## Design notes

The reasoning behind the rules, with the cases that motivated them, is in `docs/`:
- `2026-10-05-retrotransposon-evidence.md`: marker clade, CRM decisiveness, rejecting and moving satellite calls.
- `2026-10-07-short-units.md`: runs of copies for units shorter than ~60 bp.
- `2026-10-08-subtelomeres.md`: subtelomeric families across genomes; rDNA tags.
- `2026-10-08-family-roles.md`: family roles, the all-arrays GFF, LTR overlap (`in_ltr`, `te_derived`).
- `literature_cen_satellites.md`: published centromeric satellites used as known answers.

Every threshold is in `cenrep/params.py` with a one-line reason.

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
  satellite-free. CENH3 ChIP on *S. verrucosum* (same assembly): of 11 called chromosomes, 7 are within 1 Mb of
  CENH3 and 4 are 2–10 Mb off, some at high confidence.
- Under-assembled centromeres give short calls. Lengths reflect the assembly.
- In large grass genomes, a satellite array at the edge of a wide CRM core (8–19 Mb in *Aegilops*) is accepted and
  can be 2–11 Mb from CENH3 (4/14 *Ae. geniculata* chromosomes). The tolerance stays because CRM elements flank, rather
  than fill, true satellite centromeres (rice CRR beside CentO): requiring overlap would move CentO calls off
  their arrays.
- Without CRM, retrotransposon evidence only flags. Where the assembly lacks a chromosome's centromeric array, a
  non-centromeric satellite elsewhere on it can be called at high confidence (`te_discordant`; *B. napus* A06/A09
  on an AT-rich 156-bp repeat). Calls are not yet checked for consistent position between homologous genomes.
- Chromosome ends: a family is subtelomeric when most of it lies beyond the outermost genes (needs `genes/`; else
  the terminal 5% of each sequence), judged across genomes (see Family roles). An acrocentric centromere next to an
  end is kept but flagged `near_end` (*Fragaria* Chr1D, CENH3-bound, 1.0 Mb from the end). Gene models reaching
  into a subtelomere can still let a subtelomeric family through when no genome has it beyond the outermost gene
  (flagged `near_end`).
  Interstitial blocks a few hundred kb from an end, with genes beyond them, are not subtelomeric by this rule
  and can be called (*B. rapa* Chiifu v3: CL25-like AT-rich 330–350-bp repeat; `near_end`, often `gene_rich`).
- Tandem repeats that are part of retrotransposons, or derived from them, are reported as satellites. In
  retrotransposon centromeres the call is then in the right place (it agrees with the CRM core), but labelled
  `satellite`. Held-out tests (fixed in advance): *Glycine* (CENH3: no centromeric satellite) gets a 1.8-kb
  tandem family, and *Musa* A a 62-bp repeat from inside SIRE elements, as dominant families.
  - With `tes/`, `satellite_families.tsv` gives each family's share of array bp inside clade-classified LTR
    retrotransposons (`in_ltr`, `ltr_clade`), and calls on such families are flagged `te_derived`.
  - Elements without a clade are ignored: kb-unit tandem arrays look like LTR pairs to structural detectors.
  - Tandem arrays built from TE fragments (*Coffea* CRM/Tekay units, the published centromeres there) are
    mostly not inside annotated elements, so they are rarely flagged.
  - The flag is a report, not a verdict.
- Retrotransposon evidence needs ≥ 5 intact elements of a clade per chromosome. Genomes with few intact LTR-RTs
  (mosses: ~10 per chromosome) are near this limit.
