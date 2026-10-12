# Retrotransposon evidence for centromere calls (design, 2026-10-05)

Status: current in v0.8.9. Later changes are folded in: element floor 5 and sequence-based kinds (v0.8.5), LTR
overlap reporting (v0.8.8–v0.8.9, docs/2026-10-08-family-roles.md).

## Why
- Many centromeres are not satellite-based: Triticeae and *Avena* (CRW/Cereba CRM elements), *Musa* (Nanica LINE +
  CRM), perennial *Glycine*, most *Nicotiana* centromeres, half of potato's, *Cryptococcus* (Tcn elements).
- On such genomes cenrep v0.5 still calls the largest eligible satellite cluster. In *Aegilops ventricosa*, 27 of 28
  calls sit on interstitial satellite blocks 22–409 Mb from the centromere, mostly at high/medium confidence.
- Centromeric retrotransposons concentrate at one locus per chromosome, and the clade differs by lineage (Athila in
  *A. thaliana*, CRM in grasses, Ale/CRM in *Brassica*). Discovered per genome, the most concentrated clade's peak
  agrees with validated satellite calls (*A. thaliana* Col-CEN within 0.2–0.9 Mb, *B. rapa* within 0–1.8 Mb,
  *O. minuta* B within 0.3–3.6 Mb). It also exposes the *Aegilops* calls (CRM peaks at the known wheat D
  centromere positions).

## Input (optional; no CLI option)
`GENOME_DIR/tes/<genome>.gff3[.gz]`: LTR retrotransposons with a clade label. A feature counts when its type
contains `LTR_retrotransposon` and it has a `clade=` attribute (case-insensitive key; values unknown/mixture/NA are
ignored). Producers: LTRquest; EDTA or LTR_retriever + TEsorter (write the TEsorter clade as `clade=`).
Without the file, cenrep is a satellite-only tool and the honesty rule below applies.

## Rules
1. **Profile per clade and chromosome.** Element midpoints are counted in bins of length L/500 (resolution: 1/500 of
   the chromosome, at least 1 kb) and smoothed over 5 bins. The **core** is the run of bins around the maximum
   whose smoothed count is at least half the peak count (full width at half maximum).
2. **Significance.** With k of the chromosome's n elements of the clade in a core of length c, the core is
   significant if the binomial tail P(X ≥ k; n, c/L) × (L/c) < 1e-3. The L/c factor is a Bonferroni correction for
   the number of places a core could sit; 1e-3 is the conventional level. Only chromosomes with ≥ 5 elements of the
   clade are tested (v0.8.5; was 10, which left most chromosomes of the *Sphagnum troendelagicum* B haplotype
   untested although 8–17 Bryco elements cluster at one locus on each; 5 keeps a single nested cluster of 2–3
   elements from passing alone).
3. **Marker clade.** Candidates are clades with a significant core on the majority of chromosomes.
   - If CRM is a candidate, it is the marker and its evidence is **decisive**. CR chromoviruses are the
     centromere-specific LTR lineage of plants: CRM in maize, CRR in rice, CRB in *Brassica*, CRS in sugarcane.
   - Otherwise the most concentrated candidate (highest median core share) is the marker, and its evidence only
     flags. Other clades are often pericentromeric. In *B. juncea* A only Retand/Athila qualify (its centromeres are
     collapsed in the assembly), and their cores sit beside CentBr.
   - Rejected alternatives:
     - The most concentrated clade alone: *B. juncea* A loses CentBr.
     - Every candidate clade must disagree: Triticeae chromosomes give 6–8 "significant" clades, so nothing is ever
       rejected.
     - A two-fold concentration margin: it labels *Aegilops* genomes inconsistently, e.g. CRM 0.39 vs 0.23.
4. **Agreement.** A call agrees with the core if the gap between them is at most the wider of the two. Two estimates
   of one locus; centromeric retrotransposons are partly pericentromeric, so cores can sit beside the satellite array.
5. **Satellite families on trial (classify, decisive evidence only).**
   - Centromeric family kinds are tested: sequence-related families (same superfamily at ≥ 70% identity, or one
     chain of close subfamilies: CentBr1/2, CentO 155/165) count together, so a subfamily on a few chromosomes is
     not rejected alone. v0.8.5: period alone no longer makes families kin. Before, the *Populus* 148-bp
     knob repeat (no sequence similarity) was chained with the 144/156-bp centromeric satellites, and its
     rejection took them along.
   - A chromosome agrees if any member's main cluster agrees with its core. A kind whose clusters disagree on the
     majority of its chromosomes is made ineligible: "main clusters away from the CRM retrotransposon centromere
     (k/n)".
   - Families are then chosen again, until no kind disagrees. A family that was outcompeted by a rejected one gets
     its chance: the 553-bp CentT550-like satellite in *Aegilops*.
6. **Calls.**
   - A satellite or long-unit call that agrees gets flag `te_concordant`. With decisive evidence it is raised to at
     least `medium`: two independent lines agree.
   - One that disagrees gets `te_discordant` and, if the evidence is not decisive, keeps its position.
     With decisive evidence (CRM) the core replaces it (flag `satellite_away_from_core`; v0.7.0). CENH3 ChIP on
     *Ae. geniculata* (same assembly) puts the centromere on the CRM core, not on the distant array: 2U 30 Mb
     away, 7M 268 Mb away. Validated satellite genomes keep every call except *O. officinalis* C08.
   - A chromosome with no call, or with a disagreeing long-unit call, gets the core as its centromere: evidence
     `retrotransposon`, confidence `medium` (decisive) or `low`.
7. **Honesty rule without TE input.** If no family is dominant (strong, and its kind is the largest cluster on the
   majority of chromosomes), every call is capped at `low`, with flag `no_dominant_satellite`. Satellite evidence
   alone cannot tell a centromeric satellite from an interstitial block on such genomes.

## Outputs
- regions/summary: `te_clade`, `te_core_start`, `te_core_end`.
- GFF centromere attributes: `te_clade`, `te_core`.
- Family table: `dominant` column (classify); rule-5 rejections appear as ineligible, class `other`, with their reason.

## Versions (when introduced)
classify (dominant column) and centromeres (rules 1–7) change → STAGE_VERSION 0.6.0 for both; rerun from classify.

## Validation (sanity, not tuning)
- Synthetic: TE cluster at the satellite centromere (concordant); TE cluster elsewhere (family rejected, TE call
  made); uniformly spread clade (no evidence).
- Real data:
  - An earlier test set (*A. thaliana* Col-CEN, *O. minuta*, *B. rapa*): satellite calls stay, concordant.
  - *Aegilops ventricosa*: interstitial families rejected.
  - *Ae. geniculata*: calls checked against the published CENH3 ChIP peaks on the same assembly (Yunusbaev 2026).

## Real-data check (classification and calls rerun on existing arrays; rDNA flank test and long-unit scan skipped)
| genome | marker | outcome |
|---|---|---|
| *A. suecica* complex, earlier test set (Col-CEN *A. thaliana*) | Athila (flags) | unchanged: CEN180/pAa calls all kept |
| *O. minuta* B / C, *O. officinalis* C, *O. punctata* B | none / CRM | unchanged; the satellite-free C-genome CEN2 gets a retrotransposon call |
| *B. rapa* A | CRM | unchanged: 10/10 high |
| *B. juncea* A | Retand (flags) | CentBr kept |
| *B. nigra* B, *B. juncea* B | CRM | 492/71-bp families rejected; 8/8 retrotransposon calls (literature: satellite-free B genome) |
| *Brachypodium* 4 genomes | CRM / Tekay | unchanged |
| *Mimulus* 4 genomes | CRM / Angela | unchanged |
| *Aegilops ventricosa* 4 genomes | CRM | interstitial families rejected; the 553-bp CentT550-like satellite is centromeric on 5–6/7 per genome |

## CENH3 test (Ae. geniculata, Yunusbaev 2026 ChIP-seq on the same assembly)
- CRM core alone: overlaps the CENH3 peak on 14/14 chromosomes.
- v0.6.2 calls: 8/14 overlap. v0.7.0 rule (decisive core replaces a disagreeing satellite call): 10/14 overlap; the other 4
  are 2.0-10.9 Mb off, satellite arrays inside wide (8-19 Mb) CRM cores that count as agreeing.
- Rejected alternative: "midpoint of either interval inside the other" as the agreement test gives 13/14 here, but
  overrides validated CENH3-bound satellite calls (CentO-C 4/8 in O. minuta C, CentBs, Cent728, CentO).
