# cenrep parameter audit (v0.4.3)

date: 2026-10-05 | scope: the 51 `Params` fields (all exposed as CLI flags) plus the magic constants outside `params.py`
read: `cenrep/*.py`, `tests/conftest.py`, `docs/2026-09-25-centromere-design.md`, `.sdd/progress.md`, `tools/*.sh`
run: nothing from the pipeline or the test suite. The only computation was a set of read-only summaries of the v0.3.7 outputs
of the 10 allo10 complexes (`/data2/chris/tetraploid/allo10/clean/*/centromere/`) and of `.fai` lengths. No held-out complex was used.

Yardstick (owner's rules for portable tools):
- "Few options. Before adding a CLI option, ask whether it is truly needed. Best: no option and no hardcode (the value is derived
  from the data or fixed by biology). Next: an option with a principled default. Worst: a constant tuned to this dataset. Fixed
  internal constants sit in one place, each with a one-line reason. Keep `--help` short."
- "The dataset is a test case, not a training set. Never tune cutoffs until the output looks right here; prefer thresholds that
  scale with each input (e.g., relative to that genome's own background)."
- "Choose rules and thresholds by principle (conventional evidence levels, data-derived or cross-validated choices, explicit
  priors), never by which value best matches known answers."
- "Generic I/O ... Nothing hard-coded to this dataset (paths, species, chromosome names, sizes, counts)."
- "Discover, don't presuppose: learn what the method needs from each input so it works on an unstudied species."

Categories: **A** engineering knob -> internal constant | **B** fixed by biology/convention -> internal constant + reason |
**C** absolute size that should come from each input | **D** classification ratio/fraction | **E** genuine user input.
Tally: A 13, B 13, C 11, D 13, E 1. Proposal: 0 threshold flags, 6 user options (section 5).

## Evidence gathered for this audit

From the v0.3.7 outputs: 40 genomes, 371 called regions, chromosomes 10–94 Mb.
- **Satellite amount does not scale with chromosome length.** Log-log slope of centromere satellite bp against chromosome length
  = -1.6 (r = -0.45, n = 371). Lineage dominates: per-complex median centromere satellite is 88 kb in *B. juncea* and 4.7 Mb in
  *Capsella*. So a rule "fraction of chromosome length" is not supported for satellite amounts.
- **No natural break in array gaps.** Gaps between consecutive primary-family arrays (n = 30,183, pooled) have one mode at 3–10 kb
  and a smooth tail out to 30 Mb, with no second mode. A gap cutoff therefore cannot be learned from a valley. It is a convention.
- **The calls are sensitive to `cluster_gap`.** I re-clustered the same arrays with cenrep's own `make_clusters`. At 250 kb or
  1 Mb instead of 500 kb, the region length changes by more than 10% for 25–26% of the 361 satellite calls, and more than 2-fold
  for 9–10%. At 125 kb the figures are 39% and 17%, and 7 calls jump to a different cluster. Inside called regions, 74 gaps fall
  at 400–499 kb, just under the cutoff. 100 of 288 regions have another primary array 0.5–1 Mb beyond their edge.
- **Chromosome size range in the batch** (33 non-held-out complexes, 2,052 sequences): 0.73–724 Mb, about 1,000-fold.
  47 chromosomes are ≤ 2 Mb (all *Cryptococcus*). 300 are 2–20 Mb (Lemna, Populus, Sphagnum, Actinidia, Mimulus, Arabidopsis,
  fish).
- **The Cryptococcus run reached the right answer for the wrong reason.** It found 4 TRF seeds ≥ 2 kb and 0 families ≥ 20 kb, so
  every chromosome was called `none`. That answer comes from an abundance cutoff, not from any test that ran.
- **The test suite already fails on the defaults.** `tests/conftest.py` `TEST_PARAMS` overrides 14 size fields so that a 300-kb
  synthetic chromosome works. Acceptance criterion for the C rules below: the synthetic complex runs on defaults. `trf_chunk` may
  stay overridden, since it only affects speed.

## 1. Parameter table

Line numbers are per file (v0.4.3). "→" means the parameter is passed on to that line.

| name | default | stage | what it controls | where used | cat | proposal |
|---|---|---|---|---|---|---|
| trf_chunk | 5,000,000 | trf | chromosome chunk size for TRF | trf.py:90,144 | A | internal; coupled to TRF `-l 6` (trf.py:19), so derive `-l` from the chunk size or fix both together; part of the TRF cache key |
| trf_overlap | 10,000 | trf | overlap between TRF chunks | trf.py:90,144 | A | internal; must be ≥ 2 × TRF MaxPeriod (4 kb) |
| min_period | 30 | families | shortest seed period (excludes microsatellites) | families.py:24 | B | internal: "satellites, not microsatellites". Blind spot: short centromeric satellites (*Drosophila* 5–12 bp; Chang 2019 PLoS Biol). Document it |
| max_period | 2,000 | families | longest seed period | families.py:24 | B | delete as a separate value: it duplicates MaxPeriod in TRF_ARGS (trf.py:19). Raising it does nothing; lowering it opens a gap below lu_min_unit. Make one `TRF_MAX_PERIOD` |
| min_seed_span | 2,000 | families | min TRF array span to become a seed | families.py:25 | A | internal, = TRF_MAX_PERIOD (one longest unit). It only bounds BLAST work; abundance is judged later |
| fam_ident | 0.90 | families | identity to join a family (equal periods) | families.py:426 (→ :38-53) | D | internal. Reason line: "TRF consensuses are near-homogeneous; 90% splits variants > 10% apart". Note: the value was chosen to keep known subfamilies apart (CentBd/CentBs 87%, CentBr1/2 ≤ 86.9%) |
| multi_ident | 0.80 | families | identity when one seed is a k-mer of the other | families.py:426 (→ :75-83) | D | internal: "a HOR-period consensus averages divergent monomers" |
| superfam_ident | 0.70 | families | superfamily grouping **and** frame-rotation identity | families.py:400,407,479; read by classify.py:57 | D | internal. **The help text "descriptive only" is false.** The value sets each consensus frame (and so the HMMs), and superfamily membership drives co-centromeric and dominance decisions through `classify.related`. Merge with absorb_ident |
| absorb_ident | 0.70 | families | fold a k-mer family into its monomer family | families.py:432 | D | internal, = superfam_ident (one `RELATED_IDENT`) |
| dedup_ident | 0.97 | families | merge near-identical refined families | families.py:465 | B | internal. Anchor: Naish 2021 HOR rule (≤ 5 differences per 178 bp), = hor_ident, giving one `SAME_VARIANT_IDENT`. The ledger rejected 0.95 because it merged Bjuncea CentBr1/2; record the anchor, not that outcome |
| chain_ident | 0.90 | families | subfamilies chained into shared arrays | families.py:481 | D | internal, = fam_ident (same "same family" level) |
| min_fam_span | 20,000 | families | union seed span in ≥ 1 genome for a family to be annotated | families.py:439,443 | C | coverage rule R7, or at least the same prior as chrom_min_bp (R2) |
| max_families | 40 | families | cap on annotated families | families.py:441,443 | A | engineering ceiling, because annotate cost grows with #HMMs × genome size. Sized on 10 plant complexes; satellite-rich Triticeae may hit it. Keep the "dropped by max_families" log line, and use R7 |
| hmm_sample | 300 | families | monomers sampled for MSA/HMM | families.py:368,392 | A | internal (MAFFT cost) |
| hmm_sample_bp | 500,000 | families | seed-array bp scanned per genome for sampling | families.py:344,349 | A | internal; optionally derive as ≥ hmm_sample × 3 × consensus length (500 kb holds only 250 copies of a 2-kb unit) |
| min_hmm_seqs | 20 | families | min sampled monomers for a profile HMM | families.py:403,410 | B | internal: "below ~20 sequences, emission estimates are noise; use a single-sequence HMM" |
| nhmmer_evalue | 1e-5 | families (+annotate) | E-value cutoff for monomer hits (sampling and annotation) | families.py:379,381; annotate.py:166 (→ :73,80) | B | internal. Make the search space the whole genome (R4), so the cutoff no longer depends on chromosome or contig length |
| array_gap_factor | 2.0 | annotate | max gap inside an array, in monomer lengths | annotate.py:136,150 | B | internal. It already scales with monomer length; this is the style the C rows should follow |
| min_array_monomers | 3 | annotate | min monomers per array | annotate.py:136,153 | B | internal: "a tandem array = ≥ 3 copies" |
| evalue_ref_len | 0 | annotate | score short sequences as if this long | annotate.py:166 (→ :73) | C | delete via R4. Interim: set automatically for --classify-from targets to the median source chromosome length. The batch borrows nothing, so this has no effect there |
| chrom_min_bp | 20,000 | classify | **6 jobs**: family bp on a chromosome (prevalence); family genome total; min main cluster; satellite vs weak_satellite; long-unit cluster; borrowed `defines` | classify.py:76,90,109,174; centromere.py:25,47 | C | split by job. One named prior S_MIN for "centromere-sized cluster", with the rungs derived from it (R2) |
| cluster_gap | 500,000 | classify | arrays closer than this form one cluster; also the co-centromeric distance | classify.py:74,125; centromere.py:23,45 | C | density-segment rule R1. It is the most influential size in the tool (see evidence) |
| min_concentration | 0.5 | classify | median share of a family's per-chromosome bp held by its main cluster | classify.py:92-93; centromere.py:47 | D | internal: "the main cluster holds the majority" |
| max_subtel_frac | 0.5 | classify | max share of family bp in end zones | classify.py:94-95; centromere.py:47 | D | internal (majority); make the test two-ended (R3) |
| subtel_min_bp | 1,000,000 | classify | absolute floor of the end zone | classify.py:37 (→ :40-44; centromere.py:47,71-73) | C | delete the floor (R3). At L ≤ 2 Mb it makes the whole chromosome an end zone; at 2–20 Mb it covers 10–100% of the chromosome |
| subtel_frac_len | 0.05 | classify | end-zone width as a fraction of chromosome length | classify.py:37 | B | internal; the only end-zone parameter left after R3 |
| max_gene_ratio | 0.5 | classify | max gene density in main clusters vs the genome mean | classify.py:96-97; centromere.py:50 | D | internal; use one baseline (chromosome mean) everywhere (M4, section 2) |
| strong_prevalence | 0.5 | classify | **2 jobs**: `strong` label; dominant-family precedence share | classify.py:114-115,134 | D | internal (majority). Compute prevalence length-weighted, since unplaced scaffolds and organelles deflate a per-sequence count |
| strong_cluster_bp | 100,000 | classify | median main-cluster bp for `strong`/dominant | classify.py:114,134 | C | a rung of R2 (5 × S_MIN). It decides precedence, so keep it tied to S_MIN until there is evidence |
| rdna_flank | 15,000 | classify | window searched for 45S genes around arrays | pipeline.py:204 (→ families.py:181-205) | C | derive from the genome (R5). 15 kb fits plant 45S units (9–11 kb), not the 43-kb human-type units |
| rdna_arrays | 10 | classify | largest arrays tested per family | pipeline.py:204 | A | internal cost bound; R5 removes it |
| rdna_array_frac | 0.5 | classify | share of tested arrays next to 45S genes that marks an rDNA spacer | classify.py:149; pipeline.py:206 | D | internal (majority) |
| classify_from | "" | classify | contig assemblies take their centromeric families from named sources | pipeline.py:61-62; classify.py:155-182; centromere.py:95-107 | E | keep as an optional override. Default: auto. A genome that is not chromosome-level borrows from all chromosome-level genomes in the dir |
| low_sat_bp | 50,000 | centromeres | `low_satellite` flag | centromere.py:67 | C | relative: satellite bp < 1/10 of this genome's median called-centromere satellite bp (catches collapsed centromeres such as Bjuncea_A), with S_MIN as floor |
| competing_ratio | 0.5 | centromeres | `competing_cluster` flag | centromere.py:69 | D | internal |
| fallback_min_bp | 5,000 | centromeres | min cluster for a `weak_satellite` call | centromere.py:27 | C | R2: enrichment over the family's own dispersed background, or S_MIN/4 as a pure refactor |
| cen_assoc_frac | 0.5 | centromeres | share of a non-primary family inside centromeres that makes it `cen_associated` | classify.py:193 | D | internal (majority) |
| min_gap_n | 10 | centromeres | N-run length counted as an assembly gap | centromere.py:65 | B | internal: "assemblers write gaps as runs of ≥ 10 N; shorter runs are ambiguity codes" |
| diffuse_min_len | 5,000,000 | centromeres | min region length for the `diffuse` flag | centromere.py:74 | C | absorbed by R1. Fitted on one complex (only Carabica is affected) |
| diffuse_max_frac | 0.2 | centromeres | satellite fraction below which a long region is diffuse | centromere.py:74 | D | becomes ρ of R1 (one density constant) |
| lu_flank | 500,000 | centromeres | long-unit scan flank around each call | centromere.py:106 | A | **effectively dead for final outputs**: GFF/summary keep only long-unit arrays inside the call (output.py:125). Set it to 2 × lu_max_unit, or report pericentromeric long units on purpose |
| lu_window | 1,000,000 | centromeres | minimap2 self-alignment window | longunit.py:105 | A | internal (cost); must be ≥ 4 × lu_max_unit |
| lu_window_overlap | 60,000 | centromeres | overlap between windows | longunit.py:105 | A | derive = 2 × lu_max_unit (the current value already equals this) |
| lu_min_unit | 2,000 | centromeres | min long-unit period | longunit.py:66,112 | B | derive = TRF_MAX_PERIOD (hand-off from TRF) |
| lu_max_unit | 30,000 | centromeres | max long-unit period | longunit.py:66 | A | internal. There is no principled upper bound, so say that; 2 copies must fit in a window |
| lu_min_ident | 0.80 | centromeres | min self-alignment identity | longunit.py:66 | B | internal: matches the divergence the minimap2 asm20 preset aligns through (longunit.py:56) |
| contig_end_bp | 10,000 | centromeres | borrowed genomes: a call this close to a contig end is flagged | centromere.py:58,77 | C | derive from the cluster tolerance (R1): more satellite beyond the end would have joined the call. No effect on the batch |
| hor_ident | 0.97 | hor | cd-hit-est identity for monomer variant classes | hor.py:205 (→ :175) | B | internal (Naish 2021), = dedup_ident. Derive cd-hit `-n` safely (M7: < 0.88 is invalid with `-n 8`) |
| hor_min_cov | 0.90 | hor | min model coverage for HOR monomers | hor.py:192 | B | internal, = the sampling coverage at families.py:384 (one `FULL_MONOMER_COV`) |
| hor_max_len | 200 | hor | max HOR unit length (monomers) | hor.py:236 | A | internal runtime bound (O(n × max_len)) |
| tandem_min_len | 500 | output | min length of other tandem repeats reported inside centromeres | output.py:125 | A | internal reporting filter; users can filter the GFF |

### Data-derived rules for the C rows (with the evidence available)

- **R1. Cluster scale** (replaces cluster_gap, diffuse_min_len, diffuse_max_frac and contig_end_bp; also gives the co-centromeric
  distance).
  - Rule: a family's main cluster, and the centromere, is the maximal segment whose satellite fraction stays ≥ ρ. It is found
    with the maximal-scoring-subsequence algorithm (Ruzzo & Tompa 1999), scoring +(1-ρ) per bp inside arrays and -ρ per bp
    outside.
  - The tolerated gap then scales with the arrays on either side. A 500-kb ATHILA island inside 2 Mb of CEN180 is bridged; a lone
    20-kb array 500 kb away is not. This holds on a 1-Mb fungal chromosome and on a 700-Mb *Aegilops* chromosome alike.
  - ρ is one dimensionless constant; the current diffuse bound of 0.2 is the obvious candidate.
  - Why not a learned gap cutoff: the gap distribution has no valley (see evidence).
  - Unclear: the value of ρ. Evidence that would settle it: region sizes at ρ = 0.1/0.2/0.3 on the 10 validated complexes against
    published sizes (*At* 2.2–3.6 Mb; Brassica A 2.74–13.12 Mb; *B. stacei* 1.40–2.50 Mb; *O. punctata* 1.23–1.67 Mb). Use this as
    a sanity check that the rule behaves, not to pick the best-matching ρ.
- **R2. Abundance ladder** (min_fam_span, chrom_min_bp, fallback_min_bp, low_sat_bp, strong_cluster_bp). These five are one scale:
  5k / 20k / 50k / 100k = 0.25 / 1 / 2.5 / 5 × 20 kb, and min_fam_span = 20 kb.
  - Step (a), pure refactor: one named prior `S_MIN = 20 kb` ("smallest satellite cluster treated as centromere-sized"), stated
    as an assumption in the README and the GFF header, with the rungs as multiples.
  - Step (b): replace the two rungs that are really relative judgments.
    - `low_satellite`: compare with the same genome's median call.
    - `weak_satellite`: the primary family's best cluster must be enriched over that family's dispersed density on the
      chromosome (one-sided Poisson on monomer counts, conventional p < 0.001).
  - Why S_MIN stays a prior: the only data-derived candidate, scaling with chromosome or genome size, is contradicted here
    (slope -1.6). Functional CENH3-domain size does track genome size in grasses (Zhang & Dawe 2012, Chromosome Res), but
    satellite content is lineage-specific.
  - Evidence that would settle it: CENH3-domain size against satellite bp in species with small chromosomes and satellite
    centromeres (insects, fungi with satellites), from published ChIP.
- **R3. End zone.**
  - Rule: end zone = `subtel_frac_len` × L at each end, with no 1-Mb floor.
  - "Subtelomeric" means the family's bp sits in end zones at **both** ends, on the majority of chromosomes that carry it. A
    centromeric satellite on telocentric or acrocentric chromosomes sits at one end and is then not vetoed. `near_end` stays an
    informational flag.
  - Evidence that would settle the two-ended test: per-family one-end vs two-end occupancy on the 10 complexes. The known
    subtelomeric families (Brassica ~350-bp cluster, subtel_frac 0.86–0.99; TrsC 366) must come out two-ended. If they do not,
    keep the one-sided test and only drop the floor.
- **R4. E-value search space.**
  - Rule: scale each chunk's E-value by genome length, not chromosome length (annotate.py:73). Contigs and chromosomes then get
    the same bitscore cutoff, and evalue_ref_len disappears.
  - At equal E this is n_chrom-fold stricter (5–50×, about +2.3–5.6 bits). Check that *At* CEN180 monomer counts (Naish 66,131)
    and primary-array bp barely move; this is a check, not a target.
  - Interim (zero batch effect): derive evalue_ref_len for borrowed genomes only.
- **R5. rDNA spacer test.**
  - Rule: locate 45S genes genome-wide, for example by adding SSU/LSU profile HMMs to annotate's nhmmer pass, which already scans
    the genome. Measure this genome's rDNA unit as the median spacing between adjacent SSU hits. Call a family an rDNA spacer if
    the majority of its bp lies inside NOR loci (the 45S gene cluster span ± one unit).
  - This removes rdna_flank and rdna_arrays and works for any unit size.
- **R6. Long-unit geometry.** lu_min_unit = TRF_MAX_PERIOD; lu_window_overlap = 2 × lu_max_unit; lu_flank = 2 × lu_max_unit;
  lu_window = engineering constant ≥ 4 × lu_max_unit.
- **R7. Family library size.**
  - Rule: carry forward the largest clusters until they hold ≥ 95% of each genome's seed-array bp, with an engineering ceiling
    that warns when hit. This replaces both min_fam_span and max_families.
  - Unclear: whether 95% adds many families in satellite-rich genomes. Check the 10 complexes: how many families 95% needs, and
    whether any v0.3.7 centromeric family would fall outside.

## 2. Redundancies and coupled parameters

- **One parameter doing several jobs.**
  - chrom_min_bp: 6 jobs (table).
  - cluster_gap: 2 jobs, cluster merging and co-centromeric distance. lu_flank uses the same 500 kb for an unrelated job.
  - strong_prevalence: 2 jobs (`strong` label, dominance share).
  - superfam_ident: 3 jobs (superfamily table, consensus frame, co-centromeric/dominance through `related`).
  - nhmmer_evalue: sampling and annotation, over two different search spaces.
  - The end zone: family eligibility, `near_end` flag and long-unit fallback.
- **One rule split over several parameters.**
  - End zone (subtel_min_bp + subtel_frac_len).
  - Diffuse flag (diffuse_min_len + diffuse_max_frac).
  - rDNA test (rdna_flank, rdna_arrays, rdna_array_frac).
  - Long-unit detector (6 lu_* fields).
- **Same concept, same value: merge into one constant.** These merges do not change results at defaults.
  - 2,000 = max_period = TRF MaxPeriod in TRF_ARGS = lu_min_unit = min_seed_span.
  - 0.90 = fam_ident = chain_ident. 0.97 = dedup_ident = hor_ident. 0.70 = superfam_ident = absorb_ident.
  - 0.90 = hor_min_cov = sampling coverage (families.py:384).
  - 20 kb = min_fam_span = chrom_min_bp. 60 kb = lu_window_overlap = 2 × lu_max_unit.
  - 50 = MAX_SELF_SEQS (pipeline.py:50) = MAX_ALL (plot.py:27).
  - Length tolerance 5% is repeated at families.py:43,132,276 and longunit.py:89. Tolerance 15% is repeated at
    families.py:245,300 and classify.py:57.
- **Coupled pairs that break silently if only one is changed.**
  - trf_chunk and TRF `-l 6` (trf.py:19): a chunk > 6 Mb exceeds the `-l` bound.
  - trf_overlap must be ≥ 2 × MaxPeriod.
  - max_period > 2000 is ignored by TRF.
  - lu_window must be ≥ 2 × lu_max_unit.
  - hor_ident and cd-hit `-n` (hor.py:175, M7).
  - BATCH (families.py:14) changes per-round BLAST database size and so the E-values (M9). An engineering knob changes results.
- **Inconsistent baseline (M4, deferred since 0.3.0).** Gene density is genome-wide for eligibility and the long-unit fallback
  (classify.py:62-64,85; centromere.py:49-51,101), but chromosome-wide for `gene_rich` and the reported gene_ratio
  (centromere.py:82-85,117). On grass chromosomes with steep gene gradients the two disagree. Use the chromosome mean everywhere.
- **Thresholds set with a known answer in view** (from the ledger). Under `issues.md`, a rule must help ≥ 2 independent complexes,
  so each needs a principled anchor recorded instead of the outcome it fixed:
  - dedup_ident: 0.95 was rejected because it merged CentBr1/2.
  - fam_ident: chosen to keep CentBd/CentBs and CentBr1/2 apart.
  - diffuse: fitted to Carabica; only Carabica is affected.
  - dominance rule: B. juncea A09 / B. napus A06.

## 3. Assumptions that break outside large plant genomes

1. **Small chromosomes (≤ 2 Mb: fungi, protists, Drosophila dot).** The end zone max(1 Mb, 5% L) covers the whole chromosome, so
   `end_zone_bp` returns all bp (classify.py:40-44). subtel_frac is then 1 and no family is ever eligible, and every call is
   `near_end`. cluster_gap 500 kb is about half the chromosome; lu_window 1 Mb is larger than the chromosome; diffuse_min_len 5 Mb
   can never trigger. In the batch: 47 *Cryptococcus* chromosomes. Responsible: subtel_min_bp, cluster_gap, lu_*, diffuse_min_len.
2. **Small or regional centromeres** (*Cryptococcus* transposon-rich regional centromeres, Yadav 2018 PNAS; *S. pombe* inverted
   repeats; *Candida* 3–5 kb; *Plasmodium* 2–3 kb AT-rich). The abundance ladder (5k–100k) presumes Mb-scale arrays, and
   non-tandem centromeres are out of scope by design. A species with real ~10-kb satellite centromeres would get `none`.
   Responsible: min_fam_span, chrom_min_bp, fallback_min_bp, low_sat_bp, strong_cluster_bp. State the scope in the README.
3. **Telocentric or acrocentric karyotypes** (mouse; the subtelocentric/acrocentric chromosomes common in cyprinid fish). The
   centromeric satellite sits in one end zone, so it is vetoed as subtelomeric or lowered by `near_end`. Testable on the two fish
   complexes. Responsible: max_subtel_frac, subtel_*. Fix: R3.
4. **Holocentric species** (*Rhynchospora* Tyba arrays spread along chromosomes, Marques 2015 PNAS; *Luzula*, *C. elegans*,
   Lepidoptera). The concentration test fails and the result is `none`, or worse, one arbitrary cluster called "the centromere".
   There is no diagnostic for this. Responsible: min_concentration, cluster_gap. Add a dispersion flag: many evenly spaced arrays
   of one family on most chromosomes.
5. **Contig/scaffold assemblies, unplaced sequences, organelles.** Prevalence is the share of *sequences* (classify.py:82), so
   unplaced scaffolds and chrM/chrC deflate it. MAX_SELF_SEQS only warns. Short contigs get looser nhmmer cutoffs unless the user
   knows to set evalue_ref_len. classify_from must be given by hand. Responsible: strong_prevalence, evalue_ref_len, classify_from,
   contig_end_bp, MAX_SELF_SEQS. Fixes: length-weighted prevalence, R4, auto borrowing.
6. **Many chromosomes** (n > 50: some fish, ferns, butterflies, sturgeons). With MAX_SELF_SEQS = MAX_ALL = 50, a chromosome-level
   genome is treated as a contig assembly in the warning and plotted as 20 tracks. *Carassius*/*Cyprinus* have exactly 50, right
   at the edge.
7. **Short-monomer centromeric satellites** (*Drosophila* AATAT/AAGAG/dodeca, Chang 2019). These are excluded by min_period 30.
   Units of 2–30 kb are found only by the low-confidence long-unit fallback and never become a family; units > 30 kb are missed
   (lu_max_unit).
8. **Very large, retrotransposon-dominated centromeres** (Triticeae/Aveneae, 400–720-Mb chromosomes). Scattered satellite inside
   5–10-Mb CENH3 domains is split by cluster_gap 500 kb into fragments. max_families 40 may truncate satellite-rich libraries.
   Chromosomes without a satellite call are scanned whole by the long-unit detector, about 700 minimap2 windows each (cost, not
   correctness). Responsible: cluster_gap, max_families, lu_*.
9. **rDNA unit size.** rdna_flank 15 kb assumes plant 45S units. With animal-sized units (about 43 kb in human) the genes can lie
   outside the window. R5 fixes this.

## 4. Other portability problems

- **Dataset-specific gene lookup.** The window_identity fallback (genomes.py:53-58, 108-129) reads `window_identity/*/<g>_mod.gff`
  + `run.exit` + `window_identity/genomes/chr_map.tsv` and mRNA lines. Delete it. All 38 batch complexes use `genes/<g>.gff3.gz`;
  check every genome has that file first, otherwise its gene evidence silently turns off.
- **Gene GFF filter is non-standard.** genomes.py:84 filters on `locus_type`, which is neither standard nor written by
  `make_gene_gffs.py`. NCBI/Ensembl `gene` lines also include tRNA, ncRNA and pseudogenes; when `gene_biotype`/`biotype` is
  present, keep only protein_coding.
- **Narrow FASTA discovery.** Only `*.fa`/`*.fa.gz` are found (genomes.py:38); NCBI's `*.fna(.gz)` and `*.fasta` are ignored.
  The output goes inside the input dir (pipeline.py:57), and `--trf-cache` is required (cli.py:26).
- **The rRNA models are not part of the package.** `RRNA_CM` (families.py:15) resolves to the repo-root `data/`, outside the
  package. After `pip install` the CM is missing, and both rDNA tests are skipped with a warning (families.py:157-159, 189-192),
  so results silently depend on the install. Ship it as package data via importlib.resources, and fail fast if it is missing.
- **pyproject.toml is incomplete.** It has no `[build-system]`, no console script (`cenrep`), no package data and no license.
  `environment.yml` names the local env `synLTR`.
- **tools/*.sh are one-off batch helpers.** They contain absolute paths (`/data2/chris/tetraploid/allo10/clean`), a Claude session
  scratchpad path (finish_classify_rerun.sh:8), a hard-coded PID 2535367, 0.3.x version checks, complex names, and calls to a
  launcher (`run_centromere.sh`) that is not in the repo. Remove them from the repo.
- **Dataset-specific docs and dev notes.** `.sdd/` (tracked) and `docs/2026-09-25-*.md` describe allo10. The design doc says
  "exactly 4 per complex" and "fails fast when ≠ 4 genomes", which is no longer true. Move these to dev notes, and write a
  biologist-facing README.
- **Magic constants outside params.py, in no single place:**
  - trf.py: :19 TRF_ARGS, :20 POLL.
  - families.py: :14 BATCH; :63 blastn `-word_size 11 -evalue 1e-5`; :43/47/53 tolerances; :150/:182 cmsearch E 1e-5/1e-10;
    :204 100-bp min flank; :220 trim 20 bp/0.90; :309 gap 0.5; :384 0.9; :430/452 1.05; :446 RNG seed.
  - annotate.py: :56 2-Mb chunk; :96 1.9/0.5; :115 0.2.
  - classify.py: :57 15%.
  - centromere.py: :84 `gene_rich` > 1; :88 confidence ladder.
  - longunit.py: :56 asm20; :89 5%.
  - hor.py: :175-176 cd-hit options; :232 ≥ 4 classed monomers.
  - output.py: :115/125/153 0.5; :128 100-bp consensus.
  - pipeline.py: :50 MAX_SELF_SEQS.
  - plot.py: :24-28 BIN, H, MAX_TRACKS, MAX_ALL, TRACK_IN.
- **Plot bin size is absolute.** plot.py:24 `BIN = 100 kb` gives 8 bins on a 0.73-Mb chromosome and 7,240 on a 724-Mb one
  (bars thinner than a hairline in a 3.6-in panel). Use bin = longest chromosome in the panel / ~300.
- **TRF lock is single-host only.** The lock (trf.py:108-124) checks whether the PID is alive on *this* host. On a shared
  filesystem across SLURM nodes, a live lock held on another node looks stale and is stolen. Store host + PID, and only steal
  locks from the same host.
- **Checkpoint fingerprints are keyed by Params field names** (pipeline.py:96-101). Renaming, merging or moving a field changes
  every fingerprint, and existing checkpoints are refused. A "pure" refactor still forces reruns unless the names are kept or old
  keys are mapped.

## 5. Target CLI and migration order

```
cenrep run GENOME_DIR [-o OUTDIR] [-t THREADS] [--trf-cache DIR] [--classify-from T=S[+S2],...] [--force STAGE] [-v]
```

- GENOME_DIR is today's `--complex-dir`; keep that flag as an alias for the batch driver. OUTDIR defaults to
  `GENOME_DIR/centromere`, and the TRF cache to `OUTDIR/work/trf_cache`. `--classify-from` becomes an optional override of
  auto-detection.
- Everything else moves to one `cenrep/constants.py`, with one line of reason each, and is still written to `params.json` and the
  GFF header. `Params` stays as an internal dataclass, so tests keep `Params(**overrides)`.
- No `--config` file in v1. It would be a back door for tuning; add it only if an outside user needs it.

| step | changes | stages | results change? | STAGE_VERSION bump | rerun of 38 | risk |
|---|---|---|---|---|---|---|
| 0 | constants module + CLI cut (field names kept, so fingerprints are unchanged); same-value merges as aliases; superfam help fix; package data + fail-fast CM; delete tools/*.sh and the window_identity fallback (after the genes/ check); -o, optional cache, .fna/.fasta; host-aware lock | none | no (verify: synthetic run byte-identical apart from the GFF date, as in the 0.4.3 bgzip check) | none | none | low |
| 1 | R3 end zone; chromosome-mean gene baseline (M4); length-weighted prevalence; auto classify_from + derived evalue_ref_len/contig_end_bp (no batch effect); R6 long-unit geometry; relative `low_satellite`; plot bin | classify, centromeres, plot | yes: 347 chromosomes < 20 Mb in 13 non-held-out complexes (Lemna 81, Populus 60, Cryptococcus 54, Actinidia 47, Sphagnum 41, ...) | classify (+centromeres, plot) | `--force classify` (cheap) | medium |
| 2 | R1 density-segment clusters (replaces cluster_gap/diffuse); R2(b) relative rungs | classify, centromeres | yes, widely: ±2× cluster_gap already moves 25% of the 361 calls by > 10% | classify, centromeres | `--force classify` | **highest** |
| 3 | R4 genome-level E-value (drop evalue_ref_len); R5 45S HMMs in annotate | annotate | yes: monomer counts and array bp in every genome | annotate | `--force annotate` (hours on 4–5-Gb genomes) | medium |
| 4 | R7 coverage-based library (min_fam_span/max_families); min_seed_span = TRF_MAX_PERIOD (same value) | families | yes: a new HMM library changes every downstream number | families | `--force families` (BLAST + MAFFT + annotate) | high; do it only if logs show the cap or the 20-kb floor dropping a candidate (check `dropped by max_families` on Aegilops/Avena first) |

Notes:
- Ship steps 1–4 as one release, so the planned final 38-complex rerun starts once, from the earliest changed stage. If evidence
  does not demand steps 3–4, that rerun stays at `--force classify`.
- Validate each step on the 10 allo10 complexes and the reviewed development complexes against `literature_cen_satellites.md`.
  This checks that the rules are sane; it does not tune them. Look at the held-out complexes only after the code is frozen.
- Riskiest, in order:
  1. R1, which moves region boundaries, and the 0.3.7 dominance rulings rest on cluster sizes.
  2. R7, which changes everything downstream.
  3. R3, which changes eligibility on small and telocentric chromosomes.
  4. R4, which changes monomer counts.
