# cenrep — centromeric satellite / HOR annotation for the allo10 complexes

date: 2026-09-25 | status: spec (approved design, pending spec review)

## 1. Goal

For each of the 10 complex dirs in `clean/` (4 genomes each: 2 diploid progenitors +
2 subgenomes), create `<Complex>/centromere/` containing, per genome, a GFF3 of the
centromeric tandem repeats:

- centromeric satellite arrays and their monomers,
- higher-order repeats (HORs) within them,
- other tandem-repeat arrays located inside the called centromere,
- the centromere region call itself.

Scope decisions (user, 2026-09-25):

- installed tools only (no SRF/TRASH/TideCluster);
- tandem repeats only (no LTR-RTs in the GFF);
- satellite-free centromeres are reported as `none`, with no retrotransposon or
  gene-desert placement;
- a long-unit (> 2 kb) tandem detector runs in centromere regions (S5b).

Success criteria

1. 40 GFFs `<Complex>/centromere/<Genome>_<L>_centromeric_repeat.gff`, GFF3-valid
   (every `Parent` exists, 1-based closed coordinates within the sequence length).
2. A. thaliana sanity check: centromeric family ≈ 178 bp (CEN180). About 66k monomers
   genome-wide (Naish 2021: 66,131; nhmmer probe: 66,692). One 2–3 Mb centromere per
   chromosome.
3. Every complex's centromeric families agree with published monomer lengths where
   those exist. Disagreements are explained in the memo, not hidden.
4. Chromosomes with no convincing satellite centromere get `confidence=none` in the
   summary, not a forced call.
5. pytest suite on synthetic genomes passes. Checkpointed stages; a rerun resumes.

## 2. Environment and inputs

- Host: lab server. The pipeline uses at most 32 threads total. It is launched with
  `setsid nohup`, and the launch is checked for PPID 1.
- Env: mamba `synLTR`, which contains:
  - Python 3.10.18: numpy 1.26, scipy, pandas, matplotlib, pyhmmer 0.12.0,
    pyfastx 2.3.0, pywfa, pytest 8.3.4
  - trf 4.10.0-rc.2, mafft 7.525, blastn 2.16.0+, HMMER 3.4, cd-hit-est 4.8.1,
    samtools
- Genomes: `clean/<Complex>/*.fa`, exactly 4 per complex. Headers are
  `<Genome>_<L><NN>`, chromosomes only. `.fai` files are required; Cbursapastoris
  lacks them, so they are created with `samtools faidx`, the only write outside
  `centromere/`.
- Genes (evidence only): miniprot mRNA from
  `clean/<Complex>/window_identity/<run>/<Genome>_<L>_mod.gff`.
  - Only runs with `run.exit` = 0 are used.
  - seqids there are `<Genome>_<L>_chrN`; they are mapped back to the original IDs
    via `window_identity/genomes/chr_map.tsv` (stem, orig_id, chrN).
  - If there is no finished run, the genome has no gene evidence and is flagged
    `no_genes`.
- Shared genomes (md5-identical): Brapa_A (Bjuncea = Bnapus), Bnigra_B
  (Bcarinata = Bjuncea), Boleracea_C (Bcarinata = Bnapus). TRF runs once per md5.
  Everything downstream runs per complex, because families are defined per complex.

## 3. Layout

```
clean/cenrep/                 code (package `cenrep`), tests/, docs/
clean/run_centromere.sh       launcher (all complexes sequentially, or named ones)
clean/centromere_cache/trf/   <md5>.trf.tsv.gz (TRF results, shared across complexes)
clean/<Complex>/centromere/
    <G>_centromeric_repeat.gff        x4   (G = <Genome>_<L>)
    <G>_centromere_summary.tsv        x4   one row per chromosome
    <G>_hor_units.tsv                 x4   one row per HOR unit
    satellite_families.tsv                 family metrics + class per genome + reason
    satellite_families.fa                  family consensus monomers
    satellite_families.hmm                 profile HMMs (HMMER3)
    centromere_overview.pdf                figure (panels A-D = 4 genomes)
    cenrep.log  cenrep.exit  cenrep.time   run log / exit code / GNU time
    work/                                  per-stage intermediates + .done markers
```

## 4. Pipeline (per complex; stages are checkpointed with `work/<stage>.done`)

### S1 trf — tandem-repeat discovery (per genome, cached by md5)

- Chromosomes are split into 5 Mb chunks with 10 kb overlap. They are streamed with
  pyfastx, and one chunk at a time is written to a temp FASTA.
- Each chunk runs `trf chunk.fa 2 7 7 80 10 50 2000 -l 6 -h -ngs`. The chunks run in
  a process pool (≤ 32 workers across the 4 genomes).
- Parameter rationale:
  - Benson's recommended weights 2/7/7, PM 80, PI 10.
  - Minscore 50 keeps weak or short repeats; they are filtered downstream.
  - MaxPeriod 2000 is the TRF maximum. It is needed for long monomers (e.g. the
    ~728 bp Erythranthe satellite) and lets TRF report HOR-scale periods.
  - `-l 6`: chunks are ≤ 5 Mb, so no array can be longer.
- Parsed into `start end period copies consensus_size pct_match pct_indel score
  entropy consensus`. Coordinates are shifted to the chromosome and rows are
  de-duplicated in overlaps (same period ±1 and span overlap ≥ 50% → keep the higher
  score). The flank and repeat_seq columns are dropped.
- Output: `centromere_cache/trf/<md5>.trf.tsv.gz`, then symlinked into `work/trf/`.
- Measured: 4 Mb of A. thaliana CEN1 takes 58 s and 2.5 GB RSS (at `-l 10`).
  Expect ~1–2 h wall time for the 37 unique genomes at 32 workers.

### S2 families — satellite family library (per complex)

1. **Seeds.** TRF rows from all 4 genomes with 30 ≤ period ≤ 2000 and span ≥ 2 kb.
   Each seed has a consensus, a period and a span.
2. **Rotation-invariant, multimer-aware similarity** (blastn 2.16):
   `-task blastn -word_size 11 -evalue 1e-5 -dust no -soft_masking false
   -max_target_seqs 100000`. The query is the seed monomer; the subject database is
   every seed's consensus written twice in a row (a dimer). Short seed S matches long
   seed L at identity threshold T (period ratio k = pL/pS within 5% of an integer
   k ≥ 1) if:
   - there is an HSP with pident ≥ T covering ≥ 80% of S, and
   - S's HSPs tile ≥ 80% of L's dimer.

   Both strands are allowed.
3. **Greedy centroid clustering (no chaining).** Seeds are sorted by per-complex span,
   largest first. Each seed joins the first existing family whose representative it
   matches as above; otherwise it founds a new family. The representative is the
   founding seed. If a later member's period is k-fold smaller and it matches, the
   family period becomes the base period.
   - **Family threshold T = 90%.** Seeds are TRF array consensuses, which are much
     more homogeneous than single monomers. At 90%, published subfamilies that must
     stay separate for allopolyploid comparisons are kept apart:
     - CentBd vs CentBs, ~87% identical (19 SNPs + 1 indel)
     - CentBr1 vs CentBr2, ≤ 86.9%
     - CEN180 vs pAa, 58–80%
     - pAge1 / pAge2 / pAa (~30% divergent)
     - CentO vs CentO-C
   - **Superfamily (descriptive only):** the same greedy rule at T = 70% over the
     family representatives, reported in `satellite_families.tsv`. For example,
     CentBd and CentBs fall in one superfamily.
4. **Families carried forward:** union span ≥ 20 kb in at least one genome, capped at
   the top 40 families by max per-genome span. Every other seed is kept only as
   evidence for the S7 `tandem_repeat` loci.
5. **Consensus and HMM, two passes:**
   - (a) Build a single-sequence HMM from the representative, rotated so that it
     starts at its lexicographically smallest 16-mer.
   - (b) Run nhmmer only on that family's TRF seed intervals, and sample up to 300
     near-full-length hits (≥ 90% model coverage), taking them from the largest
     arrays of all 4 genomes in proportion.
   - (c) Align them with `mafft --auto --thread N`. The consensus is the column
     majority (columns with ≥ 50% gaps are dropped), rotated by the same smallest-16-mer
     rule.
   - (d) Build the profile HMM with the pyhmmer Builder from that MSA.
6. **Family ID:** a provisional `famNNN` in S2–S6. S7 assigns the final
   `<Complex>_sat<NN>_<consensus_len>`, where NN ranks families by total annotated bp
   over the complex's 4 genomes (from S3). `satellite_families.tsv` keeps both IDs.

### S3 annotate — genome-wide monomer calls (per genome)

- pyhmmer nhmmer of all family HMMs against one chromosome at a time (streamed, never
  the whole genome in memory), with `cpus=32` and E ≤ 1e-5.
  - Measured: CEN180 single-sequence HMM vs the whole A. thaliana genome took 62 s on
    8 CPUs and gave 66,692 hits.
- Any two hits (same or different family) overlapping by > 50% of the shorter are
  resolved in favour of the higher bitscore.
- Per monomer: family, strand, bitscore, model coverage, and identity to the
  consensus (exact matches ÷ alignment columns, from the nhmmer alignment).
- **Arrays:** consecutive monomers of one family on one strand with gap ≤ 2 × monomer
  length form an array. Arrays need ≥ 3 monomers. A strand switch or a longer gap
  (e.g. a TE insertion) splits the array. Monomers not in arrays are dropped.
- Output: `work/annotate/<G>.monomers.tsv.gz` and `<G>.arrays.tsv`.

### S4 classify — which families are centromeric (per genome)

Per family, over the genome's arrays:

| metric | definition |
|---|---|
| `bp` | total array span |
| `prevalence` | fraction of chromosomes with ≥ 20 kb of the family |
| `concentration` | median over those chromosomes of (bp in the largest cluster ÷ bp on the chromosome). Clusters are arrays merged with gaps ≤ 500 kb. |
| `subtel_frac` | fraction of `bp` within max(1 Mb, 5% of chromosome length) of either end |
| `gene_ratio` | gene (mRNA) density inside the main clusters ÷ chromosome-wide density (NA if `no_genes`) |

Rules. The call is made per chromosome, because published centromeric satellites can
be chromosome-specific: A. lyrata AlyCEN168 dominates CEN2/5/6/8 and AlyCEN179
dominates CEN1/3/4/7, A. halleri pAge1/pAge2/pAa vary by chromosome and strain, and
Brassica C03/C04/C08/C09 carry only CentBr2. A genome-wide prevalence cut would drop
these.

- **eligible** (family level): bp ≥ 20 kb, concentration ≥ 0.5, subtel_frac < 0.5,
  and gene_ratio < 0.5 (the gene test is skipped when NA). Prevalence is reported but
  not required.
- **chromosome's centromeric family:** the eligible family with the largest main
  cluster (≥ 20 kb) on that chromosome.
- **co-centromeric:** another eligible family whose main cluster lies within 500 kb
  of the chosen cluster. This covers layered centromeres such as a CentBr1 core
  flanked by CentBr2.
- **centromeric (primary):** the union, over all chromosomes, of chosen and
  co-centromeric families. `satellite_families.tsv` lists, per family, the
  chromosomes it defines.
- **Evidence strength per family:** `strong` if prevalence ≥ 0.5 and the median main
  cluster is ≥ 100 kb, else `weak`. Weak primaries are reported in the summary and
  checked by hand against the literature.
- **cen_associated:** a non-primary family with ≥ 50% of its bp inside centromere
  regions (from S5).
- **other:** everything else.

Expected true negatives (literature, 2026-09-25 research), where there must be no
forced satellite call:

- B. nigra B genome and the B subgenomes of B. juncea / B. carinata: ALE/CRM/CRB
  centromeres flanked by pBN35 islands; only a 28-copy CEN176 relic.
- All Coffea: CRC-dense, no short satellite.
- O. officinalis CEN2, and probably O. minuta C02.
- B. distachyon Bd5: only ~40 kb of CentBd.

These chromosomes are expected to come out `low`/`none` or with a `weak` primary. The
memo must say so explicitly.

Every metric and the reason for each call go into `satellite_families.tsv`. The
thresholds are CLI parameters.

### S5 centromeres — region call (per chromosome)

- Primary-family arrays merged with gaps ≤ 500 kb form clusters. The centromere is the
  cluster with the most bp; its span runs from the first array start to the last
  array end.
- Flags:
  - `low_satellite`: primary bp < 50 kb
  - `competing_cluster`: the second cluster has ≥ 0.5 × the main cluster's bp
  - `near_end`: the region lies within max(1 Mb, 5% of chromosome length) of an end
  - `gaps=n`: runs of ≥ 10 N inside the region, with their total bp
  - `gene_rich`: gene density inside the region > chromosome mean
- Confidence: `high` = none of low_satellite / competing_cluster / near_end; `medium`
  = one of them; `low` = two or more.
- Fallback (tandem evidence only):
  1. If there is no primary cluster of ≥ 20 kb, the largest primary-family cluster of
     ≥ 5 kb gives the region at `confidence=low`, with `evidence=weak_satellite`.
     cen_associated families can't be used here: they are defined only after the
     regions exist.
  2. Otherwise, the S5b long-unit arrays found on that chromosome are tested with the
     S4 eligibility rules (concentration, subtel_frac, gene_ratio). If they pass, the
     region is `confidence=low`, `evidence=long_unit_tandem`.
  3. If neither applies, the chromosome is `none`: a summary row, no GFF `centromere`
     feature, and all of that chromosome's primary-family arrays with
     `in_centromere=false`.
- **No retrotransposon or gene-desert placement** (user, 2026-09-25: satellite-free
  centromeres → `none`). Gene density is only used as a veto in the eligibility test.

### S5b longunit — tandem arrays with units > 2 kb (per genome)

- **Why:** TRF's MaxPeriod is 2000. Coffea heterochromatin has CRM-derived tandem
  arrays with ~2.7 kb units (Scalabrin 2024: 690 kb on Chr7c).
- **Where** (user, 2026-09-25: "in cen regions"):
  - every centromere region from S5, ±500 kb;
  - plus, for chromosomes S5 leaves without a satellite-based call, the whole
    chromosome. Without this, the detector would never run on Coffea. Its arrays can
    only rescue a call through S5 fallback step 2.
- **Method:** each target is cut into 1 Mb windows with 60 kb overlap. Each window is
  aligned to itself with `minimap2 -x asm20 -c -DP --dual=no` (`-DP` drops the
  trivial self-diagonal).
  - Keep alignments whose shift d = |t_start − q_start| is between 2,000 and 30,000 bp,
    whose length is ≥ d (at least two adjacent copies), identity ≥ 80%, and whose two
    copies are on the same strand.
  - Each alignment implies an array covering [min start, max end] with unit ≈ d.
  - Arrays that overlap and whose units agree within 5% are merged. The array's
    unit is the median d; copies = span / unit.
- **Output:** `tandem_repeat` features with `period=<unit>`, `unit_source=selfaln`,
  copies and mean identity. The ≥ 500 bp, inside-centromere and not-already-covered
  rules are the same as for TRF loci.

### S6 hor — higher-order repeats (per genome, per centromere)

1. **Monomers:** primary-family monomers with model coverage ≥ 90%, taken from the
   whole genome. Minus-strand monomers are reverse-complemented so all are in
   consensus orientation.
2. **Variant classes:** `cd-hit-est -c 0.97 -n 10 -G 1 -aS 0.9 -aL 0.9 -g 1 -r 0
   -d 0 -M 0 -T 32` on all of them. 97% identity is ≈ ≤ 5 differences per 178 bp,
   the HOR identity rule of Naish 2021. Classes are genome-wide, so a class shared
   between chromosomes is visible.
3. **Class strings:** each centromere becomes the string of its arrays' class IDs,
   read 5'→3' in consensus orientation (minus-strand arrays in reverse genomic order),
   with a unique separator between arrays. Partial and filtered-out monomers become
   unique separators too.
4. **Maximal repeats:** a suffix array plus LCP (Kasai) on the integer string, then
   LCP-interval enumeration. A maximal repeat is a substring of ≥ 2 monomers that
   occurs ≥ 2 times and is left- and right-diverse.
5. **Valid HOR units:** a unit must have ≥ 2 distinct classes, and must not contain
   two full copies of a shorter period (smallest period p must satisfy 2p > length).
   So `abcdabcd` is rejected, but `aba` and `abcdabc` are allowed. Runs of one class
   are homogeneous monomers, not a HOR.
6. **Greedy assignment:** units are ranked by achievable coverage, then length. The
   achievable coverage is the number of non-overlapping occurrences (taken left to
   right) × length. An occurrence is accepted only if none of its monomers is already
   assigned, and a unit needs ≥ 2 accepted occurrences.
   - Why coverage, not "longest first": in a tandem `abcd`×14 array, the longest valid
     repeat is `abcdabc` (1.75 copies). Picking it first would leave every 8th monomer
     orphaned and hide the 4-monomer tandem HOR.
   - Unit lengths are capped at 200 monomers (`--hor-max-len`), which bounds the
     runtime. Longer duplications are reported in 200-monomer pieces.
7. **Tandem HOR runs:** accepted occurrences of one unit that are back-to-back are
   merged, recording `tandem_copies`.
8. **Output:** per-monomer `variant_class` and `hor` attributes, HOR occurrence
   features, `<G>_hor_units.tsv` (unit, n_monomers, occurrences, class string, max
   tandem copies, chromosomes), and per-array/centromere `hor_fraction`.

### S7 output — GFF3 + tables (per genome)

`##gff-version 3`, `##sequence-region` for each chromosome, and a comment header with
the cenrep version, date and parameters. Source column = `cenrep`. Sorted by seqid,
then start, with parents before children.

| type | parent | key attributes |
|---|---|---|
| `centromere` | — | ID=cen_<chr>; confidence; flags; main_family; satellite_bp; sat_fraction; n_arrays; hor_fraction; n_genes; gaps |
| `satellite_DNA` | centromere if inside, else none | ID=arr_...; family; family_class=centromeric\|cen_associated\|other; monomer_len; n_monomers; mean_identity; hor_fraction; in_centromere=true\|false |
| `satellite_monomer` | its array | ID=mon_...; family; identity; model_cov; variant_class; hor (unit ID or .) — score column = nhmmer bitscore |
| `higher_order_repeat` | its array | ID=hor_...; unit; unit_monomers; unit_occurrences; tandem_copies |
| `tandem_repeat` | centromere | ID=tr_...; period; copies; pct_match; consensus (≤ 100 bp, else truncated + consensus_len) |

Which arrays go in:

- Primary-family arrays: all of them, genome-wide, with `in_centromere` set.
- cen_associated and other families: only arrays overlapping a centromere.
- `tandem_repeat`: TRF rows merged into loci wherever they overlap. A locus is
  reported if it is ≥ 500 bp, lies inside a centromere, and overlaps no
  `satellite_DNA` feature by ≥ 50%. Its period and consensus come from the locus's
  top-scoring row.

Other tables:

- `<G>_centromere_summary.tsv` columns: chrom, chrom_len, cen_start, cen_end, cen_len,
  confidence, flags, main_family, main_family_bp, cen_sat_bp, sat_fraction, n_arrays,
  largest_array_bp, n_monomers, mean_identity, hor_monomers, hor_fraction,
  top_tandem_hor, n_gaps, gap_bp, n_genes, gene_ratio, other_tandem_bp.
- `satellite_families.tsv`: one row per family × genome (metrics + class + reason).

### S8 plot — `centromere_overview.pdf` (per complex)

- Panels A–D, one per genome, labelled in the upper-left, with no title.
- Each chromosome is one horizontal track. The track shows:
  - bars of primary-family bp per 100 kb bin, one Okabe–Ito colour per family;
  - bars for cen_associated bp;
  - a line for gene density;
  - an outlined box for the centromere region.
- Vector PDF, 300 dpi for any raster elements, fonts sized for a 2-column journal page.

## 5. CLI

```
python -m cenrep run --complex-dir clean/<Complex> --threads 32 \
    --trf-cache clean/centromere_cache/trf [--force STAGE] [-v]
```

- The default log is minimal (start, stage milestones, done, errors); `-v` adds
  per-chunk and per-family progress plus sanity checks.
- Fails fast when there are ≠ 4 genomes, a tool is missing, a FASTA is corrupt or
  `.fai` is missing (Cbursapastoris is handled by the launcher), or TRF/nhmmer exits
  non-zero.
- Warns and skips on: an unparsable TRF line, a family with < 20 sampled monomers (it
  keeps the single-sequence HMM and is flagged `hmm=single`), or a missing gene GFF.
- `run_centromere.sh [Complex...]` runs complexes sequentially, each at 32 threads, in
  one detached `setsid nohup` job. It writes `cenrep.{log,exit,time}` per complex
  (GNU time -v).

## 6. Testing (pytest, synthetic fixtures; no real data)

`tests/conftest.py` builds a random 3-chromosome genome (~300 kb each, seeded RNG)
with:

- a 150 bp centromeric satellite array (2% mutations) in the middle of each chromosome;
- HOR structure: a 4-monomer unit repeated in tandem 15× on chromosome 1, and a
  3-monomer unit duplicated non-tandemly on chromosome 2. It has 3 monomers because the
  consensus rotation can shift monomer boundaries by a partial monomer, which leaves
  ≥ 2 consistent monomers per copy;
- a 90 bp subtelomeric satellite decoy at both ends of every chromosome;
- a 600 bp (AT)n microsatellite and a 500 bp N gap inside chromosome 3's centromere;
- a minus-strand block of the centromeric satellite;
- on chromosome 3, a 3 kb-unit tandem array (6 copies, 3% mutations) inside the
  centromere, for S5b;
- a 4th, satellite-free chromosome, which must come out `none`.

Unit tests:

- TRF `-ngs` parsing and chunk offset/overlap de-duplication;
- the rotation/multimer similarity rule (rotated copy → same family; 3× multimer →
  same family with base period; unrelated → new family);
- array chaining (gap/strand splits);
- the suffix array, LCP and maximal repeats against a brute-force oracle on random
  small strings;
- HOR validity rules (non-primitive, single-class) and greedy assignment;
- GFF3 validity (Parents exist, coordinates in range, parents contain children);
- classify rules (decoy → `other`, planted satellite → primary).

End-to-end tests on the synthetic genome (skipped if trf/nhmmer are missing):

- centromere regions recovered within ±1 monomer of the planted boundaries;
- the tandem HOR found with period 4 and ≥ 14 copies;
- the chromosome 2 duplication found (unit ≥ 2 monomers, 2 non-adjacent occurrences);
- the (AT)n locus reported as `tandem_repeat` with period 2;
- the decoy classified `other`;
- the gap flagged;
- the 3 kb-unit array reported with period within 5% of 3,000;
- the satellite-free chromosome `none`.

## 7. Validation on real data (before calling it done)

1. A. thaliana (Asuecica/Athaliana_T): family ≈ 178 bp, ~66k monomers, 5 high-confidence
   centromeres of 2–3 Mb.
2. Cross-check every complex's primary families against the literature (research
   2026-09-25). The table is saved as `cenrep/docs/literature_cen_satellites.md`.

   | complex | expected |
   |---|---|
   | Asuecica | CEN180 ~178 bp (T) + AthCEN159 minor; pAa ~179 bp (A) |
   | Akamchatica | AlyCEN168 (pAge-type) + AlyCEN179 (pAa-type), chromosome-specific; A. halleri pAa/pAge1/pAge2 ~180 bp |
   | Cbursapastoris | C. rubella satellite 168 bp; nothing published for C. grandiflora, C. orientalis, C. bursa-pastoris |
   | Brassica A/C | CentBr1/CentBr2 176 bp; C03/C04/C08/C09 CentBr2-only; A03/A05 CentBr2 |
   | Brassica B | no satellite (see S4 negatives) |
   | Bhybridum | CentBd 156 bp (D), CentBs 157 bp (S); Bd5 ~40 kb only |
   | Ominuta | CentO 165 + 155 bp (B, O. punctata); CentO-C 126 bp (C); TrsC/CentO-C2 366 bp is subtelomeric too, so not a marker; O. officinalis CEN7 carries CentO and CEN2 no satellite |
   | Msookensis | Cent728 (728 bp), all chromosomes; 27–41 Mb per genome |
   | Carabica | no short satellite; > 2 kb CRM-derived units (e.g. 2,683 bp, 690 kb on Chr7c) |

   A. thaliana CEN178 template (TRASH README) for a direct identity check:
   `AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG`
3. Subgenome vs progenitor consistency: the same family class and similar centromere
   positions on homoeologous chromosomes. Differences are reported, not forced to
   agree.
4. Look at every `centromere_overview.pdf` and each low/none confidence chromosome.
5. Shared diploids: the same TRF cache is used in both complexes, and the centromere
   calls agree (they can differ only through family IDs and complex-specific
   families).

## 8. Out of scope

LTR-RTs (CRM/Athila) in the GFF or as positional evidence; gene-desert-only
centromere placement; CENH3 or other epigenetic evidence; cross-complex
family unification; naming families after known satellites automatically (done by
hand in the memo from the literature check).
