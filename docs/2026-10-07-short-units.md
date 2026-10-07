# Short satellite units: significance of runs of copies (design, 2026-10-07)

## Problem (issue I009)
- Families with units shorter than ~60 bp were invisible to the annotation: nhmmer reported no significant hit.
  - *Avena magna* × *longiglumis*: 9 of 40 families (38–51 bp) had 0 bp annotated in every genome. Their TRF seeds
    span up to 159 Mb in one genome.
  - Same pattern in the *Aegilops* 44–46-bp families and in *Avena agadiriana* (41 bp). Over 5 complexes,
    families ≤ 45 bp recovered a median 0% of their own seed arrays, 45–60 bp 5%, ≥ 60 bp ~100%.
  - The same failure hit model building: monomer sampling found < 20 significant hits, so these families got
    single-sequence models ("only 0 monomers sampled").
- Cause: HMMER entropy-weights DNA models to ~0.6 bits per position. A 39-bp model scores at most ~17–19 bits.
  E ≤ 1e-5 over a 500-Mb chromosome needs ~47 bits, so no single copy can ever pass.
- Why it matters: short centromeric satellites exist, e.g. the 41–42-bp chicken CNM repeat. Undetectable families
  also waste slots of the 40-family cap (I003).

## Rule: a run of copies is significant even when no copy is
- The evidence for a satellite is the tandem array, not one copy. After one hit, a random hit starts within the
  array gap W with probability ~W × E_i / L:
  - W = 2 units + 1, the same gap that joins monomers into arrays (`array_gap_factor`).
  - E_i is the hit's E-value for the chromosome; L is the chromosome length.
- A run of r whole copies (hits over ≥ 90% of the model, the existing whole-monomer rule; same family and strand;
  gaps ≤ 2 units) therefore has E = ∏ E_i × (W / L)^(r−1).
- A hit is kept if its own E or its run's E is ≤ 1e-5, the existing threshold.
- Fragments count only on their own E:
  - An array is made of whole copies.
  - Weak partial hits of long models cluster in repeat-rich DNA, so they are not independent, and a 2-kb unit
    allows a 3.9-kb window.
  - With fragments allowed, *A. thaliana* chr1 gained ~100-bp pieces: CEN180-like family 340 → 593 monomers for
    +13% bp; 1.5–1.9-kb families 0 → 63 hits.
- Weak evidence never overrides strong evidence. The multimer rule (a long hit covered ≥ 50% by hits of
  ≥ 1.9×-shorter models is redescribed by the short unit) counts only short copies that are significant on their own.
  - Run-only copies fill sequence that nothing else explains. Alone-significant hits always outscore them, so they
    also lose every overlap contest.
  - The new annotation is therefore the old one plus additions.
  - Tried first: counting run-only copies, then coverage per short unit, then ≥ 2 copies. On a 10-Mb *Avena*
    window these erased the 87-bp and 130-bp composite arrays (87 bp = a 48-bp plus a 39-bp unit). The weak
    short copies then lost the overlap contest, and per-family runs broke, so annotated bp in arrays fell 10–31% in
    dense 300-kb stretches.
- Hits are kept up to HMMER's default reporting E of 10 per search (pyhmmer's TopHits also holds weaker ones);
  significance is decided afterwards.
- Model sampling uses the same rule, so short units get profile models.
- No new threshold. The rule is uniform for all families: long units, which are significant alone, are unchanged.
- Chunk-overlap duplicates are removed before runs are scored, so no copy counts twice.

## Checks (not tuned)
- Own-array recall of the batch models on 100 kb of each family's largest seed array, at chromosome-scale E,
  with nhmmer time. "runs" counts runs of any hits; with whole copies only, recall is 0.70–0.97, except fam004 at
  0.00 (its 87-bp-period arrays fit the 39-bp model only partially; they were 0 before too):
  | family | unit | per hit | runs | time / 100 kb |
  |---|---|---|---|---|
  | *Avena* fam004 (arrays of period 87) | 39 | 0.00 | 0.66 | 0.1 s |
  | *Avena* fam015 | 40 | 0.00 | 0.92 | 0.2 s |
  | *Avena* fam018 | 39 | 0.62 | 0.98 | 0.1 s |
  | *Avena* fam020 | 40 | 0.00 | 0.91 | 0.2 s |
  | *Avena* fam040 | 41 | 0.51 | 0.77 | 0.1 s |
  | *Ae. geniculata* fam013 / fam026 | 46 / 44 | 0.00 / 0.93 | 0.98 / 0.98 | 0.3 / 0.9 s |
  | long controls: *Avena* fam037, fam011; *A. thaliana* fam001 | 87 / 346 / 176 | 0.98 / 1.00 / 0.85 | unchanged | 0.6–2.6 s |
- False positives: every model above against 1 Mb of gene-rich *A. thaliana* chr1 (5–6 Mb) and against the same
  sequence shuffled. 0 hits accepted (0–6 weak hits reported per search). The synthetic tests do the same on
  random sequence.

## Rejected alternative: tandem-copy models (tried first)
- Each unit < 150 bp was modelled as k copies (k = ceil(150 / unit)) and hits were cut back into units.
- Recall was similar (0.89–0.99 at 3–4 copies).
- But nhmmer slows sharply on arrays with multi-copy models, because the model aligns equally well at every unit
  shift. 39-bp unit on 100 kb of its array: 1 copy 0.1 s, 2 copies 2.8 s, 4 copies 43 s. Other short families:
  2–80 s per 100 kb. Long models take 0.4–1.9 s.
- With hundreds of Mb of short-unit arrays per *Avena* C genome, that is many CPU-hours. Rebuilding one family
  model took 83 min.
- The run rule reaches the same recall at the speed of a single-copy model.

## Other changes
- Per-family random seed for monomer sampling (`random.Random(i)`): a family's model no longer depends on how
  many hits earlier families drew.
- Stage versions: families 0.5.0, annotate 0.5.0. `array_gap_factor` moves to the families stage fingerprint
  (used by sampling). Every complex reruns from families in the final sweep; TRF caches stay valid.
