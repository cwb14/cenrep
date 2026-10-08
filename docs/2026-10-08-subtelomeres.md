# Subtelomeres and rDNA tags across genomes (design, 2026-10-08; issues I002, I016, I018)

## Subtelomeric families
- **v0.8.2 and earlier:** end zone = max(1 Mb, 5% of length) at each end, one-ended.
  - Chromosomes ≤ 2 Mb were all end zone (*Cryptococcus*).
  - The CENH3-bound Type II centromere of *Fragaria × ananassa* Chr1D (Jin 2025), 1.0 Mb from the end, was
    vetoed. The chromosome had no call.
- **v0.8.3:** 5% end zone, no floor; bp counted only where a family sits at both ends of a chromosome.
  - Fixed *Fragaria* Chr1D.
  - Let one-sided subtelomeric repeats through: rice TrsC on *O. malampuzhaensis* chr06B, and *Solanum* 182-bp
    blocks at chromosome starts (41/72 *Solanum* calls changed in a replay).
- **v0.8.4:** the textbook subtelomere, the stretch between a chromosome end and its outermost gene.
  - A family is subtelomeric if most of its bp lies beyond the outermost gene loci (`terminal_bp`).
  - Fallback with < 2 gene loci on a sequence: the 5% end zone.
  - A centromere has a gene-bearing arm on both sides even next to an end. Subtelomeric repeats have none
    beyond them, at one end or both.
  - Checks:
    - *Fragaria* Chr1D Type II array: 5 genes beyond it, kept.
    - TrsC: first gene at the array's end, rejected.
    - *Solanum* 182-bp blocks: 0 genes beyond, rejected.
  - Replay (classify onward) against the batch outputs:
    - *O. malampuzhaensis* 1/48 calls changed (chr06B: TrsC → CentO inside the CRM peak).
    - *Fragaria* 2/84 (Chr1D gains its Type II centromere).
    - *Lemna* 2/126; *Solanum* 6/72 (was 41 with v0.8.3); *Poa* 0/28.
  - Residual (I020): small relatives of a subtelomeric family can still be called near an end when gene models
    reach into the subtelomere (4 *S. acaule* SG3 calls, all near_end).
- The `near_end` flag still marks calls within 5% of an end.
- **v0.8.5: one verdict per family across the complex (vote).** A family is subtelomeric in every genome when it is
  subtelomeric (`terminal_bp` share ≥ 0.5) in at least half of the chromosome-level genomes that carry ≥ 20 kb of it;
  otherwise in none. One family is one sequence, and a truncated assembly end can make a centromeric array look
  terminal in one genome.
  - *Populus* clone 741 chr14A/14D/09D: the arrays end the assembled chromosome (0 genes beyond, vs 19–28 in the other
    assemblies), so the Palv78/107/156 families were subtelomeric there alone. Restored by the vote.
  - TrsC (*Oryza*) and the *Solanum* 182-bp block (6/6 genomes) stay rejected.
  - Replay (classify onward, with I021 and I022): *Populus* 13/95 calls changed vs v0.8.4 (true families labelled),
    *O. punctata* 0/48, *O. malampuzhaensis* 1/48 and *Fragaria* 2/84 (both as v0.8.4), *Lemna* 2/126, *Poa* 0/28,
    *Solanum* 16/72 vs v0.8.2 (CENH3 within 1 Mb unchanged at 7/11).
  - Rejected alternative: pooling bp over the complex lets one large genome decide for all.

## rDNA spacer tags
- A family tagged as a 45S rDNA spacer repeat in one genome (most of its largest arrays within 15 kb of a 45S gene)
  is excluded in every genome of the complex: one family is one sequence.
- rRNA genes are often missing from an assembly, so their absence next to the arrays is weak evidence.
- *Lemna*: Lj7182 chr19T was called on the 45S spacer repeat (1/6 arrays flanked in that genome; 5/9 and 6/10 in two
  others). It now has no call.

## Rejected alternative
- Pooling the two-ended share over the complex fixes TrsC (two-ended in the C genomes). It does not fix
  one-sided subtelomeric repeats (*Solanum*), so the definition itself had to change.

Versions: classify 0.9.0 (v0.8.5), centromeres 0.9.0.
