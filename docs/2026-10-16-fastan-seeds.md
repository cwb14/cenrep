# FasTAN seed engine (design, 2026-10-16)

Status: current in v0.9.0 (default engine; `--seed-engine trf` keeps TRF).

## Problem
- Stage 1 (TRF on 5-Mb chunks) was the least predictable cost of a run. TRF's time explodes on some satellite-rich
  chunks: *Mimulus sookensis* spent 2.8 h in it (10,014 s at 16 threads), and one 1-Mb piece of a *Carassius* chunk
  took 196 min alone. FasTAN (Myers) did the same 5-Mb *Carassius* chunk in 1.6 s.
- FasTAN reports each tandem array's interval, mean unit and % differences between adjacent copies, and (with
  `-p`) the cuts between copies, but no unit sequence. Seeds need one: the families stage clusters seed consensuses
  at 90% identity (BLAST, monomer vs dimer) and founds each model on a whole-unit consensus.

## Engine
- Same chunking, cache, locks and chunk-overlap dedupe as TRF (`seeds.py`); each engine is a module (`fastan.py`,
  `trf.py`) with its tools, cache key, duplicate tolerance and `chunk_rows`.
- Per chunk: `FAtoGDB` (FasTAN masks only from a GDB), `FasTAN -m -p -T1`, `ONEview -h` (text .1ano).
  FasTAN's search is fixed (no options): units 2 bp – 30 kb, ≥ 1.95 copies, adjacent copies aligned at ≥ ~70%.
- **Units > 2 kb are dropped**, as TRF's maximum period was 2 kb: the long-unit detector (minimap2 self-alignment)
  keeps that range, and masks the seed rows first.
- **Consensus**: abPOA of up to 30 copies spread along the array, most-frequent-base mode (gap-majority columns
  dropped: the rule of `families.msa_consensus`). 30 copies make a column's majority call wrong with p < 1e-3
  even at 20% divergence (binomial tail). Copies off the unit by > 15% (`RELATED_LEN_TOL`) are misparsed (a missed
  or extra cut) and left out. One usable copy is used as is; none (near tandems, two copies with a gap): the first
  unit from the genome.
  - Checked against TRF's own monomer consensus of the same *A. arenosa* arrays (516 seeds): mean identity
    0.995 (median 1.000). One raw copy instead: 0.957, min 0.849, which would break the 90% clustering.
  - Synthetic copies at 10–20% divergence with indels: consensus 99.4–100% identical to the true unit.
- **The FasTAN 0.8 release ends every array one copy early.** It writes the end of the lagging copy of the
  first-wave alignment (`bepos`); upstream fixed it to the array end (`aepos`) on 2026-08-09, after the release that
  bioconda packages. cenrep probes the installed binary at stage start (ten copies of a 150-bp unit) and, for an old
  build, adds one mean unit to each end, stopping at an assembly gap. On a 5-Mb *A. arenosa* chunk the corrected
  ends match the fixed build's: median error 0 bp, mean 1.25 bp, max 30 bp (1,108 arrays).
- Seed-table columns as TRF's: `pct_match` = 100 − FasTAN's % differences, `score` = matched bp (only ranks
  overlapping rows), `pct_indel` and `entropy` = NaN (not measured). Chunk-overlap duplicates allow a 5% period
  difference (`LEN_TOL`): a FasTAN unit is a mean over the copies aligned.
- The TRF cache key is unchanged, so TRF caches of earlier versions stay valid, including the old default
  `OUTDIR/work/trf_cache` (used when it exists); FasTAN caches are `*.fastan.tsv.gz`.
  Stage `trf` is renamed `seeds`; checkpoints of either engine are never reused by the other.

## Behaviour that differs from TRF
- FasTAN reports the monomer where TRF often reports a multimer: *A. arenosa* pAa seeds at 179 bp 4.61 Mb (TRF
  4.74 Mb), but TRF also reports 2.1 Mb at 358 bp and 3.3 Mb at 537 bp.
- FasTAN breaks long arrays into pieces of ~10–80 kb that abut (0–3 bp apart once the last copy is added); seeds
  are intervals, so this costs nothing downstream, but a non-family array can give a few abutting `tandem_repeat`
  features. Its shortest unit is 2 bp, so homopolymer runs come out at 2–6-bp periods.
- Missed: arrays below ~70% adjacent-copy identity. *A. arenosa*: 4 of 768 TRF seeds (11.6 kb), all 71–74%.
- Fewer short 2-copy repeats (26,357 vs 44,589 rows in *A. arenosa*). The long-unit detector N-masks seed rows, so
  it sees more: *Cryptococcus* polyploid D chr09D gets a long-unit call (6.3-kb unit, 3.9 copies, 99.9%) that TRF's
  2-copy sub-repeats inside each unit (66/108/252 bp) had masked.

## Checks
- Synthetic complex (test suite): every end-to-end truth holds for both engines (calls within 150 bp of the planted
  centromeres, decoy, HOR, duplication, gap, microsatellite, 3-kb long unit, resume, bgzip).
- Real genomes: see the benchmark below.

## Benchmark against TRF (v0.8.9 batch outputs, not rerun)
Five smallest complexes of the batch, same inputs (genes/, tes/) and the batch's software versions, so the seed
engine is the only difference; 100 threads (4 × 24 + 4), each run capped at 55 min
(`_cenrep_batch/bench_fastan/`, memo there).
- **Seed stage: 20–30 s** per complex at 24 threads, vs TRF 239 s (*A. suecica*, 25 threads), 279 s
  (*Brachypodium*, 16), 314 s (*A. kamchatica*, 16); 6 vs 3 s on the 74-Mb *Cryptococcus* (the FasTAN self-test and
  4 threads). *A. arenosa* alone: 90 vs ~540 CPU-s.
- **Calls** (Cryptococcus, *A. suecica*, Allo738, *Brachypodium*; 86 sequences called by either run):
  | | agree |
  |---|---|
  | same evidence | 84/86 |
  | same confidence | 83/86 |
  | calls overlap (both called) | 85/85 |
  | Jaccard ≥ 0.9 | 77/85 (median 1.000) |
  | main family consensus ≥ 95% identical | 77/80 |
  - *A. arenosa* / *A. suecica* Ar / Allo738 Ar: pAa on 24/24 chromosomes at the same place; pAa copies within 0.1%
    (26,297 vs 26,319; 26,908 vs 26,917; 16,074 vs 16,077).
  - *Brachypodium*, CENH3-validated copy numbers (Chen 2024): CentBs in BhS 86,247 (TRF 87,173; published 86,896);
    CentBd in BhD 18,198 (TRF 18,120; published 17,838).
- **Differences** (all reproducible from `calls.tsv`):
  - CentBd is split into 3 families (sat03/08/09_156, one chain group, 94.9–96.8% identical) instead of 1. On
    chr04/05 CentBd forms a dimeric HOR (adjacent monomers ~80% identical); TRF reports the region at 156, 310 and
    620 bp at once, FasTAN at one period per stretch, so a 155-bp piece mixing both monomer types founds its own
    cluster (84% to the main representative, < 90%). After refinement the two models differ by 5%, under the 97%
    merge rule. Calls are unchanged (BhD Jaccard 1.0); labels differ on 2 BhD and 2 Bd21 chromosomes.
  - Bd21 v3.0 (under-assembled; 0.4 Mb of CentBd in all): boundaries shift on 4 chromosomes (Jaccard 0.23–0.76),
    chr02Bd moves from a 39-kb, TE-derived 1,549-bp family to the CRM core (both medium), chr05Bd drops to medium
    (low_satellite; Chen 2024: ~40 kb of satellite there, CENH3 on TEs).
  - TAIR10 chr02 (CEN arrays mostly unassembled): called on CEN180 (27.7 kb) instead of a 97-bp family (55.8 kb), same
    place, medium instead of high.
  - *Cryptococcus* D chr09: one more long-unit call (see above).
- **Run time** (not attributable beyond the seed stage: the batch ran at 16 threads under other load):
  *Brachypodium* 10.9 min (TRF run 29), *A. suecica* 16.1, Allo738 17.4 (33).
- ***A. kamchatica* hit the 55-min cap** in annotation, 3 of 4 genomes done (seeds 30 s, families 9.4 min, then 8–15
  min per genome vs 5–6 in its TRF run). The FasTAN family set holds 15 families with ≥ 1-kb units vs 11: fewer candidates competed for
  the 40-family cap (34 vs 48), so 22-kb long-unit families got in. Long models cost 2–19 s per 5 Mb vs 0.2 s for a
  179-bp one, which explains part of the slowdown; the rest is not explained (no controlled rerun). Seeds there:
  FasTAN covers 84–88% of TRF seed bp, TRF 98–99% of FasTAN's; monomers in the 3 annotated genomes within 0.8%
  (106,716 vs 105,899; 105,743 vs 104,909; 112,443 vs 112,236).
