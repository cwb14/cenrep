# Family roles and the all-arrays GFF (design, 2026-10-08; v0.8.6)

User request: report what cenrep finds along the way (rDNA, subtelomeres, knobs), and every array of every family.

## Roles (`classify.assign_roles`, one per family and genome, first match wins)
The roles reuse the tests classify already runs; there are no new thresholds.
- `rDNA`: `rdna_tag` = 5S_rRNA (Rfam on the consensus, families stage) or 45S_spacer (largest arrays within 15 kb of
  45S genes; the tag from one genome holds in all, as for eligibility).
- `centromeric`, `cen_associated`: classes from classify / centromeres.
- `subtelomeric`: the classify verdict, including the complex vote (docs/2026-10-08-subtelomeres.md).
- `minor`: < chrom_min_bp (20 kb) in this genome.
- `interstitial`: concentration >= min_concentration (one main cluster holds most of the family on each chromosome
  that carries >= 20 kb), not centromeric. Includes CRM-vetoed blocks (the *Populus* 148-bp knob, *Aegilops*
  interstitial blocks). It is not called "knob": without cytology a knob cannot be told from other heterochromatic
  blocks.
- `dispersed`: the rest (spread along chromosomes, or < 20 kb on every chromosome).
- `unassessed`: borrowed (contig) genomes, beyond rDNA / centromeric / cen_associated.

## Outputs
- `satellite_families.tsv`: new column `role` (absent rows: `absent`).
- `<genome>_satellite_arrays.gff`: every array (`work/annotate/<genome>.arrays.tsv`), `satellite_DNA`, attributes
  ID, family (final name), role, rdna (when tagged), monomer_len, n_monomers, mean_identity, in_centromere.
  The centromeric_repeat GFF is unchanged (centromeric families and arrays overlapping calls, plus monomers and HORs).
- Log: role counts per genome.

## Not reported
Telomeres (period < 30 bp is outside the search), rRNA genes themselves, long-unit arrays (work/centromeres/*.longunit.tsv),
and tandem repeats outside the <= 40 families kept per run.

Versions: classify 0.9.2, centromeres 0.9.1, output 0.4.0.
