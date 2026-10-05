# Published centromeric satellites — validation reference (web research 2026-09-25)

Sources were read directly (PMC full text, supplements, GenBank records).
[inference] = the researcher's reasoning, not a published statement.

## Arabidopsis

| genome | satellite(s) | monomer | evidence | size / HOR | exceptions |
|---|---|---|---|---|---|
| A. thaliana | CEN180 (AthCEN178; pAL1/AtCon), minor AthCEN159 (was CEN160) | ~178 bp, 39% GC | CENH3 ChIP-seq (Naish 2021) | Col-CEN 66,131 copies (11,848–15,613 per chr); HORs mean 2.41 monomers, max 60 (10,689 bp), chromosome-specific. 66 accessions: arrays mean 2.91 Mb per chr (1.11–6.49) (Wlodzimierz 2023) | CEN5 and a CEN4 block ATHILA-rich; CEN3 = 2 blocks on opposite strands |
| A. arenosa | pAa | 166–179 bp, AT-rich; 178.7±1.7 bp, 73% identity to Ath (Hall 2005) | FISH, all 32 cen (Kamm 1995) | HORs present | — |
| A. lyrata | AlyCEN168 + AlyCEN179 (pAge / pAa families) | ~168 / ~179 bp | FISH; HiFi arrays; no ChIP | similar to Ath | chromosome-specific: AlyCEN168 dominates CEN2/5/6/8, AlyCEN179 CEN1/3/4/7 (Wlodzimierz 2023) |
| A. halleri | pAa, pAge1, pAge2 | ~180 bp; ~10% divergence within, ~30% between families; ~50 bp shared core | FISH | — | chromosome-specific and strain-polymorphic (Kawabe & Nasuda 2005, 2006) |
| A. suecica | CEN180 (T) + pAa (A) | 58–80% identical to each other | genome-specific FISH; both CENH3s load on both subgenomes | — | pAL1 reorganized in the hybrid, pAa not (Kamm) |
| A. kamchatica | nothing published; expect the lyrata + halleri families | | | | |

[inference] GenBank clone lengths from Kawabe & Charlesworth 2007: ge1 ≈ 168 bp multiples, ge2 ≈ 175 bp, pAa ≈ 178 bp. So AlyCEN168 is most likely pAge-type.

A. thaliana CEN178 template, verbatim from the TRASH README (labelled 178, 177 nt as published):
`AGTATAAGAACTTAAACCGCAACCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCATATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGATGTCATGTGTATGATTG`

## Capsella

- C. rubella: unnamed satellite, 168.4±1.0 bp, 32% GC. FISH hybridized every DAPI-bright region; HOR-like; no similarity to CEN180 (Hall 2005; AY640635–AY640681).
- C. grandiflora, C. orientalis, C. bursa-pastoris: nothing published. Recent assemblies (Kasianova 2025; Penin 2024) don't characterize one.

## Brassica

- **A genome** (B. rapa, and the A subgenomes of B. napus and B. juncea):
  - CentBr1 (8 chromosomes) + CentBr2 (A03, A05) = "CEN176", 176 bp, ≤ 86.9% identity to each other. Plus the CRB retrotransposon.
  - Evidence: ChIP (Wang 2011); CENH3 ChIP-seq (Chen 2025).
  - Size: 7 T2T assemblies, 3.19 M copies; arrays mean 6.85 Mb per chromosome (2.74–13.12). CEN176 is 65% of centromere sequence.
- **C genome** (B. oleracea, and the C subgenomes of B. napus and B. carinata):
  - CentBr1 + CentBr2, 176 bp, 40% GC.
  - Layered: a CentBr1 core with CentBr2 flanks on 5 chromosomes; C03/C04/C08/C09 carry CentBr2 only.
- **B genome** (B. nigra, and the B subgenomes of B. juncea and B. carinata):
  - No satellite. Cores are ALE (Ty1/copia) + CRM + CRB, flanked by islands of the B-specific pBN35 (X16588, 496 bp).
  - A 28-copy CEN176 relic sits at one B-centromere boundary.
  - Evidence: CENH3 ChIP (Wang 2019; Perumal 2020). B-subgenome centromeres average 5.46 Mb (Chen 2025).
  - B-derived ChIP clusters spread onto A chromosomes in B. juncea, but not onto C chromosomes in B. carinata.

## Brachypodium

- **B. distachyon** — CentBd, 156 bp, 49% GC. ChIP (Li 2018; Chen 2024).
  - Bd21: 2.18 Mb in total, 535–3,890 copies per chromosome; CENH3 domains 0.74–1.25 Mb; HORs present.
  - **Bd5 has only ~40 kb of satellite; its CENH3 domain is TE-rich.**
- **B. stacei** — CentBs, 157 bp, 19 SNPs + 1 indel different from CentBd.
  - Bst99: 15.98 Mb, 7,652–14,348 copies per chromosome; centromeres 1.40–2.50 Mb.
  - Arrays also occur outside the centromeres on most chromosomes.
- **B. hybridum** — both satellites, subgenome-specific: 17,838 D copies, 86,896 S copies.
  - BhS centromeres 0.80–2.40 Mb, BhD 0.55–1.30 Mb.
  - Some interior CentBd arrays lack CENH3.
- CRBd1/2: KF040483/KF040484.

## Coffea (C. arabica, C. canephora, C. eugenioides)

- No short satellite has been reported. Putative centromeres are dense in CRC retrotransposons, with "complete absence of satellite arrays" (de Castro Nunes 2018). No CENH3 ChIP, so this is putative.
- C. arabica: no tandem units shorter than 2 kb. Heterochromatin is made of arrays with units longer than 2 kb, derived from CRM/Tekay/Athila. Example: a 690 kb array of a 2,683 bp CRM-derived monomer on Chr7c (Scalabrin 2024).
- Coffea_sat11 (82 bp) is pericentromeric, not centromeric (Cintra 2022).

## Oryza

- **O. punctata (BB):** CentO, a major 165 bp form (10 nt insertion) + a 155 bp subfamily. CENH3 immunostaining shows 12 signals (Zhang 2005). Centromeres 1.23–1.67 Mb; ~17.55 Mb in total.
- **O. officinalis (CC):**
  - CentO-C (C1), 126 bp.
  - TrsC/CentO-C2, 366 bp — also subtelomeric, so not a centromere marker.
  - CentO-C is found on only 10 centromeres: **CEN7 carries CentO and CEN2 has no satellite** (Bao 2006; Lee 2005).
- **O. minuta (BBCC):** nothing published. Expect CentO (B) + CentO-C (C).
  - [inference] Wang 2025 PNAS did not place centromeres on C02 or Ct02, which fits a satellite-free CEN2.
- An ~80 bp region is conserved among CentO-C1, CentO and maize CentC.

## Erythranthe / Mimulus

- **E. guttata:** Cent728, 728 bp, 37% GC. FISH + drive-locus genetics (Fishman & Saunders 2008); arrays sit at the centre of every chromosome (Lovell 2025). No ChIP.
  - Size: IM62 v3 41.32 Mb, IM767 31.83 Mb (BLAST hits merged when > 1.2 kb).
  - The driving LG11 D haplotype has 2 dense Cent728 arrays flanking a ~12 Mb gene-rich, Cent728-depleted region.
  - No GenBank record for the consensus.
- **E. nasuta:** Cent728, 27.13 Mb.
- **M. sookensis:** nothing published; single origin from guttatus × nasutus.
- A 198 bp M. guttatus repeat listed in Melters 2013 Table S2 is unverified.

## References (DOIs)

Naish 2021 Science 10.1126/science.abi7489; Wlodzimierz 2023 Nature 10.1038/s41586-023-06062-z;
Kamm 1995 PMB 10.1007/BF00037014; Hall 2005 Genetics 10.1534/genetics.104.038208;
Kawabe & Nasuda 2005 MGG 10.1007/s00438-004-1081-x; Kawabe & Nasuda 2006 Genetica 10.1007/s10709-005-0718-y;
Kawabe & Charlesworth 2007 JME 10.1007/s00239-006-0097-8; Berr 2006 Plant J 10.1111/j.1365-313X.2006.02912.x;
Comai 2003 Chrom Res 10.1023/a:1022883709060; Karimi-Ashtiyani 2024 PMB 10.1007/s11103-024-01474-5;
Melters 2013 Genome Biol 10.1186/gb-2013-14-1-r10; Wang 2011 Chromosoma 10.1007/s00412-011-0315-z;
Wang 2019 Chromosoma 10.1007/s00412-019-00701-z; Perumal 2020 Nat Plants 10.1038/s41477-020-0735-y;
Chen W 2025 Nat Plants 10.1038/s41477-025-02131-5; Li 2018 Plant J 10.1111/tpj.13832;
Chen C 2024 Genome Biol 10.1186/s13059-024-03206-w; de Castro Nunes 2018 Front Plant Sci 10.3389/fpls.2018.00175;
Scalabrin 2024 Nat Commun 10.1038/s41467-023-44449-8; Cintra 2022 Genome 10.1139/gen-2021-0045;
Zhang W 2005 Plant Physiol 10.1104/pp.105.064147; Lee 2005 PNAS 10.1073/pnas.0503863102;
Bao 2006 MGG 10.1007/s00438-006-0103-2; Wang X 2025 PNAS 10.1073/pnas.2424854122;
Fishman & Saunders 2008 Science 10.1126/science.1161406; Lovell 2025 Mol Ecol Resour 10.1111/1755-0998.70012
