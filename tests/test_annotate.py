import pandas as pd
import pytest

from cenrep.annotate import MON_COLS, chain_arrays, drop_multimer_hits, resolve_overlaps
from conftest import have


def _mon(rows):
    return pd.DataFrame(rows, columns=MON_COLS)


def test_resolve_overlaps_keeps_best_and_adjacent():
    rows = [("c", 0, 150, "+", "f1", 100.0, 1e-20, 1, 150, 150, 0.9),
            ("c", 2, 150, "+", "f2", 80.0, 1e-15, 1, 148, 150, 0.85),   # same monomer, other family: dropped
            ("c", 145, 295, "+", "f1", 90.0, 1e-18, 1, 150, 150, 0.9)]   # 5-bp overlap only: kept
    out = resolve_overlaps(_mon(rows))
    assert out.family.tolist() == ["f1", "f1"] and out.start.tolist() == [0, 145]


def test_chain_arrays_splits_on_gap_and_strand(test_params):
    rows = [("c", i * 150, (i + 1) * 150, "+", "f", 90.0, 1e-20, 1, 150, 150, 0.9) for i in range(5)]
    rows += [("c", 1150, 1300, "+", "f", 90.0, 1e-20, 1, 150, 150, 0.9)]          # gap 400 > 300: singleton
    rows += [("c", 2000 + i * 150, 2150 + i * 150, "-", "f", 90.0, 1e-20, 1, 150, 150, 0.8) for i in range(3)]
    mon, arr = chain_arrays(_mon(rows), {"f": 150}, test_params)
    assert sorted(arr.n_monomers.tolist()) == [3, 5]
    assert set(arr.strand) == {"+", "-"} and len(mon) == 8


@pytest.mark.skipif(not have("trf", "blastn", "mafft"), reason="tools missing")
def test_annotate_synthetic_centromere(synth_annot, synth_complex, synth_families):
    mon, arr = synth_annot
    fam, _ = synth_families
    cen_fam = fam[(fam.consensus_len - 150).abs() <= 3].fam_id.iloc[0]
    s, e = synth_complex[1]["genomes"]["Pa_A"]["Pa_A01"]["cen"]
    a = arr[(arr.chrom == "Pa_A01") & (arr.family == cen_fam)]
    assert abs(a.start.min() - s) <= 150 and abs(a.end.max() - e) <= 150
    assert (a.strand == "-").any()                                    # the inverted block
    assert len(mon[(mon.chrom == "Pa_A01") & (mon.family == cen_fam)]) >= 250
    assert mon.mon_idx.tolist() == list(range(len(mon)))


@pytest.mark.skipif(not have("trf", "blastn", "mafft"), reason="tools missing")
def test_chunked_search_matches_whole_chromosome(synth_complex, synth_families):
    """Speed fix: nhmmer parallelises over query models only, so long models left 31 of 32 threads idle.
    Chunked search (chunks x all models on a thread pool, E-values rescaled to the chromosome) must give
    the same monomer calls as one whole-chromosome search."""
    import pyhmmer
    import pysam
    from cenrep.annotate import search_chrom, search_genome
    from cenrep.genomes import discover_genomes
    fam, wd = synth_families
    with pyhmmer.plan7.HMMFile(wd / "families.hmm") as f:
        hmms = list(f)
    g = [x for x in discover_genomes(synth_complex[0]) if x.name == "Pa_A"][0]
    with pysam.FastaFile(str(g.fasta)) as fa:
        whole = resolve_overlaps(_mon(search_chrom(hmms, "Pa_A01", fa.fetch("Pa_A01"), 4, 1e-5)))
    chunked = search_genome(hmms, g, 4, 1e-5, chunk=40_000, overlap=2_000)["Pa_A01"]
    key = lambda d: set(zip(d.start, d.end, d.strand, d.family))
    a, b = key(whole), key(chunked)
    assert len(a ^ b) <= 0.01 * len(a), (len(a), len(b), len(a ^ b))


def test_drop_multimer_hits_prefers_shortest_unit():
    """Real data: composite/multimer families (CentBd x7 + spacer = 1549 bp; CentO 155 x2.3 = 353 bp) whose long
    HMM hits outscore the monomer hits on the monomer arrays. A hit is dropped when hits of families with a
    >= 1.9x shorter model cover >= 50% of it; a long hit without such cover is kept."""
    rows = [("c", 0, 1549, "+", "long", 900.0, 1e-200, 1, 1549, 1549, 0.9)]
    rows += [("c", 170 + 157 * i, 327 + 157 * i, "+", "mono", 120.0, 1e-30, 1, 157, 157, 0.92) for i in range(7)]
    rows += [("c", 5000, 6549, "+", "long", 900.0, 1e-200, 1, 1549, 1549, 0.9)]        # no monomer cover
    rows += [("c", 5100, 5257, "+", "mono", 120.0, 1e-30, 1, 157, 157, 0.92)]           # 10% cover only
    out = drop_multimer_hits(_mon(rows))
    assert (out.family == "mono").sum() == 8
    assert out[out.family == "long"].start.tolist() == [5000]


def test_phase_shifted_hits_are_not_double_counted():
    """Real data (A. thaliana, 4 phase-shifted CEN180 subfamilies): hits overlapping each of two accepted
    neighbours by <= 50% survived the pairwise clash rule, inflating monomer counts by ~27% and putting
    overlapping monomers into HOR strings. A hit is accepted only if <= 20% of it is covered by already
    accepted (higher-scoring) hits."""
    rows = [("c", 0, 178, "+", "A", 100.0, 1e-30, 1, 178, 178, 0.95),
            ("c", 178, 356, "+", "A", 100.0, 1e-30, 1, 178, 178, 0.95),
            ("c", 89, 267, "+", "B", 95.0, 1e-28, 1, 178, 178, 0.93),        # half-phase shifted: 100% covered
            ("c", 150, 285, "+", "C", 60.0, 1e-15, 20, 154, 178, 0.9),       # bridging partial
            ("c", 350, 528, "+", "B", 90.0, 1e-27, 1, 178, 178, 0.92)]       # 6 bp overlap: kept
    out = resolve_overlaps(_mon(rows))
    assert list(zip(out.start, out.family)) == [(0, "A"), (178, "A"), (350, "B")]


def test_chain_arrays_across_related_families(test_params):
    """Interleaved monomers of related subfamilies form one array (majority family) instead of the minority ones
    being dropped as singletons."""
    fams = ["A", "A", "B", "A", "A", "B", "A"]
    rows = [("c", 155 * i, 155 * i + (165 if f == "B" else 155), "+", f, 90.0, 1e-20, 1, 155, 155, 0.9)
            for i, f in enumerate(fams)]
    mon, arr = chain_arrays(_mon(rows), {"A": 155, "B": 165}, test_params, groups={"A": "A", "B": "A"})
    assert len(arr) == 1 and arr.n_monomers.iloc[0] == 7 and arr.family.iloc[0] == "A" and len(mon) == 7
    mon, arr = chain_arrays(_mon(rows), {"A": 155, "B": 165}, test_params)
    assert (mon.family == "B").sum() == 0                                  # per-family chaining drops them


@pytest.mark.skipif(not have("trf", "blastn", "mafft"), reason="tools missing")
def test_evalue_ref_len_scales_short_sequences(synth_complex, synth_families):
    """Contigs: E-values scale with sequence length, so a 1-Mb contig got a cutoff ~8-13 bits looser than a
    100-Mb chromosome. evalue_ref_len scales every shorter sequence as if it had that length."""
    import pyhmmer
    from cenrep.annotate import search_genome
    from cenrep.genomes import discover_genomes
    fam, wd = synth_families
    with pyhmmer.plan7.HMMFile(wd / "families.hmm") as f:
        hmms = list(f)
    g = [x for x in discover_genomes(synth_complex[0]) if x.name == "Pa_A"][0]
    run = lambda n: search_genome(hmms, g, 4, 1e-5, chunk=100_000, overlap=2_000, ref_len=n)["Pa_A01"]
    base, same, big = run(0), run(1_000), run(3_000_000_000)
    key = lambda d: dict(zip(zip(d.start, d.end, d.family), d.evalue))
    assert key(same) == key(base)                                    # shorter than the chromosome: no effect
    b, s = key(base), key(big)
    shared = set(b) & set(s)
    assert shared and all(abs(s[k] / b[k] / 10_000 - 1) < 1e-6 for k in shared if b[k] > 0)
    assert (big.evalue <= 1e-5).all()
