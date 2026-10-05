import numpy as np
import pytest

from cenrep.families import (RRNA_CM, absorb_multimers, cluster_seeds, dedup_families, frame_rotation, msa_consensus,
                             chain_groups, pad_hit, rdna_tags, seed_match, trim_overhang)
from conftest import have
from synth import mutate, rseq


def test_seed_match_rotation_same_period():
    assert seed_match(150, 150, [(97.0, 1, 150, 40, 189)], 0.9, 0.8)


def test_seed_match_rejects_low_identity_same_period():
    assert not seed_match(150, 150, [(88.0, 1, 150, 40, 189)], 0.9, 0.8)


def test_seed_match_multimer_needs_tiling_and_integer_ratio():
    tiled = [(85.0, 1, 150, 1 + 150 * i, 150 * (i + 1)) for i in range(6)]   # 6 copies on a 900-bp dimer
    assert seed_match(150, 450, tiled, 0.9, 0.8)
    assert not seed_match(150, 450, tiled[:1], 0.9, 0.8)      # covers 1/3 of the circle
    assert not seed_match(150, 400, tiled, 0.9, 0.8)          # 400/150 is not an integer ratio


def test_pad_hit_extends_to_full_model():
    # + strand: hit covers model 3..148 of 150 at target 101..246 -> pad 2 left, 2 right
    assert pad_hit(101, 246, "+", 3, 148, 150, 1000) == (99, 248)
    # - strand: model 1 sits at the high coordinate; model 3..148 at 101..246 -> 2 right, 2 left
    assert pad_hit(101, 246, "-", 3, 148, 150, 1000) == (99, 248)
    # clipped to the sequence
    assert pad_hit(1, 146, "+", 5, 150, 150, 1000) == (1, 146)


def test_msa_consensus_drops_gappy_columns():
    assert msa_consensus(["AC-GT", "ACTGT", "AC-GA"]) == "ACGT"


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_cluster_seeds_rotation_multimer_unrelated(tmp_path):
    rng = np.random.default_rng(5)
    a, b = rseq(rng, 150), rseq(rng, 150)
    seqs = {"A": a,
            "Arot": mutate(rng, a[40:] + a[:40], 0.02),
            "A3": "".join(mutate(rng, a, 0.05) for _ in range(3)),
            "B": b,
            "Asub": mutate(rng, a, 0.14)}
    assign = cluster_seeds(seqs, ["A", "Arot", "A3", "B", "Asub"], tmp_path, 2, 0.90, 0.80)
    assert assign["Arot"] == "A" and assign["A3"] == "A"
    assert assign["B"] == "B" and assign["Asub"] == "Asub"


@pytest.mark.skipif(not have("trf", "blastn", "mafft"), reason="tools missing")
def test_build_families_on_synthetic(synth_families):
    fam, workdir = synth_families
    near = lambda n: fam[(fam.consensus_len - n).abs() <= 3]
    assert len(near(150)) == 1 and near(150).hmm_mode.iloc[0] == "profile"   # CEN + its HOR multimer
    assert len(near(90)) == 1                                                  # subtelomeric decoy
    assert len(near(600)) == 0                                                 # multimer merged into CEN
    assert (workdir / "families.hmm").stat().st_size > 0


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_absorb_multimers_folds_divergent_hor_multimers(tmp_path):
    """TRF often reports HOR-rich arrays at the HOR period; when some monomers of the HOR copy are too
    divergent for the 80% multimer rule, the multimer becomes its own family (real data: CEN180 x5/x7).
    Absorption at superfamily identity folds it back into its monomer family."""
    rng = np.random.default_rng(8)
    a, b = rseq(rng, 180), rseq(rng, 180)
    cons = {"A": a,
            "A5": "".join(mutate(rng, a, 0.12 if i % 2 else 0.25) for i in range(5)),   # some copies ~75%
            "B": b,
            "B2": "".join(mutate(rng, b, 0.22) for _ in range(2)),
            "C": rseq(rng, 900)}                                                     # unrelated 900 bp
    parent = absorb_multimers(cons, tmp_path, 2, 0.70)
    assert parent == {"A5": "A", "B2": "B"}


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_monomer_must_not_join_multimer_rep_below_family_identity(tmp_path):
    """Review C1: TRF often calls the largest array at a multimer period, so the first representative is a
    dimer; a monomer of a related but distinct subfamily (87%, CentBd/CentBs-like) must not join it at the
    80% multimer threshold. Result must not depend on seed order."""
    rng = np.random.default_rng(12)
    a = rseq(rng, 160)
    seqs = {"A2": a + mutate(rng, a, 0.03), "B1": mutate(rng, a, 0.13), "A1": mutate(rng, a, 0.02)}
    for order in (["A2", "B1", "A1"], ["A1", "B1", "A2"]):
        assign = cluster_seeds(seqs, order, tmp_path / "_".join(order), 2, 0.90, 0.80)
        assert assign["B1"] == "B1", order
        assert assign["A1"] == assign["A2"]


@pytest.mark.skipif(not have("cmsearch", "cmemit"), reason="infernal missing")
def test_rdna_tags_flags_5s_rdna_units(tmp_path):
    """Review C2: pericentromeric 5S rDNA arrays look like satellites; tag them from the Rfam 5S model."""
    import subprocess
    out = subprocess.run(["cmemit", "-c", "--seed", "1", str(RRNA_CM)], capture_output=True, text=True).stdout
    five_s = "".join(out.split(">5_8S_rRNA")[0].splitlines()[1:]).upper().replace("U", "T")
    rng = np.random.default_rng(3)
    unit = rseq(rng, 200) + five_s + rseq(rng, 180)                   # ~500-bp 5S rDNA repeat unit
    tags = rdna_tags({"f5s": unit, "sat": rseq(rng, 178)}, tmp_path, 2)
    assert tags == {"f5s": "5S_rRNA", "sat": "."}


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_frame_rotation_aligns_related_family_to_reference(tmp_path):
    """Real data: related subfamilies (98-100% identical) got monomer frames 110/178 and 113/157 bp apart from
    independent lexicographic rotations, so their hits half-overlapped and 10-18% of array bp stayed unannotated.
    A related family is rotated (and oriented) into the frame of the most similar earlier family."""
    from cenrep.util import canonical_rotation, revcomp
    rng = np.random.default_rng(21)
    a = rseq(rng, 180)
    b = mutate(rng, a[70:] + a[:70], 0.05)
    for i, q in enumerate((b, revcomp(b))):
        c = frame_rotation(q, {"A": a}, tmp_path / f"q{i}", 1, 0.70)
        assert sum(x == y for x, y in zip(c[:40], a[:40])) >= 32, i
    u = rseq(rng, 180)
    assert frame_rotation(u, {"A": a}, tmp_path / "u", 1, 0.70) == canonical_rotation(u)
    assert frame_rotation(u, {}, tmp_path / "e", 1, 0.70) == canonical_rotation(u)


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_dedup_families_merges_near_identical(tmp_path):
    rng = np.random.default_rng(22)
    x = rseq(rng, 178)
    cons = {"f1": x, "f2": mutate(rng, x[50:] + x[:50], 0.02), "f3": mutate(rng, x, 0.13), "f4": rseq(rng, 178)}
    assert dedup_families(cons, ["f1", "f2", "f3", "f4"], tmp_path, 1, 0.95) == {"f2": "f1"}
    # real data: CentBr1 vs CentBr2 consensus ~95% identical must stay separate at the 0.97 default
    from cenrep.params import Params
    pair = {"b1": x, "b2": mutate(np.random.default_rng(23), x, 0.05)}
    assert dedup_families(pair, ["b1", "b2"], tmp_path / "p", 1, Params().dedup_ident) == {}


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_trim_overhang_cuts_terminal_redundancy(tmp_path):
    """Real data (Msookensis): a 796-bp family = Cent728 unit + the first 68 bp of the next unit (terminally
    redundant) outscored the real 730-bp model on every Cent728 array once both shared a frame. Such a consensus
    is trimmed to its true period; a normal monomer, an internal motif repeat and a short chance match are not."""
    rng = np.random.default_rng(31)
    u = rseq(rng, 730)
    assert trim_overhang(u + mutate(rng, u[:68], 0.02), tmp_path / "a", 1) == u
    assert trim_overhang(u, tmp_path / "b", 1) == u
    w = u[:500] + u[:40] + u[540:]                                     # internal 40-bp duplication at 500
    assert trim_overhang(w, tmp_path / "d", 1) == w
    v = rseq(rng, 300)
    assert trim_overhang(v + v[:12], tmp_path / "c", 1) == v + v[:12]


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_chain_groups_join_close_subfamilies_only(tmp_path):
    """Real data (Oryza): CentO-165 monomers interleaved among CentO-155 were dropped as per-family singletons.
    Arrays are chained per chain group: subfamilies >= 90% identical (length within 15%); CentBd/CentBs-like
    87% pairs and unrelated families stay apart."""
    rng = np.random.default_rng(41)
    a = rseq(rng, 155)
    b = mutate(rng, a[:80] + "ATTTATAGGC" + a[80:], 0.01)          # 165-bp insertion variant
    c = mutate(rng, a, 0.13)
    d = rseq(rng, 155)
    g = chain_groups({"a": a, "b": b, "c": c, "d": d}, ["a", "b", "c", "d"], tmp_path, 1, 0.90)
    assert g["a"] == g["b"] == "a" and g["c"] == "c" and g["d"] == "d"


@pytest.mark.skipif(not have("cmsearch", "cmemit"), reason="infernal missing")
def test_rdna_neighbourhood_flags_45s_spacer_repeats(tmp_path):
    """Real data (P. infirma chr1A NOR): the 45S rDNA intergenic-spacer subrepeat (Poa_sat18_125) forms satellite-like
    arrays without any rRNA gene in its consensus, so rdna_tags missed it and it was called centromeric. Arrays whose
    flanks carry 45S rRNA genes (18S/5.8S/25S) are rDNA-associated; a satellite with plain flanks is not."""
    import subprocess
    import pandas as pd
    import pysam
    from cenrep.families import rdna_neighbourhood
    out = subprocess.run(["cmemit", "-c", "--seed", "1", str(RRNA_CM)], capture_output=True, text=True).stdout
    recs = {r.splitlines()[0].split()[0]: "".join(r.splitlines()[1:]) for r in out.split(">")[1:]}
    ssu = next(s for n, s in recs.items() if n.startswith("SSU")).upper().replace("U", "T")
    rng = np.random.default_rng(5)
    igs, sat = rseq(rng, 125), rseq(rng, 146)
    seq, rows = [], []
    pos = 0
    for i in range(6):                                    # 6 rDNA units: 18S gene + spacer with an IGS array
        for part, fam in ((rseq(rng, 2000), None), (ssu, None), (rseq(rng, 3000), None), (igs * 20, "igs")):
            if fam:
                rows.append((len(rows), "c1", pos, pos + len(part), "+", fam, 20, 0.9, 125))
            seq.append(part)
            pos += len(part)
    for i in range(6):                                    # satellite arrays with random flanks
        seq.append(rseq(rng, 40_000))
        pos += 40_000
        rows.append((len(rows), "c1", pos, pos + 146 * 30, "+", "sat", 30, 0.9, 146))
        seq.append(sat * 30)
        pos += 146 * 30
    seq.append(rseq(rng, 40_000))
    fa = tmp_path / "g.fa"
    fa.write_text(">c1\n" + "".join(seq) + "\n")
    pysam.faidx(str(fa))
    from cenrep.annotate import ARR_COLS
    arrays = pd.DataFrame(rows, columns=ARR_COLS)
    res = rdna_neighbourhood(fa, {"c1": len("".join(seq))}, arrays, ["igs", "sat"], tmp_path / "w", 2,
                             flank=15_000, n_arrays=5)
    assert res["igs"] == (5, 5) and res["sat"] == (5, 0)
