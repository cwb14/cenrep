"""Short satellite units: copies too short to be significant alone are kept as runs (annotate.run_significant).

A lone 40-bp model scores at most ~20 bits and can never reach E 1e-5 over a chromosome (Avena 38-51-bp families:
0 bp annotated from seeds spanning up to 159 Mb)."""
import random

import numpy as np
import pandas as pd
import pysam
import pytest

from cenrep.annotate import MON_COLS, run_significant, search_genome
from cenrep.families import build_hmm, single_hmm
from cenrep.genomes import discover_genomes
from conftest import have
from synth import mutate, rseq

U = 40
CHROM = 500_000_000          # E-values as on a 500-Mb chromosome (search_genome ref_len)


def _array(rng, unit, n, rate=0.08):
    out = []
    for _ in range(n):
        m = list(mutate(rng, unit, rate))
        if rng.random() < 0.2:                                   # occasional 1-bp indel
            i = int(rng.integers(1, len(m) - 1))
            m.insert(i, "A") if rng.random() < 0.5 else m.pop(i)
        out.append("".join(m))
    return "".join(out)


def _hits(rows, to=U):
    return pd.DataFrame([("c", s, e, st, "f", 10.0, ev, 1, to, U, 0.9) for s, e, st, ev in rows], columns=MON_COLS)


def test_run_significance():
    L = 500_000_000
    run = [(i * U, (i + 1) * U, "+", 100.0) for i in range(5)]          # 5 adjacent weak copies: E 1e10 x (81/L)^4
    lone = [(10_000, 10_040, "+", 1e-3), (20_000, 20_040, "+", 1e-7)]   # a weak and a strong single hit
    far = [(30_000, 30_040, "+", 100.0), (30_200, 30_240, "+", 100.0)]  # gap 160 > 2 units: not a run
    flip = [(40_000, 40_040, "+", 100.0), (40_040, 40_080, "-", 100.0)]  # other strand: not a run
    ok = run_significant(_hits(run + lone + far + flip), L, 1e-5, 2.0)
    assert ok.tolist() == [True] * 5 + [False, True] + [False, False] + [False, False]
    frag = run_significant(_hits(run, to=U // 2), L, 1e-5, 2.0)         # half-copies: an array is made of whole copies
    assert not frag.any()


def _genome(tmp_path, seq):
    fa = tmp_path / "Sx_A.fa"
    fa.write_text(f">c\n{seq}\n")
    pysam.faidx(str(fa))
    return discover_genomes(tmp_path)[0]


def _cov(df, s, e):
    c = np.zeros(e - s, bool)
    for a, b in zip(df.start, df.end):
        c[max(a, s) - s:max(min(b, e) - s, 0)] = True
    return c.mean()


def test_short_unit_array_is_found_as_a_run(tmp_path):
    rng = np.random.default_rng(5)
    unit = rseq(rng, U)
    arr = _array(rng, unit, 150)
    s = 300_000
    g = _genome(tmp_path, rseq(rng, s) + arr + rseq(rng, 300_000))
    e = s + len(arr)
    mon = search_genome([single_hmm("f", unit)], g, 1, 1e-5, ref_len=CHROM)["c"]
    alone = mon[mon.evalue <= 1e-5]
    assert _cov(alone, s, e) < 0.1                    # no copy is significant on its own
    assert _cov(mon, s, e) >= 0.9                     # the run is
    assert ((mon.start >= s - U) & (mon.end <= e + U)).all()               # nothing outside the array


def test_no_runs_in_random_sequence(tmp_path):
    rng = np.random.default_rng(9)
    g = _genome(tmp_path, rseq(rng, 1_000_000))
    models = [single_hmm(f"f{i}", rseq(rng, n)) for i, n in enumerate((30, 40, 60, 150))]
    assert search_genome(models, g, 1, 1e-5, ref_len=CHROM)["c"].empty


@pytest.mark.skipif(not have("blastn", "mafft"), reason="tools missing")
def test_short_unit_gets_a_profile_model(tmp_path, test_params):
    rng = np.random.default_rng(7)
    unit = rseq(rng, U)
    fa = tmp_path / "g.fa"
    arr = _array(rng, unit, 400, 0.05)
    fa.write_text(f">c\n{rseq(rng, 1000)}{arr}{rseq(rng, 1000)}\n")
    pysam.faidx(str(fa))
    hmm, cons, mode, n = build_hmm("f", unit, [("g", "c", 1000, 1000 + len(arr))], {"g": fa}, tmp_path / "w", 2,
                                   test_params, random.Random(1))
    assert mode == "profile" and n >= test_params.min_hmm_seqs and abs(hmm.M - U) <= 1
