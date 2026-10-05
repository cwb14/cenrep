"""Synthetic 4-genome complex with planted centromeric arrays, HORs, decoys, a gap and a long-unit array.

Each chromosome is L bp. Chromosomes 1-3 carry a centromere core starting at CEN_START; chromosome 4 has none.
Both ends of every chromosome carry a 90-bp subtelomeric decoy array (40 copies).
"""
from pathlib import Path

import numpy as np
import pysam

L = 300_000
CEN_START = 130_000
MON = 150


def rseq(rng, n: int) -> str:
    return "".join(rng.choice(list("ACGT"), size=n, p=[0.31, 0.19, 0.19, 0.31]))


def mutate(rng, s: str, rate: float) -> str:
    a = list(s)
    for i in np.where(rng.random(len(a)) < rate)[0]:
        a[i] = rng.choice([b for b in "ACGT" if b != a[i]])
    return "".join(a)


def rc(s: str) -> str:
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def build_complex(outdir: Path, seed: int = 11) -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    shared = np.random.default_rng(seed)
    cen, dec = rseq(shared, MON), rseq(shared, 90)
    hor4 = [mutate(shared, cen, 0.08) for _ in range(4)]
    dup3 = [mutate(shared, cen, 0.08) for _ in range(3)]
    u3k = rseq(shared, 3000)
    truth = {"cen": cen, "dec": dec, "genomes": {}}
    for gi, gname in enumerate(["Pa_A", "Pb_B", "Tx_A", "Tx_B"]):
        rng = np.random.default_rng(seed + 100 + gi)

        def plain(n):
            return "".join(mutate(rng, cen, 0.08) for _ in range(n))

        def decoy():
            return "".join(mutate(rng, dec, 0.02) for _ in range(40))

        chroms, t = {}, {}
        for ci in range(1, 5):
            name = f"{gname}{ci:02d}"
            if ci == 1:
                hor = "".join(mutate(rng, v, 0.002) for _ in range(15) for v in hor4)
                parts = [("sat", plain(60)), ("hor", hor), ("sat", plain(60)), ("inv", rc(plain(60))), ("sat", plain(20))]
            elif ci == 2:
                d1 = "".join(mutate(rng, v, 0.002) for v in dup3)
                d2 = "".join(mutate(rng, v, 0.002) for v in dup3)
                parts = [("sat", plain(100)), ("dup", d1), ("sat", plain(60)), ("dup", d2), ("sat", plain(100))]
            elif ci == 3:
                lu = "".join(mutate(rng, u3k, 0.03) for _ in range(6))
                parts = [("sat", plain(120)), ("micro", "AT" * 300), ("sat", plain(120)),
                         ("gap", "N" * 500), ("lu", lu), ("sat", plain(60))]
            else:
                parts = []
            core = "".join(s for _, s in parts)
            left = decoy()
            left += rseq(rng, CEN_START - len(left))
            right = decoy()
            seq = left + core + rseq(rng, L - CEN_START - len(core) - len(right)) + right
            assert len(seq) == L
            feats, pos = {}, CEN_START
            for kind, s in parts:
                feats.setdefault(kind, []).append((pos, pos + len(s)))
                pos += len(s)
            t[name] = {"cen": (CEN_START, CEN_START + len(core)) if parts else None, **feats}
            chroms[name] = seq
        fa = outdir / f"{gname}.fa"
        with open(fa, "w") as f:
            for n, s in chroms.items():
                f.write(f">{n}\n")
                for i in range(0, len(s), 60):
                    f.write(s[i:i + 60] + "\n")
        pysam.faidx(str(fa))
        truth["genomes"][gname] = t
    return truth
