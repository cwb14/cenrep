import pysam


def test_synthetic_complex_layout(synth_complex):
    cdir, truth = synth_complex
    fas = sorted(cdir.glob("*.fa"))
    assert len(fas) == 4 and all((f.parent / (f.name + ".fai")).exists() for f in fas)
    t = truth["genomes"]["Pa_A"]
    assert t["Pa_A01"]["cen"] == (130_000, 130_000 + 260 * 150)
    assert t["Pa_A04"]["cen"] is None
    s, e = t["Pa_A03"]["gap"][0]
    with pysam.FastaFile(str(cdir / "Pa_A.fa")) as fa:
        assert set(fa.fetch("Pa_A03", s, e)) == {"N"}
