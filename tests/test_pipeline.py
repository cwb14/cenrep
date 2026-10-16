from dataclasses import fields

import pytest

from cenrep.params import Params
from cenrep.util import CenrepError
from cenrep.pipeline import STAGE_PARAMS, STAGES


def test_every_parameter_belongs_to_exactly_one_stage():
    mapped = [x for s in STAGES for x in STAGE_PARAMS[s]]
    assert sorted(mapped) == sorted(f.name for f in fields(Params)) and len(mapped) == len(set(mapped))


def test_every_parameter_has_a_reason():
    from cenrep.params import REASON
    assert set(REASON) == {f.name for f in fields(Params)} and all(REASON.values())


def _complex(d, names, nseq=1):
    for n in names:
        chroms = [f"{n}_{i}" for i in range(nseq)]
        (d / f"{n}.fa").write_text("".join(f">{c}\nACGT\n" for c in chroms))
        (d / f"{n}.fa.fai").write_text("".join(f"{c}\t4\t{i}\t4\t5\n" for i, c in enumerate(chroms)))
    return d


@pytest.mark.parametrize("spec, msg", [
    ("X=a", "unknown genome"), ("a=X", "unknown genome"), ("a=a", "borrows from itself"),
    ("a=b,b=c", "is itself borrowed"), ("a=b,a=c", "more than once"), ("a", "TARGET=SOURCE"),
    ("a=", "TARGET=SOURCE"), ("=b", "TARGET=SOURCE"), ("a=b+", "TARGET=SOURCE")])
def test_classify_from_is_validated(tmp_path, spec, msg):
    from cenrep.pipeline import Pipeline
    cdir = _complex(tmp_path, "abc")
    with pytest.raises(CenrepError, match=msg):
        Pipeline(cdir, 1, tmp_path / "cache", Params(classify_from=spec))


def test_classify_from_is_normalised_and_sources_come_first(tmp_path):
    from cenrep.pipeline import Pipeline
    cdir = _complex(tmp_path, ["Ctg", "a", "b"])
    pl = Pipeline(cdir, 1, tmp_path / "cache", Params(classify_from=" Ctg = b+a+b "))
    assert pl.p.classify_from == "Ctg=a+b" and pl.borrow == {"Ctg": ["a", "b"]}
    assert [g.name for g in pl.genomes] == ["a", "b", "Ctg"]


def test_seed_engine_changes_every_fingerprint(tmp_path):
    """Seeds feed every stage: a checkpoint made with one engine is never reused with the other."""
    from cenrep.pipeline import Pipeline
    cdir = _complex(tmp_path, ["a"])
    a = Pipeline(cdir, 1, tmp_path / "cache", Params())
    b = Pipeline(cdir, 1, tmp_path / "cache", Params(seed_engine="trf"))
    assert all(a.fingerprint(s) != b.fingerprint(s) for s in STAGES)


def test_default_cache_keeps_using_an_older_runs_trf_cache(tmp_path):
    """Runs of cenrep < 0.9 cached TRF in OUTDIR/work/trf_cache (hours per large genome); a rerun must find it."""
    from cenrep.pipeline import Pipeline
    cdir = _complex(tmp_path, ["a"])
    assert Pipeline(cdir, 1, None, Params()).cache == cdir / "centromere" / "work" / "seed_cache"
    (cdir / "centromere" / "work" / "trf_cache").mkdir()
    assert Pipeline(cdir, 1, None, Params(seed_engine="trf")).cache == cdir / "centromere" / "work" / "trf_cache"


def test_unknown_seed_engine_fails_clearly(tmp_path):
    from cenrep.pipeline import Pipeline
    with pytest.raises(CenrepError, match="seed engine 'nope'"):
        Pipeline(_complex(tmp_path, ["a"]), 1, tmp_path / "cache", Params(seed_engine="nope"))


def test_classify_from_changes_classify_fingerprint_only(tmp_path):
    from cenrep.pipeline import Pipeline
    cdir = _complex(tmp_path, ["Ctg", "a"])
    a = Pipeline(cdir, 1, tmp_path / "cache", Params())
    b = Pipeline(cdir, 1, tmp_path / "cache", Params(classify_from="Ctg=a"))
    assert a.fingerprint("annotate") == b.fingerprint("annotate")
    assert a.fingerprint("classify") != b.fingerprint("classify")


def test_many_sequences_warn_unless_borrowed(tmp_path, caplog, monkeypatch):
    from cenrep.pipeline import Pipeline
    from cenrep.util import log
    monkeypatch.setattr(log, "propagate", True)
    cdir = _complex(tmp_path, ["Ctg", "a"], nseq=60)
    Pipeline(cdir, 1, tmp_path / "cache", Params(classify_from="Ctg=a"))
    assert [r for r in caplog.records if "60 sequences" in r.message and r.message.startswith("a:")]
    assert not [r for r in caplog.records if r.message.startswith("Ctg:")]
