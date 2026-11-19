"""H27's species sets: who is in each, that a build without a confidence level is refused, and that the default
dataset of every lifespan analysis follows the set named in the environment and nothing else."""
import pandas as pd
import pytest

from mitoage_pipelines.longevity import aa_pgls, species_sets as ss

ROWS = [  # slug, level, source, sample, kind
    ("a", "high", "AnAge", "large", "captive"),
    ("b", "high", "AnAge", "small", "wild"),
    ("c", "caveat", "AnAge", "medium", "captive"),
    ("d", "caveat", "FishBase", None, "wild"),
    ("e", "low", "AnAge", "huge", "estimated"),
    ("f", "unreliable", "AnAge", "large", "captive"),
    ("g", None, "AnAge", "large", "wild"),
    ("h", "caveat", "FishBase", "medium", "wild"),          # a sample size from another source than AnAge
]


def table() -> pd.DataFrame:
    sp = pd.DataFrame(ROWS, columns=["slug", "data_confidence", "lifespan_source", "lifespan_sample_size", "lifespan_record_type"])
    sp["record_status"] = "new"
    sp["lifespan_quality"] = ["acceptable", "high", "acceptable", "acceptable", "low", "acceptable", "acceptable", "acceptable"]
    return sp


def in_set(name: str) -> list[str]:
    sp = table()
    return sp.slug[ss.members(sp, name)].tolist()


def test_members_of_each_set():
    assert in_set("not_unreliable") == ["a", "b", "c", "d", "e", "h"]   # not f (unreliable), not g (no level)
    assert in_set("high") == ["a", "b"]
    assert in_set("well_sampled") == ["a", "c", "e"]                    # AnAge and medium, large or huge; not h
    assert in_set("captive") == ["a", "c"]
    assert in_set("wild") == ["b", "d", "h"]
    for name in ss.SETS:                                                # every set lies within not_unreliable
        assert set(in_set(name)) <= set(in_set("not_unreliable"))


def test_a_missing_level_is_in_no_set_and_an_unknown_level_or_set_is_refused():
    for name in ss.SETS:
        assert "g" not in in_set(name) and "f" not in in_set(name)
    with pytest.raises(ValueError):
        ss.members(table(), "reliable")
    bad = table()
    bad.loc[0, "data_confidence"] = "medium"
    with pytest.raises(SystemExit, match="no level"):
        ss.members(bad, "high")


def test_a_build_without_the_level_is_refused():
    with pytest.raises(SystemExit, match="no data_confidence"):
        ss.members(table().drop(columns="data_confidence"), "not_unreliable")


def test_active_reads_the_environment(monkeypatch):
    monkeypatch.delenv(ss.ENV, raising=False)
    assert ss.active() is None
    monkeypatch.setenv(ss.ENV, " ")
    assert ss.active() is None
    monkeypatch.setenv(ss.ENV, "high")
    assert ss.active() == "high"
    monkeypatch.setenv(ss.ENV, "High")
    with pytest.raises(SystemExit, match="not a species set"):
        ss.active()


def test_default_dataset_follows_the_set_and_is_unchanged_without_one(tmp_path, monkeypatch):
    sp = table()
    sp.to_parquet(tmp_path / "species.parquet")
    monkeypatch.setattr(aa_pgls.paths, "BUILD", tmp_path)
    monkeypatch.delenv(ss.ENV, raising=False)
    plain = sp.slug[aa_pgls.default_dataset(sp)].tolist()
    assert plain == ["a", "b", "c", "d", "f", "g", "h"]                       # the entry's own rule: e is low and estimated
    for name in ss.SETS:
        monkeypatch.setenv(ss.ENV, name)
        got = sp.slug[aa_pgls.default_dataset(sp)].tolist()
        assert got == [s for s in plain if s in in_set(name)]           # the rule and the set, nothing else
    monkeypatch.setenv(ss.ENV, "high")
    with pytest.raises(SystemExit, match="without their slug"):
        aa_pgls.default_dataset(sp.drop(columns="slug"))
    # a frame that carries only the columns an analysis reads is restricted from the build's table, not from itself
    narrow = sp[["slug", "record_status", "lifespan_quality", "lifespan_record_type"]]
    assert narrow.slug[aa_pgls.default_dataset(narrow)].tolist() == ["a", "b"]
