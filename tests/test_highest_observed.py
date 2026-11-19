"""H27's second lifespan: which observation replaces the chosen value, and that nothing else of a species moves."""
import pandas as pd
import pytest

from mitoage_pipelines.longevity import highest_observed as ho


def species() -> pd.DataFrame:
    return pd.DataFrame({
        "id": [1, 2, 3, 4, 5, 6], "slug": list("abcdef"), "class": ["Aves", "Aves", "Mammalia", "Mammalia", "Aves", "Aves"],
        "max_longevity_years": [10.0, 10.0, 10.0, 10.0, 10.0, 10.0], "lifespan_source": ["AnAge"] * 6,
        "lifespan_source_record_id": ["h1", "h2", "h3", "h4", "h5", "h6"], "lifespan_record_type": ["captive"] * 6,
        "lifespan_quality": ["acceptable"] * 6, "lifespan_sample_size": ["medium"] * 6, "adult_weight_g": [1.0] * 6,
    })


def obs(rows) -> pd.DataFrame:
    chosen = [(i, "AnAge", f"h{i}", 10.0, "captive", "acceptable", True) for i in range(1, 7)]
    return pd.DataFrame(chosen + rows, columns=["species_id", "source", "source_record_id", "value_years", "record_type", "quality", "is_chosen"])


def test_the_highest_qualifying_observation_replaces_the_chosen_value():
    sp = species()
    out, changed = ho.highest(sp, obs([
        (1, "EURING", "e1", 15.0, "wild", "acceptable", False),         # higher, observed, acceptable: taken
        (1, "USGS_BBL", "u1", 12.0, "wild", "acceptable", False),       # higher but not the highest
        (2, "FishBase", "f2", 40.0, "estimated", "low", False),         # an estimate: never taken
        (2, "FishBase", "f2b", 41.0, "estimated", "acceptable", False), # nor an estimate its source rates acceptable
        (3, "FishBase", "f3", 40.0, "observed", "questionable", False), # questionable: never taken
        (4, "MitoAge_build1", "m4", 40.0, "compiled", None, False),     # the 2015 database's own value: never taken
        (4, "MitoAge_build1", "m4b", 41.0, "compiled", "acceptable", False),
        (5, "EURING", "e5", 10.0, "wild", "acceptable", False),         # a tie keeps the chosen value
        (6, "EURING", "e6", 9.0, "wild", "high", False),                # lower: never taken
    ]))
    pd.testing.assert_frame_equal(sp, species())                          # the table handed in is not written to
    assert out.max_longevity_years.tolist() == [15.0, 10.0, 10.0, 10.0, 10.0, 10.0]
    assert out.lifespan_source.tolist() == ["EURING", "AnAge", "AnAge", "AnAge", "AnAge", "AnAge"]
    row = out.iloc[0]
    assert (row.lifespan_source_record_id, row.lifespan_record_type, row.lifespan_quality) == ("e1", "wild", "acceptable")
    assert pd.isna(row.lifespan_sample_size) and out.lifespan_sample_size[1:].tolist() == ["medium"] * 5
    assert changed.to_dict("records") == [{"slug": "a", "class": "Aves", "chosen_years": 10.0, "chosen_source": "AnAge",
                                           "highest_years": 15.0, "highest_source": "EURING", "ratio": 1.5}]
    # nothing but the lifespan columns moved, and no species was lost or reordered
    same = [c for c in species().columns if c not in list(ho.COLUMNS.values()) + ["lifespan_sample_size"]]
    pd.testing.assert_frame_equal(out[same], species()[same])


def test_no_lifespan_gets_lower_and_no_observation_changes_nothing():
    sp = species()
    out, changed = ho.highest(sp, obs([]))
    pd.testing.assert_frame_equal(out, sp)
    assert changed.empty
    out, _ = ho.highest(sp, obs([(i, "EURING", f"e{i}", v, "wild", "high", False) for i, v in zip(range(1, 7), [1, 5, 9.99, 10, 10.01, 99])]))
    assert (out.max_longevity_years >= sp.max_longevity_years).all()
    assert out.max_longevity_years.tolist() == [10, 10, 10, 10, 10.01, 99]


def test_main_writes_a_copy_and_never_over_the_build(tmp_path):
    build = tmp_path / "build"
    (build / "alignments").mkdir(parents=True)
    (build / "alignments" / "x.txt").write_text("x")
    species().to_parquet(build / "species.parquet")
    obs([(1, "EURING", "e1", 15.0, "wild", "acceptable", False)]).to_parquet(build / "lifespan_observations.parquet")
    with pytest.raises(SystemExit, match="must not be written over"):
        ho.main(build, build, tmp_path / "t.tsv")
    ho.main(build, tmp_path / "copy", tmp_path / "t.tsv")
    assert pd.read_parquet(build / "species.parquet").max_longevity_years.tolist() == [10.0] * 6        # the original is untouched
    assert pd.read_parquet(tmp_path / "copy" / "species.parquet").max_longevity_years.tolist() == [15.0] + [10.0] * 5
    assert (tmp_path / "copy" / "alignments" / "x.txt").read_text() == "x"
    assert (tmp_path / "copy" / "lifespan_observations.parquet").is_symlink()
    assert pd.read_csv(tmp_path / "t.tsv", sep="\t").slug.tolist() == ["a"]
