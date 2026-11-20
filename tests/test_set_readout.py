"""H27's readout: the reading of one result, and that rows are matched by their keys and never by their order."""
import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines.longevity import set_readout as sr


def test_reading_of_a_result_that_passed():
    r = lambda est, se, q: sr.reading(1.0, 0.01, est, se, q)
    assert r(0.9, 0.2, 0.01) == "holds"
    assert r(-0.9, 0.2, 0.01) == "reverses"
    assert r(0.8, 0.3, 0.20) == "holds, less precise"          # 0.8 +- 0.59 contains 1.0
    assert r(0.2, 0.3, 0.20) == "moves"                        # 0.2 +- 0.59 does not
    assert r(-0.2, 0.3, 0.20) == "moves"
    assert r(-0.1, 0.9, 0.80) == "unclear"                     # the other sign, and an interval that contains 1.0
    assert r(np.nan, 0.3, 0.20) == "not assessed"
    assert r(0.8, 0.3, np.nan) == "not assessed"
    assert r(0.8, np.nan, 0.20) == "not assessed"
    assert r(0.9, 0.2, 0.05) == "holds, less precise"          # the threshold is q below 0.05, not at it
    # the edge of the interval counts as inside
    assert sr.reading(1.0, 0.01, 0.5, 0.5 / sr.Z, 0.5) == "holds, less precise"
    assert sr.reading(1.0, 0.01, 0.5, 0.5 / sr.Z - 1e-9, 0.5) == "moves"


def test_reading_of_a_result_that_did_not_pass():
    assert sr.reading(1.0, 0.30, 0.9, 0.2, 0.01) == "new on the restricted set"
    assert sr.reading(1.0, 0.30, 0.9, 0.2, 0.40) == "not passed on either"
    assert sr.reading(1.0, np.nan, 0.9, 0.2, 0.01) == "new on the restricted set"
    assert sr.reading(1.0, 0.05, 0.9, 0.2, 0.40) == "not passed on either"


def write(root, rows):
    f = root / "output" / "longevity" / "rates.tsv"
    f.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["class", "gene", "normalisation", "includes_saturated", "n", "slope", "se", "q"]).to_csv(f, sep="\t", index=False)


def test_rows_are_matched_by_key_and_a_missing_row_is_not_assessed(tmp_path, monkeypatch):
    monkeypatch.setattr(sr, "TABLES", {"rates.tsv": sr.TABLES["rates.tsv"]})
    write(tmp_path / "all", [("Aves", "ND1", "dS_raw", True, 90, 1.0, 0.1, 0.01), ("Aves", "ND2", "dS_raw", True, 90, -2.0, 0.1, 0.01),
                             ("Aves", "ND3", "dS_raw", True, 90, 0.1, 0.1, 0.50)])
    write(tmp_path / "high", [("Aves", "ND2", "dS_raw", True, 70, -1.9, 0.1, 0.02), ("Aves", "ND1", "dS_raw", True, 70, 0.1, 0.1, 0.60)])  # other order, ND3 absent
    d = sr.main(tmp_path, tmp_path / "out", ["high", "wild"]).set_index(["set", "result"])
    assert d.loc[("high", "Aves | ND1 | dS_raw | True"), "reading"] == "moves"
    assert d.loc[("high", "Aves | ND2 | dS_raw | True"), "reading"] == "holds"
    assert d.loc[("high", "Aves | ND2 | dS_raw | True"), "n"] == 70 and d.loc[("high", "Aves | ND2 | dS_raw | True"), "n_all"] == 90
    assert d.loc[("high", "Aves | ND3 | dS_raw | True"), "reading"] == "not assessed"
    assert (d.loc["wild"].reading == "not assessed").all() and len(d.loc["wild"]) == 3        # a set with no table at all
    assert bool(d.loc[("high", "Aves | ND1 | dS_raw | True"), "interval_excludes_original"]) is True
    assert bool(d.loc[("high", "Aves | ND2 | dS_raw | True"), "interval_excludes_original"]) is False
    s = pd.read_csv(tmp_path / "out" / "summary.tsv", sep="\t").set_index("set")
    assert s.loc["high", ["holds", "moves", "not assessed"]].tolist() == [1, 1, 1]


def test_refusals(tmp_path, monkeypatch):
    monkeypatch.setattr(sr, "TABLES", {"rates.tsv": sr.TABLES["rates.tsv"]})
    with pytest.raises(SystemExit, match="unrestricted run has no"):
        sr.compare(tmp_path, ["high"])
    write(tmp_path / "all", [("Aves", "ND1", "dS_raw", True, 90, 1.0, 0.1, 0.01)] * 2)
    with pytest.raises(SystemExit, match="do not name one row each"):
        sr.compare(tmp_path, ["high"])
