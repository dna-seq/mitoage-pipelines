"""H27's sample-size covariate (longevity/sample_size.py), on synthetic data whose answer is known."""
import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines.longevity import sample_size as ssz, species_sets


def test_sample_size_is_anages_only_and_ordered():
    sp = pd.DataFrame({"slug": list("abcdefg"),
                       "lifespan_source": ["AnAge"] * 5 + ["FishBase", "AnAge"],
                       "lifespan_sample_size": ["tiny", "small", "medium", "large", "huge", "large", None]})
    out = ssz.sample_size(sp)
    assert out.slug.tolist() == list("abcde")                    # another source's and a missing one are not in
    assert out.ss.tolist() == [0.0, 1.0, 2.0, 3.0, 4.0]


def test_an_unknown_sample_size_is_refused_not_dropped():
    sp = pd.DataFrame({"slug": ["a", "b"], "lifespan_source": ["AnAge", "AnAge"], "lifespan_sample_size": ["large", "enormous"]})
    with pytest.raises(SystemExit, match="enormous"):
        ssz.sample_size(sp)
    with pytest.raises(SystemExit, match="lifespan_source"):
        ssz.sample_size(sp.drop(columns="lifespan_source"))


def synthetic(n=400, seed=1):
    """A response made by sample size alone, and a lifespan that goes with sample size."""
    rng = np.random.default_rng(seed)
    blocks = rng.integers(0, 20, n)
    C = 0.5 * (blocks[:, None] == blocks[None, :]) + 0.5 * np.eye(n)   # species of one block share half their path
    ss = rng.integers(0, 5, n).astype(float)
    ll = 0.3 * ss + rng.normal(0, 0.3, n)
    mass = rng.normal(0, 1, n)
    y = 2.0 * ss + 0.5 * mass + rng.normal(0, 0.5, n)
    return C, y, np.column_stack([ll, mass]), ss


@pytest.mark.parametrize("lam", [1.0, None])
def test_the_term_takes_what_sample_size_made(lam):
    C, y, X, ss = synthetic()
    r = ssz.both(C, y, X, ss, lam)
    assert r["slope_without"] > 2 and r["p_without"] < 1e-10      # lifespan carries the sample size
    assert abs(r["slope_with"]) < 4 * r["se_with"] and abs(r["slope_with"]) < 0.3
    assert abs(r["ss_slope"] - 2.0) < 4 * r["ss_se"] and r["ss_p"] < 1e-10


def test_a_term_unrelated_to_the_response_leaves_the_slope():
    C, _, X, ss = synthetic(seed=2)
    rng = np.random.default_rng(3)
    ll = rng.normal(0, 1, len(ss))                               # lifespan independent of sample size
    y = 1.5 * ll + rng.normal(0, 0.5, len(ss))
    r = ssz.both(C, y, np.column_stack([ll, X[:, 1]]), ss, 1.0)
    assert abs(r["slope_with"] - r["slope_without"]) < 0.02
    assert abs(r["slope_with"] - 1.5) < 4 * r["se_with"]
    assert r["ss_p"] > 0.001


def test_reading_and_q_within_the_entries_groups():
    rows = []
    for model, p_with in (("mass", 1e-9), ("mass+composition", 0.6)):
        for k in range(20):
            hit = k == 0
            rows.append({"analysis": "H1", "class": "Mammalia", "method": "pgls", "competitor": "", "model": model, "feature": f"f{k}",
                         "n": 100, "slope_without": 1.0, "se_without": 0.1, "p_without": (1e-9 if model == "mass" else 1e-6) if hit else 0.5 + k / 100,
                         "slope_with": 1.0 if model == "mass" else 0.1, "se_with": 0.1,
                         "p_with": p_with if hit else 0.5 + k / 100, "ss_p": 0.01 if hit else 0.9})
    res = ssz.read(pd.DataFrame(rows))
    hit = res[res.feature == "f0"].set_index("model")
    assert hit.loc["mass", "q_without"] == pytest.approx(20 * 1e-9)   # 20 tests in the group, not 40
    assert hit.loc["mass", "reading"] == "holds"
    assert hit.loc["mass+composition", "reading"] == "moves"          # same sign, not passed, interval excludes 1.0
    s = ssz.summary(res).set_index("model")
    assert s.loc["mass", ["tests", "pass_without", "holds", "moves", "term_p_below_005"]].tolist() == [20, 1, 1, 0, 1]
    assert s.loc["mass+composition", ["pass_without", "holds", "moves"]].tolist() == [1, 0, 1]
    assert s.loc["mass+composition", "median_slope_ratio"] == pytest.approx(0.1)


def test_the_models_of_h11_all_hold_a_lifespan_term():
    m = ssz.h11_models()
    assert set(m) == {"maturity", "gestation", "birth_weight", "all"}
    assert all("ll" in cols for specs in m.values() for _, cols in specs)
    assert [name for name, _ in m["all"]] == ["lifespan", "all"]      # the two slow-fast axes are not refitted


def test_another_species_set_is_refused(monkeypatch, tmp_path):
    for value in ("", "high"):
        monkeypatch.setenv(species_sets.ENV, value)
        with pytest.raises(SystemExit, match="not_unreliable"):
            ssz.main(tmp_path, tmp_path / "c", tmp_path / "f", tmp_path / "t", tmp_path / "o")
