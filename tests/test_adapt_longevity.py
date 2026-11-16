"""H26's trait (adaptations/longevity.py): the residual of lifespan on mass, its cuts, and what reaches H21's code.
No test prints an estimate of the analysis: these run before the pre-registered run does."""
import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.adaptations import clades, longevity, pairs
from mitoage_pipelines.adaptations.traits import COLUMNS

need = pytest.mark.skipif(not all(p.is_file() for p in (paths.BUILD / "species.parquet", paths.DATED_TREE)),
                          reason="needs `download build` and `longevity dated-tree`")


def made_up(n_per_order=(8, 6, 3), cls="Mammalia", seed=1) -> pd.DataFrame:
    """A species table with the columns the trait reads; lifespan follows mass with noise and an order effect."""
    rng = np.random.default_rng(seed)
    rows = []
    for k, n in enumerate(n_per_order):
        for i in range(n):
            lm = rng.uniform(1, 6)
            ll = 0.2 * lm + 0.3 * k + rng.normal(0, 0.15)
            rows.append({"slug": f"o{k}-s{i}", "class": cls, "order": f"Order{k}", "family": f"Fam{k}",
                         "adult_weight_g": 10 ** lm, "max_longevity_years": 10 ** ll, "record_status": "original",
                         "lifespan_quality": "acceptable", "lifespan_record_type": "captivity"})
    return pd.DataFrame(rows)


def test_state_follows_the_two_cuts_and_leaves_the_middle_and_the_missing_unknown():
    r = pd.Series([0.3, -0.2, np.nan, 0.05, 0.9, -0.5, 0.1, 0.0, np.nan, 0.4])
    state, how = longevity.state_of(r)
    known = r.dropna()
    hi, mid = known.quantile(0.75), known.quantile(0.5)
    assert set(state[state.eq(1).fillna(False)].index) == set(known[known >= hi].index)
    assert set(state[state.eq(0).fillna(False)].index) == set(known[known < mid].index)
    middle = known[(known >= mid) & (known < hi)].index
    assert len(middle) > 0 and state[middle].isna().all()
    assert how[middle].str.contains("intermediate").all()
    assert state[r.isna()].isna().all() and (how[r.isna()] == "no residual").all()
    # a value exactly on a cut: the upper cut is 1, the median is not 0
    assert state[known[known == hi].index].eq(1).all()
    assert not state[known[known == mid].index].eq(0).fillna(False).any()


def test_state_of_nothing_is_all_unknown():
    state, how = longevity.state_of(pd.Series([np.nan, np.nan]))
    assert state.isna().all() and (how == "no residual").all()


def test_residual_is_what_mass_does_not_explain():
    sp = made_up()
    r = longevity.residuals(sp, set(sp.slug), "Mammalia")
    assert set(r.slug) == set(sp.slug)
    assert abs(r.residual.sum()) < 1e-9
    assert abs(np.corrcoef(r.residual, r.log10_mass)[0, 1]) < 1e-9
    slope, intercept = np.polyfit(r.log10_mass, r.log10_lifespan, 1)
    assert np.allclose(r.residual, r.log10_lifespan - (intercept + slope * r.log10_mass))


def test_residual_in_order_is_centred_in_each_order_and_absent_in_small_orders():
    sp = made_up()
    r = longevity.residuals(sp, set(sp.slug), "Mammalia")
    counts = r.order.value_counts()
    small = r.order.map(counts) < longevity.MIN_ORDER
    assert small.any() and r.residual_in_order[small].isna().all()
    assert r.residual_in_order[~small].notna().all()
    assert (r[~small].groupby("order").residual_in_order.sum().abs() < 1e-9).all()
    big = r[~small]
    centred_mass = big.log10_mass - big.groupby("order").log10_mass.transform("mean")
    assert abs((big.residual_in_order * centred_mass).sum()) < 1e-9
    # the order effect is in the class-wide residual and not in the one within the order
    assert big.groupby("order").residual.mean().abs().max() > 0.05


def test_only_species_with_every_ingredient_get_a_residual():
    sp = made_up()
    sp.loc[0, "adult_weight_g"] = 0
    sp.loc[1, "adult_weight_g"] = np.nan
    sp.loc[2, "lifespan_quality"] = "low"
    sp.loc[3, "lifespan_record_type"] = "estimated"
    sp.loc[4, "record_status"] = "retired"
    usable = set(sp.slug) - {sp.slug[5]}
    other = made_up((4,), cls="Aves", seed=2).assign(slug=lambda d: "bird-" + d.slug)
    both = pd.concat([sp, other], ignore_index=True)
    r = longevity.residuals(both, usable | set(other.slug), "Mammalia")
    assert set(r.slug) == set(sp.slug[6:])


def test_a_species_without_a_residual_is_unknown_never_zero():
    sp = made_up()
    sp.loc[0, "adult_weight_g"] = np.nan
    sp.loc[1, "record_status"] = "retired"
    birds = made_up((9,), cls="Aves", seed=3).assign(slug=lambda d: "bird-" + d.slug)
    both = pd.concat([sp, birds], ignore_index=True)
    t, r = longevity.classify(both, set(both.slug))
    assert list(t.columns) == COLUMNS
    live = both[both.record_status.ne("retired")]
    for trait in (longevity.TRAIT, longevity.TRAIT_IN_ORDER):
        rows = t[t.trait.eq(trait)]
        assert sorted(rows.slug) == sorted(live.slug)          # one row per live species, the retired one has none
    row = t[(t.slug == sp.slug[0]) & t.trait.eq(longevity.TRAIT)].iloc[0]
    assert pd.isna(row.state) and "adult mass" in row.how
    # each class is cut on its own residuals
    for cls in longevity.CLASSES:
        rc = r[r["class"].eq(cls)].set_index("slug").residual
        tc = t[t.trait.eq(longevity.TRAIT) & t["class"].eq(cls)].set_index("slug").state
        assert set(tc[tc.eq(1).fillna(False)].index) == set(rc[rc >= rc.quantile(0.75)].index)
        assert set(tc[tc.eq(0).fillna(False)].index) == set(rc[rc < rc.quantile(0.5)].index)
    small = t[t.trait.eq(longevity.TRAIT_IN_ORDER) & t.order.eq("Order2")]
    assert small.state.isna().all() and small.how.str.contains("fewer than").all()


@pytest.fixture(scope="module")
def build():
    species, tips, dist = longevity.load(paths.BUILD, paths.DATED_TREE)
    t, r = longevity.classify(species, tips & dist.aligned)
    return species, tips, dist, t, r


@need
def test_on_the_build_the_states_are_the_cuts_of_the_residuals(build):
    species, tips, dist, t, r = build
    assert set(r["class"]) == set(longevity.CLASSES)
    sp = species.set_index("slug")
    assert set(r.slug) <= tips & dist.aligned
    assert sp.loc[r.slug, "adult_weight_g"].gt(0).all() and sp.loc[r.slug, "lifespan_quality"].isin(["high", "acceptable"]).all()
    for cls in longevity.CLASSES:
        rc = r[r["class"].eq(cls)].set_index("slug").residual
        tc = t[t.trait.eq(longevity.TRAIT) & t["class"].eq(cls)].set_index("slug").state
        assert len(rc) >= 40
        assert set(tc[tc.notna()].index) <= set(rc.index)
        assert set(tc[tc.eq(1).fillna(False)].index) == set(rc[rc >= rc.quantile(0.75)].index)
        assert set(tc[tc.eq(0).fillna(False)].index) == set(rc[rc < rc.quantile(0.5)].index)
        assert sorted(tc.index) == sorted(species.slug[species["class"].eq(cls)])


@need
def test_on_the_build_trait_clades_and_pairs_hold_only_their_class_and_state(build):
    species, tips, dist, t, r = build
    cls_of = species.set_index("slug")["class"]
    for cls in longevity.CLASSES:
        scope = species.slug[species["class"].eq(cls)]
        k = t[t.trait.eq(longevity.TRAIT) & t.state.notna() & t.slug.isin(set(scope))]
        state = {s: int(v) for s, v in zip(k.slug, k.state)}
        root = clades.prune(clades.read(paths.DATED_TREE), set(state))
        tc = clades.trait_clades(root, state)
        in_clades = [s for c in tc for s in clades.tips(c)]
        assert sorted(in_clades) == sorted(s for s, v in state.items() if v == 1)      # every 1 in exactly one clade
        p = pairs.run(paths.DATED_TREE, t, longevity.TRAIT, tips & dist.aligned, dist, scope)
        assert len(p) > 0
        assert len(p) == len(tc) and p.sp_trait.is_unique                                # one row per origin
        paired = p[p.sp_other.notna()]
        assert all(state[s] == 1 for s in p.sp_trait) and all(state[s] == 0 for s in paired.sp_other)
        assert set(cls_of[list(p.sp_trait) + list(paired.sp_other)]) == {cls}
        taken = p[p.taken]
        assert len(taken) > 0 and taken.sp_other.is_unique      # pairs taken share no branch, so no relative either


@need
def test_selection_prep_writes_one_class_and_the_trait_clades_first(build, tmp_path, monkeypatch):
    species, tips, dist, t, r = build
    t.to_csv(tmp_path / "traits.tsv", sep="\t", index=False)
    res = longevity.selection_prep(paths.BUILD, paths.DATED_TREE, tmp_path, "Aves")
    assert list(res.draw) == list(range(len(res))) and len(res) == 21
    d0 = tmp_path / "Aves" / "selection" / longevity.TRAIT / "draw00"
    got = set((d0 / "tips.txt").read_text().split())
    birds = set(species.slug[species["class"].eq("Aves")])
    assert got <= birds
    ones = set(t.slug[t.trait.eq(longevity.TRAIT) & t.state.eq(1).fillna(False) & t["class"].eq("Aves")])
    assert ones & got == ones & (tips & dist.aligned) and len(ones & got) > 0
    with pytest.raises(SystemExit):
        longevity.selection_prep(paths.BUILD, paths.DATED_TREE, tmp_path, "Reptilia")
