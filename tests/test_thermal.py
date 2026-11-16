"""H20: the traits, FishBase matching, sister pairs, fits and readout of `longevity thermal`, on the downloaded
build, its dated tree and the FishBase tables of `download fishbase-temperature` (MITOAGE_FISHBASE overrides
where they are)."""
import math
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
from Bio import Phylo
from statsmodels.stats.multitest import multipletests

from mitoage_pipelines import paths
from mitoage_pipelines.longevity import thermal as th
from mitoage_pipelines.longevity.rates import SATURATED_DS, distance

FB = Path(os.environ.get("MITOAGE_FISHBASE", paths.B2_FISHBASE / "v26.06"))
FOURFOLD = paths.LONGEVITY_OUT / "fourfold.tsv"
pytestmark = pytest.mark.skipif(
    not all(p.is_file() for p in (paths.BUILD / "species.parquet", paths.DATED_TREE, FOURFOLD, FB / "stocks.parquet")),
    reason="needs `download build`, `longevity dated-tree`, `longevity fourfold` and `download fishbase-temperature`")


@pytest.fixture(scope="module")
def species():
    sp = pd.read_parquet(paths.BUILD / "species.parquet")
    return sp[sp.record_status.ne("retired")].reset_index(drop=True)


@pytest.fixture(scope="module")
def traits(species):
    return th.traits(species, FB)


@pytest.fixture(scope="module")
def dist(species):
    return th.Distances(paths.BUILD, species)


@pytest.fixture(scope="module")
def tree():
    return Phylo.read(paths.DATED_TREE, "newick")


@pytest.fixture(scope="module")
def fitted(species, traits, dist):
    return th.fit_pairs(species, traits[0], dist, paths.DATED_TREE)


def fish_members(traits, dist, tree, cls="Actinopterygii"):
    tr = traits[0]
    tips = {t.name for t in tree.get_terminals()}
    return set(tr[tr["class"].eq(cls) & tr.inv_kt.notna() & tr.slug.isin(tips & dist.aligned)].slug)


class TestTemperatures:
    def test_inv_kt_is_one_over_boltzmann_times_kelvin(self, traits):
        tr = traits[0].dropna(subset=["inv_kt"])
        assert len(tr) > 0
        assert np.allclose(tr.inv_kt * th.K_EV * tr.temp_k, 1.0)

    def test_a_body_temperature_out_of_range_is_unknown(self, species):
        sp = species.copy()
        mam = sp.index[sp["class"].eq("Mammalia") & sp.body_temperature_k.notna()]
        sp.loc[mam[0], "body_temperature_k"] = 400.0
        sp.loc[mam[1], "body_temperature_k"] = 250.0
        tr, _ = th.traits(sp, FB)
        assert tr.loc[mam[:2], "inv_kt"].isna().all()
        assert tr.loc[mam[2:], "inv_kt"].notna().all()

    def test_fish_temperature_is_the_median_of_its_stocks(self, traits):
        tr = traits[0]
        st = pd.read_parquet(FB / "stocks.parquet", columns=["SpecCode", "TempPref50"])
        st = st[st.TempPref50.between(*th.WATER_C)]
        want = st.groupby("SpecCode").TempPref50.median()
        fish = tr[tr.fishbase_speccode.notna()]
        got = fish.temp_k - 273.15
        exp = fish.fishbase_speccode.astype(int).map(want)
        assert got.notna().sum() == exp.notna().sum() > 100
        assert np.allclose(got.dropna(), exp.dropna())

    def _copy(self, st, tmp_path):
        for t in ("species", "synonyms"):
            shutil.copy(FB / f"{t}.parquet", tmp_path / f"{t}.parquet")
        st.to_parquet(tmp_path / "stocks.parquet")
        return th.fishbase_temperatures(tmp_path)

    def _code(self, traits, st):
        # in FishBase v26.06 no species has two stocks with a preferred temperature, so the tests add stocks
        codes = set(traits[0].fishbase_speccode.dropna().astype(int))
        return int(st[st.SpecCode.isin(codes) & st.TempPref50.notna()].SpecCode.iloc[0])

    def test_several_stocks_give_their_median(self, traits, tmp_path):
        st = pd.read_parquet(FB / "stocks.parquet")
        code = self._code(traits, st)
        v = float(st[st.SpecCode.eq(code)].TempPref50.dropna().iloc[0])
        extra = st[st.SpecCode.eq(code) & st.TempPref50.notna()].iloc[[0, 0]].copy()
        extra["TempPref50"] = [v + 1, v + 10]
        assert self._copy(pd.concat([st, extra], ignore_index=True), tmp_path).loc[code, "tpref"] == pytest.approx(v + 1)

    def test_a_stock_temperature_out_of_range_is_not_used(self, traits, tmp_path):
        st = pd.read_parquet(FB / "stocks.parquet")
        code = self._code(traits, st)
        st.loc[st.SpecCode.eq(code) & st.TempPref50.notna(), "TempPref50"] = 105.0
        assert np.isnan(self._copy(st, tmp_path).loc[code, "tpref"])

    def test_fish_size_is_log_length_and_others_log_mass(self, traits, species):
        tr = traits[0].set_index("slug")
        fish = tr["class"].isin(th.FISHES)
        assert np.allclose(tr.loc[fish, "size"].dropna(), tr.loc[fish, "llen"].dropna())
        other = ~fish & tr.lm.notna()
        assert np.allclose(tr.loc[other, "size"], tr.loc[other, "lm"])


class TestFishBaseNames:
    def test_an_exact_match_names_that_species(self, traits):
        match = traits[1]
        sp = pd.read_parquet(FB / "species.parquet", columns=["SpecCode", "Genus", "Species"]).set_index("SpecCode")
        build = pd.read_parquet(paths.BUILD / "species.parquet", columns=["slug", "scientific_name"]).set_index("slug")
        ex = match[match.how.eq("exact name")]
        assert len(ex) > 100
        named = (sp.Genus + " " + sp.Species).str.strip()
        assert (ex.SpecCode.astype(int).map(named).to_numpy() == ex.slug.map(build.scientific_name).str.strip().to_numpy()).all()

    def _synonyms(self):
        sp = pd.read_parquet(FB / "species.parquet", columns=["Genus", "Species"])
        syn = pd.read_parquet(FB / "synonyms.parquet", columns=["SpecCode", "SynGenus", "SynSpecies", "Status"])
        syn = syn.assign(name=(syn.SynGenus + " " + syn.SynSpecies).str.strip())
        return syn[~syn.name.isin(set((sp.Genus + " " + sp.Species).str.strip()))]

    def _one(self, name):
        return pd.DataFrame({"slug": ["x"], "scientific_name": [name], "lifespan_source": ["AnAge"],
                             "lifespan_source_record_id": [""]})

    def test_a_synonym_of_two_species_stays_unmatched(self):
        syn = self._synonyms()
        syn = syn[syn.Status.isin(th.SYNONYM_STATUS)]
        n = syn.groupby("name").SpecCode.nunique()
        name = n[n > 1].index[0]
        got = th.fishbase_codes(FB, self._one(name)).iloc[0]
        assert got.how == "ambiguous synonym" and pd.isna(got.SpecCode)

    def test_a_misapplied_name_is_no_match(self):
        syn = self._synonyms()
        good = set(syn[syn.Status.isin(th.SYNONYM_STATUS)].name)
        bad = syn[syn.Status.eq("misapplied name") & ~syn.name.isin(good)]
        got = th.fishbase_codes(FB, self._one(bad.name.iloc[0])).iloc[0]
        assert got.how == "no match" and pd.isna(got.SpecCode)

    def test_a_match_that_disagrees_with_the_builds_record_stays_unmatched(self, species, traits):
        match = traits[1].set_index("slug")
        fish = species[species.lifespan_source.eq("FishBase")
                       & species.lifespan_source_record_id.astype(str).str.match(th.BUILD_FB_ID.pattern)]
        fish = fish[fish.slug.map(match.how).eq("exact name")]
        assert len(fish) > 10
        row = fish.iloc[[0]].copy()
        code = int(th.BUILD_FB_ID.match(row.lifespan_source_record_id.iloc[0]).group(1))
        assert th.fishbase_codes(FB, row).iloc[0].SpecCode == code
        row["lifespan_source_record_id"] = f"species.LongevityWild:{code + 1}"
        got = th.fishbase_codes(FB, row).iloc[0]
        assert got.how == "conflicts with the build's FishBase record" and pd.isna(got.SpecCode)


class TestPairs:
    def test_pairs_share_no_branch_and_no_species(self, traits, dist, tree):
        members = fish_members(traits, dist, tree)
        pairs = th.disjoint_pairs(tree, members)
        names = [s for a, b, _ in pairs for s in (a, b)]
        assert len(names) == len(set(names)) and set(names) <= members
        parents = {c: p for p in tree.find_clades() for c in p.clades}
        leaf = {t.name: t for t in tree.get_terminals()}
        used = set()
        for a, b, age in pairs:
            mrca = tree.common_ancestor(leaf[a], leaf[b])
            edges = set()
            for s in (a, b):
                node = leaf[s]
                while node is not mrca:
                    edges.add(id(node))
                    node = parents[node]
            assert not edges & used
            used |= edges
            assert age == pytest.approx((tree.distance(mrca, leaf[a]) + tree.distance(mrca, leaf[b])) / 2)

    def test_no_two_unpaired_members_could_still_be_paired(self, traits, dist, tree):
        members = fish_members(traits, dist, tree)
        pairs = th.disjoint_pairs(tree, members)
        paired = {s for a, b, _ in pairs for s in (a, b)}
        left = members - paired
        # under every node at most one member goes up unpaired, so the whole tree leaves at most one over
        assert len(left) <= 1

    def test_the_rate_is_the_distance_over_twice_the_age(self, traits, dist, tree, species):
        members = fish_members(traits, dist, tree)
        pairs = th.disjoint_pairs(tree, members)[:12]
        p = th.pair_table(pairs, dist, traits[0], ["inv_kt"])
        assert p.dS_per_myr.notna().sum() >= 3
        aln = dist.aln
        for (a, b, age), row in zip(pairs, p.itertuples()):
            sa, sb = ("".join(aln[g][s] for g in th.MT_GENES) for s in (a, b))
            d = distance(sa, sb, dist.code(a, b))
            assert row.dS_per_myr == pytest.approx(d["dS"] / (2 * age), nan_ok=True)
            assert row.dN_per_myr == pytest.approx(d["dN"] / (2 * age), nan_ok=True)

    def test_a_pair_trait_is_missing_when_a_sister_lacks_it(self, traits, dist):
        tr = traits[0]
        mam = tr[tr["class"].eq("Mammalia") & tr.slug.isin(dist.aligned)]
        a, b = mam[mam.inv_kt.notna()].slug.iloc[0], mam[mam.inv_kt.isna()].slug.iloc[0]
        p = th.pair_table([(a, b, 10.0)], dist, tr, ["inv_kt", "size"])
        assert np.isnan(p.inv_kt.iloc[0])
        t = tr.set_index("slug")
        assert p["size"].iloc[0] == pytest.approx((t.at[a, "size"] + t.at[b, "size"]) / 2)

    def test_the_cherries_are_the_plain_cherries_of_the_tree(self, fitted, tree, traits, dist):
        _, pairs = fitted
        ch = pairs[pairs["class"].eq("Actinopterygii") & pairs.model.eq("temperature") & pairs.pairing.eq("cherries")]
        members = fish_members(traits, dist, tree)
        want = {(a, b) for _, a, b in th.cherries(tree, members)}
        assert set(zip(ch.sp_a, ch.sp_b)) == want and len(want) > 0
        dis = pairs[pairs["class"].eq("Actinopterygii") & pairs.model.eq("temperature") & pairs.pairing.eq("disjoint")]
        assert {frozenset(x) for x in want} <= {frozenset(x) for x in zip(dis.sp_a, dis.sp_b)}

    def test_saturated_pairs_are_not_used(self, fitted):
        _, pairs = fitted
        assert pairs.saturated.any()
        assert not (pairs.usable & (pairs.dS > SATURATED_DS)).any()
        fits = fitted[0]
        counted = pairs.groupby(["class", "model", "pairing"]).usable.sum()
        f = fits.drop_duplicates(["class", "model", "pairing"]).set_index(["class", "model", "pairing"]).n_pairs
        assert (f == counted.reindex(f.index)).all()


class TestFits:
    def test_lambda_zero_is_weighted_least_squares(self, fitted, tree):
        _, pairs = fitted
        d = pairs[pairs["class"].eq("Mammalia") & pairs.model.eq("size") & pairs.pairing.eq("disjoint")
                  & pairs.usable].reset_index(drop=True)
        C, order, _ = th.class_covariance(paths.DATED_TREE, list(d.sp_a))
        d = d.set_index("sp_a").loc[order].reset_index()
        y = np.log10(d.dS_per_myr.to_numpy())
        r = th.fit(C, y, d[["size"]].to_numpy(), ["size"], lam=0.0)
        w = sm.WLS(y, sm.add_constant(d["size"].to_numpy()), weights=1 / np.diag(C)).fit()
        assert r["size"][0] == pytest.approx(w.params[1])
        assert r["size"][3] == pytest.approx(w.pvalues[1])

    def test_E_is_minus_ln10_times_the_slope_and_inside_its_interval(self, fitted):
        f = fitted[0]
        e = f[f.term.str.startswith("inv_kt") & f.status.eq("assessed")]
        assert len(e) > 0
        assert np.allclose(e.E_eV, -math.log(10) * e.slope)
        assert ((e.E_low <= e.E_eV) & (e.E_eV <= e.E_high)).all()

    def test_holm_over_the_assessed_primary_tests(self, fitted):
        ro = th.readout(fitted[0], pd.DataFrame())
        prim = ro[ro.kind.eq("primary") & ro.status.eq("assessed")]
        assert len(prim) >= 2
        assert np.allclose(prim.p_holm, multipletests(prim.p, method="holm")[1])
        assert ro[ro.kind.eq("primary") & ~ro.status.eq("assessed")].reading.eq("not assessed").all()

    def test_a_failed_size_control_leaves_its_class_unread(self, fitted):
        f = fitted[0].copy()
        ctl = (f["class"].eq("Mammalia") & f.model.eq("size") & f.response.eq("log10_dS_per_myr")
               & f.pairing.eq("disjoint"))
        f.loc[ctl, "p"] = 0.9
        ro = th.readout(f, pd.DataFrame()).set_index("test")
        for t in ("T3", "M1"):
            if ro.at[t, "status"] == "assessed":
                assert ro.at[t, "reading"].startswith("not read")
        f.loc[ctl, ["p", "slope"]] = [1e-6, -0.1]
        ro = th.readout(f, pd.DataFrame()).set_index("test")
        assert not ro.loc[["T3", "M1"]].reading.str.startswith("not read").any()

    def test_lehmann_r2_is_that_of_ordinary_least_squares(self, species, traits):
        feats = th.features(paths.BUILD, FOURFOLD)
        leh = th.lehmann(species, traits[0], feats, paths.DATED_TREE)
        row = leh[leh.model.eq("gc+metabolic") & leh.term.eq(th.GC)].iloc[0]
        tr = traits[0].merge(feats[["slug", th.GC]], on="slug")
        tips = {t.name for t in Phylo.read(paths.DATED_TREE, "newick").get_terminals()}
        keep = set(species[species["class"].eq("Mammalia") & th.default_dataset(species)].slug) & tips
        d = tr[tr.slug.isin(keep)].dropna(subset=["ll", th.GC, "lmr"])
        ols = sm.OLS(d.ll, sm.add_constant(d[[th.GC, "lmr"]])).fit()
        assert row.n == len(d)
        assert row.ols_r2 == pytest.approx(ols.rsquared)
        assert row.ols_slope == pytest.approx(ols.params[th.GC])
