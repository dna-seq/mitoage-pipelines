"""H3's two-step check (`longevity branch-rates`, `branch-pgls`) on real data.

branch_rates_fit12.json: HyPhy 2.5.104 with hyphy/local_mg94.bf on 12 of Coevol's 106 mammals (the first 12 of
mammals.ali, its first 300 codons) under the vertebrate mitochondrial code, run on 2026-11-15.
coevol_dated106.nwk: `coevol-prep --fixed-ages`'s tree, the dated tree pruned to Coevol's 106 species (Myr).
coevol_mammals.lht: Coevol's lifespan and mass for the same species."""
import json

import numpy as np
import pandas as pd
import pytest
from Bio import Phylo

from mitoage_pipelines.longevity import branch_rates as br

DATA = __import__("pathlib").Path(__file__).parent / "data"
FIT = DATA / "branch_rates_fit12.json"
DATED = DATA / "coevol_dated106.nwk"
LHT = DATA / "coevol_mammals.lht"


@pytest.fixture(scope="module")
def fitted():
    return br.table(FIT, DATED)


class TestTable:
    def test_branch_lengths_are_alpha_S_plus_beta_N(self, fitted):
        df, info = fitted
        assert len(df) == 21 and df.terminal.sum() == 12            # 12 tips, unrooted: 2n - 3 branches
        assert np.allclose(df.alpha * info["S"] + df.beta * info["N"], df.length, rtol=1e-9, atol=1e-12)
        assert 0.15 < info["synonymous_fraction"] < 0.4

    def test_dS_dN_and_their_ratio(self, fitted):
        df, info = fitted
        k = info["S"] + info["N"]
        assert np.allclose(df.dS, df.alpha * k) and np.allclose(df.dN, df.beta * k)
        ok = df.alpha > 0
        assert np.allclose(df.dNdS[ok], df.beta[ok] / df.alpha[ok])
        assert df.dNdS[~ok].isna().all() and (~ok).sum() == 2      # two short internal branches have alpha 0
        t = df[df.terminal]
        assert (t.dNdS < 1).all() and (t.dNdS > 0.01).all()       # mitochondrial proteins: strong purifying selection

    def test_durations_cover_the_dated_tree_once(self, fitted):
        """Every Myr of the dated tree pruned to the 12 species is on exactly one fitted branch: the two root
        branches HyPhy merged are summed onto the one it kept."""
        df, _ = fitted
        tree = Phylo.read(DATED, "newick")
        keep = set(df.species[df.terminal])
        for t in tree.get_terminals():
            if t.name not in keep:
                tree.prune(t)
        above_root = tree.root.branch_length or 0.0                 # Bio.Phylo's prune leaves one there
        assert df.myr.sum() == pytest.approx(tree.total_branch_length() - above_root, rel=1e-9)
        root_age = max(tree.depths().values())
        marsupial = df[df.species == "antechinus_flavipes"].iloc[0]
        assert marsupial.myr > root_age                             # its branch runs through the root
        a, b = df[df.species == "addax_nasomaculatus"].iloc[0], df[df.species == "antilocapra_americana"].iloc[0]
        assert a.myr == pytest.approx(b.myr, rel=1e-12)            # sisters: the same age
        assert np.allclose(df.dS_per_myr[df.terminal], df.dS[df.terminal] / df.myr[df.terminal])

    def test_refuses_lengths_of_another_model(self, tmp_path):
        d = json.loads(FIT.read_text())
        first = next(iter(d["branch lengths"]))
        d["branch lengths"][first]["MLE"] *= 1.5
        f = tmp_path / "fit.json"; f.write_text(json.dumps(d))
        with pytest.raises(ValueError):
            br.table(f, DATED)

    def test_refuses_species_missing_from_the_dated_tree(self, tmp_path):
        text = DATED.read_text().replace("acinonyx_jubatus", "felis_nonexistens")
        f = tmp_path / "dated.nwk"; f.write_text(text)
        with pytest.raises(ValueError, match="species are not in the dated tree"):
            br.table(FIT, f)


def tree_noise(rng, species, sd: float) -> np.ndarray:
    """Noise that evolved along the real dated tree (Brownian), one value per species in `species`' order."""
    tree = Phylo.read(DATED, "newick")
    order = {t.name: i for i, t in enumerate(tree.get_terminals())}
    C = br.shared_path_covariance(tree.root, order)
    C = C / np.diag(C).mean()
    z = rng.multivariate_normal(np.zeros(len(C)), C) * sd
    return z[[order[s] for s in species]]


def synthetic(seed: int, b_ll: float, b_lm: float, noise: float = 0.05, brownian: bool = False) -> pd.DataFrame:
    """Terminal branches of the 106 real species with real lifespan and mass; responses built with known slopes.
    The noise of dN/dS is independent per species, or with `brownian` evolved along the dated tree."""
    rng = np.random.default_rng(seed)
    tr = br.read_traits(LHT)
    e = tree_noise(rng, tr.species, noise) if brownian else rng.normal(0, noise, len(tr))
    dNdS = 10 ** (-1.2 + b_ll * tr.ll + b_lm * tr.lm + e)
    dSm = 10 ** (-1.8 - 0.1 * tr.lm + rng.normal(0, noise, len(tr)))
    return pd.DataFrame(dict(branch=tr.species, terminal=True, species=tr.species, dNdS=dNdS, dS_per_myr=dSm,
                             dS=np.where(np.arange(len(tr)) % 10 == 0, 2.5, 0.8)))


class TestPgls:
    def test_traits_are_the_real_106(self):
        tr = br.read_traits(LHT)
        assert len(tr) == 106 and tr.species.is_unique
        assert tr.ll.between(0, 2.5).all() and tr.lm.between(0, 8).all()
        g = tr[tr.species == "gorilla_gorilla"].iloc[0]
        assert g.ll == pytest.approx(np.log10(60.1)) and g.lm == pytest.approx(np.log10(139842.0))

    def test_recovers_known_slopes(self):
        res = br.pgls(synthetic(1, 0.0, 0.08), br.read_traits(LHT), DATED)
        r = lambda resp, model, term, b="all": res[(res.response == resp) & (res.model == model) & (res.term == term) & (res.branches == b)].iloc[0]
        assert r("log10_dNdS", "mass", "lm").slope == pytest.approx(0.08, abs=0.02)
        assert r("log10_dNdS", "lifespan+mass", "lm").slope == pytest.approx(0.08, abs=0.03)
        assert abs(r("log10_dNdS", "lifespan+mass", "ll").slope) < 0.1
        assert r("log10_dS_per_myr", "mass", "lm").slope == pytest.approx(-0.1, abs=0.02)
        assert (res.lam.between(0, 1)).all() and (res.n[res.branches == "all"] == 106).all()
        # the sensitivity rows leave out the branches with dS above 2 (every tenth here)
        assert (res.n[res.branches != "all"] == 106 - 11).all()

    def test_lambda_follows_the_noise(self):
        """Pagel's lambda near 0 for noise independent per species, near 1 for noise that evolved on the tree;
        either way the slope is recovered. A covariance whose rows are not the species' would not see it."""
        lam = lambda d: br.pgls(d, br.read_traits(LHT), DATED).query("response == 'log10_dNdS' and branches == 'all'").lam.iloc[0]
        assert lam(synthetic(7, 0.0, 0.08, noise=0.1)) < 0.3
        d = synthetic(8, 0.0, 0.08, noise=0.3, brownian=True)
        assert lam(d) > 0.7
        res = br.pgls(d, br.read_traits(LHT), DATED)
        s = res.query("response == 'log10_dNdS' and branches == 'all' and model == 'mass'").iloc[0]
        assert s.slope == pytest.approx(0.08, abs=2.5 * s.se)

    def test_holm_over_the_two_lifespan_tests_only(self):
        res = br.pgls(synthetic(2, 0.0, 0.08), br.read_traits(LHT), DATED)
        q = res.dropna(subset=["p_holm"])
        assert len(q) == 2 and set(q.response) == set(br.RESPONSES)
        assert (q.model == "lifespan+mass").all() and (q.term == "ll").all() and (q.branches == "all").all()
        assert np.allclose(q.p_holm, np.minimum(1, np.maximum.accumulate(np.sort(q.p.to_numpy()) * [2, 1]))[np.argsort(np.argsort(q.p.to_numpy()))])

    def test_a_branch_without_a_logarithm_is_left_out_and_counted(self):
        d = synthetic(3, 0.0, 0.08)
        d.loc[[4, 7], "dNdS"] = np.nan                 # alpha 0
        d.loc[9, "dNdS"] = 0.0                         # beta 0
        res = br.pgls(d, br.read_traits(LHT), DATED)
        w = res[(res.response == "log10_dNdS") & (res.branches == "all")]
        assert (w.n == 103).all() and (w.no_log == 3).all()

    def test_reading_follows_the_fixed_rules(self):
        res = br.pgls(synthetic(4, 0.25, 0.08), br.read_traits(LHT), DATED)
        lines = br.reading(res)
        assert "passes" in lines[0] and "agrees with Coevol" in lines[2]
        res = br.pgls(synthetic(5, -0.25, 0.08), br.read_traits(LHT), DATED)
        assert "DISAGREES" in br.reading(res)[2]
        res = br.pgls(synthetic(6, 0.0, -0.08), br.read_traits(LHT), DATED)
        lines = br.reading(res)
        assert "FAILS" in lines[0] and all("not read" in l for l in lines[2:])
