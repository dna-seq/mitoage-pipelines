"""H12: the features and the fits of `longevity structure-pgls`, on the downloaded build and the contact map."""
from collections import Counter

import numpy as np
import pandas as pd
import pytest
from statsmodels.stats.multitest import multipletests

from mitoage_pipelines import paths
from mitoage_pipelines.longevity import structure_pgls as sp_
from mitoage_pipelines.seqio import read_fasta

COMPOSITION = paths.LONGEVITY_OUT / "composition.tsv"
FOURFOLD = paths.LONGEVITY_OUT / "fourfold.tsv"
SITES = paths.MITOPORT_IFACE / "sites.tsv"
needs_build = pytest.mark.skipif(
    not (paths.BUILD / "species.parquet").is_file() or not COMPOSITION.is_file() or not FOURFOLD.is_file(),
    reason="needs `mitoage download build`, `longevity controls`, `longevity fourfold`, `longevity dated-tree`")
needs_sites = pytest.mark.skipif(not SITES.is_file(), reason="needs `download pdb 5Z62`, `mitoport contacts`, `mitoport burial`")
# Kyte & Doolittle 1982, the scale GRAVY averages
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8,
          P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)


@pytest.fixture(scope="module")
def birds():
    sp, _ = sp_.select(paths.BUILD, COMPOSITION, FOURFOLD, paths.DATED_TREE, "Aves")
    C, order, _ = sp_.class_covariance(paths.DATED_TREE, list(sp.slug))
    return sp.set_index("slug").loc[order].reset_index(), C


@pytest.fixture(scope="module")
def alns():
    return {g: read_fasta(paths.BUILD / "alignments" / f"{g}.aa.fasta.gz") for g in sp_.COX}


@needs_build
@needs_sites
class TestContactMismatch:
    def test_mismatch_is_counted_against_the_classes_own_consensus(self, birds, alns):
        """By hand, with Counter, over the birds only: the share of filled contact columns that differ from the
        modal residue of the birds, weighted by the number of columns of each gene."""
        sp, _ = birds
        slugs = list(sp.slug)
        cols = sp_.contact_columns(pd.read_csv(SITES, sep="\t"))
        got = sp_.contact_mismatch(alns, slugs, cols)
        for s in slugs[:: max(1, len(slugs) // 15)]:
            num = den = 0.0
            for g in sp_.COX:
                c = cols["contact"][g]
                diff = filled = 0
                for j in c:
                    column = [alns[g][x][j] for x in slugs if alns[g][x][j] != "-"]
                    if not column:          # a column every bird leaves empty has no consensus and counts nowhere
                        continue
                    mode = Counter(sorted(column)).most_common(1)[0][0]
                    r = alns[g][s][j]
                    if r != "-":
                        filled += 1; diff += r != mode
                num += diff / filled * len(c); den += len(c)
            assert got.loc[s, "contact"] == pytest.approx(num / den, abs=1e-12)

    def test_the_test_variable_is_contact_minus_control(self, birds, alns):
        sp, _ = birds
        got = sp_.contact_mismatch(alns, list(sp.slug), sp_.contact_columns(pd.read_csv(SITES, sep="\t")))
        assert np.allclose(got["contact minus control"], got.contact - got.control)
        assert got.contact.between(0, 1).all() and got.control.between(0, 1).all()

    def test_the_contact_map_has_298_contacts_and_as_many_controls(self):
        cols = sp_.contact_columns(pd.read_csv(SITES, sep="\t"))
        assert sum(len(v) for v in cols["contact"].values()) == 298
        assert {g: len(v) for g, v in cols["contact"].items()} == {g: len(v) for g, v in cols["control"].items()}


@needs_build
class TestHydropathy:
    def test_total_gravy_is_the_kyte_doolittle_mean_of_the_amino_acid_content(self):
        """The build computes GRAVY from codon counts; composition.tsv counts residues in the alignments. Two
        independent paths to the same number: they agree to within 0.05 everywhere and 0.002 in the median."""
        g = sp_.gravy_table(paths.BUILD)["gravy:total_pc"]
        c = pd.read_csv(COMPOSITION, sep="\t").set_index("slug")
        aa = c[[f"aa_{a}" for a in KD]]
        est = (aa * pd.Series({f"aa_{a}": v for a, v in KD.items()})).sum(axis=1) / aa.sum(axis=1)
        d = (g - est).dropna().abs()
        assert len(d) > 2000
        assert d.max() < 0.05 and d.median() < 0.002

    def test_one_column_per_region_and_no_retired_record(self):
        g = sp_.gravy_table(paths.BUILD)
        assert list(g.columns) == [f"gravy:{r}" for r in sp_.GRAVY_REGIONS]
        sp = pd.read_parquet(paths.BUILD / "species.parquet", columns=["slug", "record_status"])
        assert not set(g.index) & set(sp.slug[sp.record_status.eq("retired")])


@needs_build
class TestFits:
    def test_a_slope_is_per_standard_deviation_whatever_the_unit(self, birds):
        sp, C = birds
        g = sp_.gravy_table(paths.BUILD).reindex(sp.slug)[["gravy:total_pc"]]
        a = pd.DataFrame(sp_.fit_features(sp, g, C, "Aves", "hydropathy"))
        b = pd.DataFrame(sp_.fit_features(sp, g * 37.0 - 5.0, C, "Aves", "hydropathy"))
        assert np.allclose(a.slope_per_sd, b.slope_per_sd, rtol=1e-6) and np.allclose(a.p, b.p, rtol=1e-6)

    def test_the_interval_is_the_slope_plus_minus_t_times_se(self, birds):
        from scipy import stats
        sp, C = birds
        r = pd.DataFrame(sp_.fit_features(sp, sp_.gravy_table(paths.BUILD).reindex(sp.slug), C, "Aves", "hydropathy"))
        df = r.n - 1 - (1 + r.model.map({m: len(c) for m, c in sp_.MODELS.items()}))
        h = stats.t.ppf(0.975, df) * r.se
        assert np.allclose(r.ci95_low, r.slope_per_sd - h) and np.allclose(r.ci95_high, r.slope_per_sd + h)

    def test_every_model_of_a_feature_is_fitted_on_the_same_species(self, birds):
        """With a third of the values missing, both models use the species that have one, and only those: the
        degrees of freedom of each fit (read back from its interval) follow n, not the class size."""
        from scipy import stats
        sp, C = birds
        g = sp_.gravy_table(paths.BUILD).reindex(sp.slug)[["gravy:COX1"]].copy()
        g.iloc[::3] = np.nan
        r = pd.DataFrame(sp_.fit_features(sp, g, C, "Aves", "hydropathy"))
        assert r.n.eq(int(g.notna().sum().iloc[0])).all() and r.n.lt(len(sp)).all()
        k = r.model.map({m: len(c) for m, c in sp_.MODELS.items()})
        assert np.allclose((r.ci95_high - r.slope_per_sd) / r.se, stats.t.ppf(0.975, r.n - 2 - k))

    def test_too_few_values_is_not_assessed_never_dropped(self, birds):
        sp, C = birds
        g = sp_.gravy_table(paths.BUILD).reindex(sp.slug)[["gravy:ND6"]].copy()
        g.iloc[sp_.MIN_SPECIES - 1:] = np.nan
        r = pd.DataFrame(sp_.fit_features(sp, g, C, "Aves", "hydropathy"))
        assert len(r) == len(sp_.MODELS)
        assert r.status.str.startswith("not assessed").all() and "slope_per_sd" not in r


@pytest.fixture(scope="module")
def res(tmp_path_factory):
    out = tmp_path_factory.mktemp("h12")
    sp_.main(paths.BUILD, COMPOSITION, FOURFOLD, paths.DATED_TREE, paths.MITOPORT_IFACE, out,
             classes=("Aves", "Amphibia", "Insecta"))
    return pd.read_csv(out / "structure_pgls.tsv", sep="\t")


@needs_build
@needs_sites
class TestRun:
    def test_q_is_benjamini_hochberg_within_class_family_and_model(self, res):
        t = res[res.status.eq("tested")]
        assert t.groupby(["class", "family", "model"]).ngroups >= 4
        for _, d in t.groupby(["class", "family", "model"]):
            assert np.allclose(d.q, multipletests(d.p, method="fdr_bh")[1], rtol=1e-4)   # the table keeps 5 digits
        assert res.q[res.status.ne("tested")].isna().all()

    def test_a_class_too_small_is_written_as_not_assessed(self, res):
        ins = res[res["class"].eq("Insecta")]
        assert len(ins) == len(sp_.FAMILIES) * len(sp_.MODELS)
        assert ins.status.str.startswith("not assessed").all()

    def test_every_family_and_variable_is_in_the_table(self, res):
        t = res[res.status.eq("tested")]
        assert set(t[t.family.eq("contacts")].feature) == set(sp_.CONTACT_VARIABLES)
        assert set(t[t.family.eq("hydropathy")].feature) == {f"gravy:{r}" for r in sp_.GRAVY_REGIONS}

    def test_a_missing_contact_map_is_refused(self, tmp_path):
        with pytest.raises(SystemExit, match="sites.tsv is missing"):
            sp_.main(paths.BUILD, COMPOSITION, FOURFOLD, paths.DATED_TREE, tmp_path, tmp_path, classes=("Aves",))
