"""H37 part A: the rate of the RNA genes by H7's relative-rate test, and the positive control (`longevity/rna_rates.py`).

Real data where the code reads it: real tRNA sequences of real mitogenomes (`data/input/cybrid/records.gb`) aligned
by MAFFT, H7's own published table (`results/review-20261118/longevity/rates.tsv`) for the comparison, and the downloaded
build for the equivalence with H7's code. Expected values are computed in the tests from the same inputs.
"""
import io
import json
import math
import re
import shutil
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio import Phylo
from statsmodels.stats.multitest import multipletests

from mitoage_pipelines import paths
from mitoage_pipelines.longevity import rates, rna_genes as rg, rna_rates as rr

RECORDS = paths.CYBRID_GB / "records.gb"
PUBLISHED = Path(__file__).resolve().parents[1] / "results" / "review-20261118" / "longevity" / "rates.tsv"
needs_records = pytest.mark.skipif(not RECORDS.is_file(), reason="needs `mitoage download cybrid`")
needs_mafft = pytest.mark.skipif(shutil.which("mafft") is None, reason="needs MAFFT on PATH")
COMP = paths.LONGEVITY_OUT / "composition.tsv"
FOUR = paths.LONGEVITY_OUT / "fourfold.tsv"
needs_build = pytest.mark.skipif(
    not all(p.is_file() for p in (paths.BUILD / "species.parquet", COMP, FOUR, paths.DATED_TREE)),
    reason="needs the build, `longevity controls`, `fourfold` and `dated-tree`")

SPECIES = {"NC_012920": "human", "NC_001643": "chimp", "NC_001644": "bonobo", "NC_011120": "gorilla",
           "NC_001665": "rat", "NC_005089": "mouse", "NC_007936": "hamster"}
LIFESPAN = dict(human=122.5, chimp=59.4, bonobo=55.0, gorilla=55.4, rat=3.8, mouse=4.0, hamster=3.9)
MASS = dict(human=62e3, chimp=45e3, bonobo=39e3, gorilla=140e3, rat=300, mouse=20, hamster=100)
NEWICK = "((((human:5,chimp:7):1,bonobo:7):5,gorilla:12):80,((mouse:12,rat:14):20,hamster:33):60);"


def jc(p):
    return -0.75 * math.log(1 - 4 * p / 3)


@pytest.fixture(scope="module")
def tree_file(tmp_path_factory):
    f = tmp_path_factory.mktemp("tree") / "t.nwk"
    f.write_text(NEWICK)
    return f


@pytest.fixture(scope="module")
def alignments():
    """tRNA-Val, tRNA-Leu(UUR) and tRNA-Lys of the seven real records, aligned by MAFFT as the build aligns."""
    if shutil.which("mafft") is None or not RECORDS.is_file():
        pytest.skip("needs MAFFT and the cybrid records")
    texts = [t.lstrip() + "\n//\n" for t in RECORDS.read_text().split("\n//\n") if t.strip()]
    recs = {re.search(r"^VERSION\s+(\S+)\.\d+", t, re.M).group(1): t for t in texts}
    genes = ["tRNA-Val", "tRNA-Leu(UUR)", "tRNA-Lys"]
    seqs = {g: {} for g in genes}
    for acc, name in SPECIES.items():
        got = rg.extract_record(recs[acc], "Chordata")["sequences"]
        for g in genes:
            seqs[g][name] = got[g]
    aligned, _ = rg.align_genes(seqs, 3)
    return {g: {k: v.upper() for k, v in ((r.split("\n")[0], r.split("\n")[1]) for r in t.split(">")[1:])} for g, t in aligned.items()}


def species_frame(names):
    return pd.DataFrame({"slug": names, "class": "Mammalia", "ll": [math.log10(LIFESPAN[n]) for n in names],
                         "lm": [math.log10(MASS[n]) for n in names]})


class TestDistance:
    @needs_records
    def test_identical_rows_are_zero_apart_and_the_jukes_cantor_correction_is_applied(self, alignments):
        s = alignments["tRNA-Val"]["human"].replace("-", "")[:60]
        assert rr.nt_distance(s, s) == {"dist": 0.0, "sites": 60}
        t = list(s)
        for i in range(0, 60, 10):                                        # six substitutions in 60 columns: p = 0.1
            t[i] = "A" if t[i] != "A" else "C"
        got = rr.nt_distance(s, "".join(t))
        assert got["sites"] == 60 and got["dist"] == pytest.approx(jc(0.1))
        assert got["dist"] > 0.1                                          # the correction, not the p-distance

    @needs_records
    def test_gaps_and_ambiguity_codes_are_not_columns_and_the_floor_is_thirty_columns(self, alignments):
        s = alignments["tRNA-Val"]["human"].replace("-", "")[:40]
        t = "N" * 5 + "-" * 5 + s[10:]
        got = rr.nt_distance(s, t)
        assert got["sites"] == 30 and got["dist"] == 0.0
        short = rr.nt_distance(s, "N" * 6 + s[6:])
        assert short["sites"] == 34 and not np.isnan(short["dist"])
        below = rr.nt_distance(s, "-" * 11 + s[11:])
        assert below["sites"] == 29 and np.isnan(below["dist"])

    def test_a_distance_at_saturation_is_unknown_not_a_number(self):
        s = "ACGT" * 15
        t = "CATG" * 15                                                    # every column differs: p = 1 > 0.75
        got = rr.nt_distance(s, t)
        assert got["sites"] == 60 and np.isnan(got["dist"])

    @needs_records
    def test_both_sisters_are_measured_on_the_columns_all_three_rows_resolve(self, alignments):
        a, b, o = (alignments["tRNA-Leu(UUR)"][k] for k in ("human", "chimp", "bonobo"))
        o = o[:3] + "N" + o[4:]                                           # one column the outgroup does not resolve
        b = b[:8] + "N" + b[9:]                                           # and one the second sister does not
        keep = [i for i in range(len(a)) if all(r[i] in "ACGT" for r in (a, b, o))]
        assert 3 not in keep and 8 not in keep and a[8] in "ACGT"
        da, db = rr.nt_outgroup_distances(a, b, o)
        for got, x in ((da, a), (db, b)):
            p = sum(x[i] != o[i] for i in keep) / len(keep)
            assert got["sites"] == len(keep) and got["dist"] == pytest.approx(jc(p) if p else 0.0)

    def test_rows_of_different_lengths_are_refused(self):
        with pytest.raises(ValueError, match="equal lengths"):
            rr.nt_outgroup_distances("ACGT", "ACG", "ACGT")


class TestContrasts:
    @needs_records
    @needs_mafft
    def test_a_cherry_is_the_difference_of_its_sisters_distances_to_the_outgroup(self, alignments, tree_file):
        sp = species_frame(list(SPECIES.values()))
        got = rr.nt_contrasts(sp, tree_file, {"tRNA-Val": alignments["tRNA-Val"]})
        assert sorted(zip(got.sp_a, got.sp_b, got.outgroup)) == [("human", "chimp", "bonobo"), ("mouse", "rat", "hamster")]
        for r in got.itertuples():
            ali = alignments["tRNA-Val"]
            keep = [i for i in range(len(ali[r.sp_a])) if all(ali[k][i] in "ACGT" for k in (r.sp_a, r.sp_b, r.outgroup))]
            d = lambda x: jc(sum(ali[x][i] != ali[r.outgroup][i] for i in keep) / len(keep))
            assert r.sites == len(keep)
            assert r.d_dist == pytest.approx(d(r.sp_a) - d(r.sp_b))
            assert r.d_ll == pytest.approx(math.log10(LIFESPAN[r.sp_a]) - math.log10(LIFESPAN[r.sp_b]))
            assert r.d_lm == pytest.approx(math.log10(MASS[r.sp_a]) - math.log10(MASS[r.sp_b]))
            assert r.time_my == (6.0 if r.sp_a == "human" else 13.0)
            assert r.gene == "tRNA-Val"

    @needs_records
    @needs_mafft
    def test_a_gene_a_sister_lacks_makes_no_cherry_and_each_gene_has_its_own_species(self, alignments, tree_file):
        sp = species_frame(list(SPECIES.values()))
        lys = {k: v for k, v in alignments["tRNA-Lys"].items() if k != "chimp"}
        got = rr.nt_contrasts(sp, tree_file, {"tRNA-Lys": lys, "tRNA-Val": alignments["tRNA-Val"]})
        assert Counter(got.gene) == {"tRNA-Val": 2, "tRNA-Lys": 1}
        lys = got[got.gene == "tRNA-Lys"]
        assert list(zip(lys.sp_a, lys.sp_b)) == [("mouse", "rat")]
        nothing = rr.nt_contrasts(sp[sp.slug.isin(["human", "mouse"])], tree_file, {"tRNA-Val": alignments["tRNA-Val"]})
        assert nothing.empty

    @needs_records
    @needs_mafft
    def test_a_species_with_no_lifespan_is_not_in_the_frame_and_so_in_no_cherry(self, alignments, tree_file):
        sp = species_frame([n for n in SPECIES.values() if n != "rat"])
        got = rr.nt_contrasts(sp, tree_file, {"tRNA-Val": alignments["tRNA-Val"]})
        assert list(zip(got.sp_a, got.sp_b)) == [("human", "chimp")]

    @needs_records
    @needs_mafft
    def test_the_pair_is_flagged_saturated_when_a_distance_passes_one(self, alignments, tree_file):
        sp = species_frame(list(SPECIES.values()))
        al = dict(alignments["tRNA-Val"])
        ref = al["bonobo"]
        al["chimp"] = ref                                              # the sister that equals the outgroup
        other = {"A": "C", "C": "G", "G": "T", "T": "A"}
        al["human"] = "".join(other[ch] if ch in other and i % 5 < 3 else ch for i, ch in enumerate(ref))
        got = rr.nt_contrasts(sp, tree_file, {"g": al})
        row = got[got.sp_a == "human"].iloc[0]
        keep = [i for i in range(len(ref)) if ref[i] in "ACGT"]
        p = sum(al["human"][i] != ref[i] for i in keep) / len(keep)
        assert 0.5 < p < 0.75 and row.dist_a == pytest.approx(jc(p)) and row.dist_a > 1 and row.dist_b == 0.0
        assert row.saturated and not got[got.sp_a == "mouse"].iloc[0].saturated


class TestTable:
    def frame(self):
        rng = np.random.default_rng(7)
        rows = []
        for gene, slope in (("g1", 0.9), ("g2", 0.0), ("g3", -0.5), ("g4", 0.1)):
            for i in range(30):
                dll = rng.normal(0, 0.3)
                rows.append({"class": "Mammalia", "gene": gene, "d_ll": dll, "d_lm": rng.normal(0, 1), "time_my": 20.0 + i,
                             "d_dist": slope * dll + rng.normal(0, 0.1), "saturated": i >= 27})
        return pd.DataFrame(rows)

    def test_holm_q_and_the_partial_correlation_follow_the_tests_within_each_group(self):
        c = self.frame()
        t = rr.test_table(c, ["g1", "g2", "g3", "g4"], rr.NT_NORMALISATIONS, sort=False)
        for (cls, norm, sat), g in t.groupby(["class", "normalisation", "includes_saturated"]):
            assert len(g) == 4
            assert np.allclose(g.p_holm, multipletests(g.p, method="holm")[1], rtol=1e-9, atol=0)
            assert np.allclose(g.q, multipletests(g.p, method="fdr_bh")[1], rtol=1e-9, atol=0)
            cov = len(rr.NT_NORMALISATIONS[norm].get("covariates", []))
            assert np.allclose(g.r, g.t / np.sqrt(g.t ** 2 + g.n - 1 - cov), rtol=1e-9, atol=0)
        mid = t[(t["class"] == "Mammalia") & (t.normalisation == "d_raw") & ~t.includes_saturated].sort_values("p")
        assert not np.allclose(mid.p_holm, mid.q, rtol=1e-9, atol=0)          # the two corrections differ on this frame
        by = t[(t["class"] == "Mammalia") & (t.normalisation == "d_raw") & ~t.includes_saturated].set_index("gene")
        assert by.loc["g1", "slope"] > 0.5 > 0.1 > abs(by.loc["g2", "slope"]) and by.loc["g3", "slope"] < -0.3
        assert by.loc["g1", "n"] == 27 and t[(t["class"] == "Mammalia") & t.includes_saturated & (t.normalisation == "d_raw")].n.eq(30).all()
        direct = rates.test_contrasts(c[(c.gene == "g1") & ~c.saturated], y="d_dist")
        assert by.loc["g1", "p"] == pytest.approx(direct["p"]) and by.loc["g1", "sign_p"] == pytest.approx(direct["sign_p"])
        assert by.loc["g1", "p_holm"] <= 0.05 < by.loc["g2", "p_holm"]

    def test_a_gene_with_too_few_cherries_has_no_row(self):
        c = self.frame()
        c = c[~((c.gene == "g2") & (c.index % 30 >= rates.MIN_CONTRASTS - 1))]
        t = rr.test_table(c, ["g1", "g2", "g3", "g4"], rr.NT_NORMALISATIONS, sort=False)
        assert "g2" not in set(t.gene) and {"g1", "g3", "g4"} == set(t.gene)

    def test_the_effect_size_is_the_partial_correlation_and_unknown_when_the_test_is(self):
        assert rr.partial_r(2.0, 8) == pytest.approx(2 / math.sqrt(12))
        assert rr.partial_r(-2.0, 8) == pytest.approx(-2 / math.sqrt(12))
        assert np.isnan(rr.partial_r(float("nan"), 8)) and np.isnan(rr.partial_r(1.0, 0))


class TestControl:
    @pytest.mark.skipif(not PUBLISHED.is_file(), reason="needs results/review-20261118")
    def test_a_table_equal_to_the_published_one_is_reproduced_and_one_changed_cell_is_not(self):
        ref = pd.read_csv(PUBLISHED, sep="\t")
        assert len(ref) == 655
        same = rr.compare(ref.copy(), ref)
        assert same["reproduced"] and same["rows_that_differ"] == 0
        near = ref.copy()
        near.loc[3, "slope"] = near.loc[3, "slope"] * (1 + 1e-9)             # below the precision the table is written at
        assert rr.compare(near, ref)["reproduced"]
        moved = ref.copy()
        moved.loc[3, "slope"] = ref.loc[3, "slope"] * 1.001
        got = rr.compare(moved, ref)
        assert not got["reproduced"] and got["rows_that_differ"] == 1 and got["differing_columns"] == {"slope": 1}

    @pytest.mark.skipif(not PUBLISHED.is_file(), reason="needs results/review-20261118")
    def test_a_missing_row_an_extra_row_and_a_missing_column_are_each_a_failure(self):
        ref = pd.read_csv(PUBLISHED, sep="\t")
        got = rr.compare(ref.iloc[1:].copy(), ref)
        assert not got["reproduced"] and got["only_in_reference"] == 1 and got["only_in_mine"] == 0
        extra = pd.concat([ref, ref.iloc[[0]].assign(gene="NEW")], ignore_index=True)
        got = rr.compare(extra, ref)
        assert not got["reproduced"] and got["only_in_mine"] == 1
        got = rr.compare(ref.drop(columns=["q"]), ref)
        assert not got["reproduced"] and got["columns_missing_in_mine"] == ["q"]

    def test_the_gate_refuses_a_missing_control_another_build_and_a_failed_control(self, tmp_path):
        data, out = tmp_path / "data", tmp_path / "out"
        data.mkdir(), out.mkdir()
        (data / "species.parquet").write_bytes(b"one build")
        with pytest.raises(SystemExit, match="rna-control"):
            rr.require_control(data, out)
        ok = {"reproduced": True, "species_sha256": rr._sha(data / "species.parquet")}
        (out / rr.CONTROL_FILE).write_text(json.dumps({**ok, "species_sha256": "0" * 64}))
        with pytest.raises(SystemExit, match="another build"):
            rr.require_control(data, out)
        (out / rr.CONTROL_FILE).write_text(json.dumps({**ok, "reproduced": False, "rows_that_differ": 4, "only_in_mine": 0, "only_in_reference": 1}))
        with pytest.raises(SystemExit, match="cannot reproduce H7"):
            rr.require_control(data, out)
        (out / rr.CONTROL_FILE).write_text(json.dumps({**ok, "reproduced": "yes"}))      # only the boolean true passes
        with pytest.raises(SystemExit):
            rr.require_control(data, out)
        (out / rr.CONTROL_FILE).write_text(json.dumps(ok))
        assert rr.require_control(data, out)["reproduced"] is True

    def test_no_rna_rate_is_made_without_the_control(self, tmp_path):
        with pytest.raises(SystemExit, match="rna-control"):
            rr.rna_rates(tmp_path, tmp_path / "c", tmp_path / "f", tmp_path / "t", tmp_path, tmp_path)


class TestOnTheBuild:
    @needs_build
    def test_the_control_reproduces_h7s_own_table_and_stops_the_pipeline_when_it_does_not(self, tmp_path, monkeypatch):
        sp = rates.select(paths.BUILD, COMP, FOUR, paths.DATED_TREE)
        tips = [t.name for t in Phylo.read(paths.DATED_TREE, "newick").get_terminals()]
        order = {s: i for i, s in enumerate(tips)}
        part = sp[sp["class"].isin(["Mammalia", "Aves"])].assign(i=lambda d: d.slug.map(order)).sort_values("i").head(500).drop(columns="i").reset_index(drop=True)
        monkeypatch.setattr(rates, "select", lambda *a, **k: part)        # one subset for H7's own code and for the control
        h7 = tmp_path / "h7"
        rates.main(paths.BUILD, COMP, FOUR, paths.DATED_TREE, h7)
        assert len(pd.read_csv(h7 / "rates.tsv", sep="\t")) > 100
        out = tmp_path / "out"
        v = rr.control(paths.BUILD, COMP, FOUR, paths.DATED_TREE, h7 / "rates.tsv", out)
        assert v["reproduced"] is True and v["rows_that_differ"] == 0 and v["only_in_mine"] == 0 and v["rows_mine"] == v["rows_reference"]
        assert (out / "h7_control_rates.tsv").read_bytes() == (h7 / "rates.tsv").read_bytes()     # the same bytes, not only the same numbers
        assert rr.require_control(paths.BUILD, out)["reproduced"] is True
        ref = pd.read_csv(h7 / "rates.tsv", sep="\t")
        ref.loc[5, "slope"] = ref.loc[5, "slope"] * 1.01
        ref.to_csv(h7 / "tampered.tsv", sep="\t", index=False)
        bad = rr.control(paths.BUILD, COMP, FOUR, paths.DATED_TREE, h7 / "tampered.tsv", out)
        assert bad["reproduced"] is False and bad["rows_that_differ"] == 1
        with pytest.raises(SystemExit, match="cannot reproduce H7"):
            rr.require_control(paths.BUILD, out)
        with pytest.raises(SystemExit, match="cannot reproduce H7"):
            rr.rna_rates(paths.BUILD, COMP, FOUR, paths.DATED_TREE, out, out)
        # and the command says so with its exit code, so that a job that depends on it does not start
        from typer.testing import CliRunner
        from mitoage_pipelines.cli import app
        run = lambda ref: CliRunner().invoke(app, ["longevity", "rna-control", "--reference", str(ref), "--out", str(out)])
        assert run(h7 / "tampered.tsv").exit_code == 1
        assert run(h7 / "rates.tsv").exit_code == 0


    @needs_build
    def test_the_codon_mode_makes_h7s_contrasts_on_a_clade_of_mammals(self):
        sp = rates.select(paths.BUILD, COMP, FOUR, paths.DATED_TREE)
        tips = [t.name for t in Phylo.read(paths.DATED_TREE, "newick").get_terminals()]
        order = {s: i for i, s in enumerate(tips)}
        mam = sp[sp["class"] == "Mammalia"].assign(i=lambda d: d.slug.map(order)).sort_values("i").head(160).drop(columns="i").reset_index(drop=True)
        mine = rr.codon_contrasts(mam, paths.DATED_TREE, rr.codon_alignments(paths.BUILD))
        theirs = rates.build_contrasts(paths.BUILD, mam, paths.DATED_TREE)
        assert len(theirs) > 14 * 20
        pd.testing.assert_frame_equal(mine.reset_index(drop=True), theirs.reset_index(drop=True))

    @needs_build
    def test_the_cherries_are_those_of_the_species_every_gene_has(self):
        sp = rates.select(paths.BUILD, COMP, FOUR, paths.DATED_TREE)
        tips = [t.name for t in Phylo.read(paths.DATED_TREE, "newick").get_terminals()]
        order = {s: i for i, s in enumerate(tips)}
        mam = sp[sp["class"] == "Mammalia"].assign(i=lambda d: d.slug.map(order)).sort_values("i").head(120).drop(columns="i").reset_index(drop=True)
        al = rr.codon_alignments(paths.BUILD)
        full = rr.codon_contrasts(mam, paths.DATED_TREE, al)
        lost = full.sp_a.iloc[0]
        al["ATP8"] = {k: v for k, v in al["ATP8"].items() if k != lost}      # one species the ATP8 alignment lacks
        got = rr.codon_contrasts(mam, paths.DATED_TREE, al)
        assert not ({lost} & set(got.sp_a) | {lost} & set(got.sp_b) | {lost} & set(got.outgroup))
        assert 0 < len(got) < len(full)

    @needs_build
    def test_the_high_set_is_the_default_dataset_at_confidence_high_and_the_environment_restriction_is_refused(self, monkeypatch):
        every = rr.species_for(paths.BUILD, COMP, FOUR, paths.DATED_TREE, None)
        high = rr.species_for(paths.BUILD, COMP, FOUR, paths.DATED_TREE, "high")
        level = pd.read_parquet(paths.BUILD / "species.parquet", columns=["slug", "data_confidence"]).set_index("slug").data_confidence
        assert set(high.slug) == {s for s in every.slug if level[s] == "high"}
        assert 0 < len(high) < len(every)
        monkeypatch.setenv("MITOAGE_SPECIES_SET", "high")
        with pytest.raises(SystemExit, match="unset it"):
            rr.species_for(paths.BUILD, COMP, FOUR, paths.DATED_TREE, None)
