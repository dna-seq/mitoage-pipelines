"""H17: the co-evolution score and the bookkeeping of the ESM-2 scores, on cases small enough to work out by hand.
No model is loaded here: `cybrid_esm.score` takes the log-probabilities as a table."""
import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines.mitoport import cybrid, cybrid_coev as cc, cybrid_esm as ce


class TestWeights:
    def test_identical_rows_share_their_weight(self):
        w = cc.weights({"a": "MKVLAWYRTE", "b": "MKVLAWYRTE", "c": "GGGGGGGGGG"})
        assert w == {"a": 0.5, "b": 0.5, "c": 1.0}

    def test_the_threshold_and_gaps(self):
        # a and b agree at 8 of the 10 columns both fill: exactly 0.8, clustered; c fills 5 columns and agrees with
        # each at 3 of them (0.6): alone. A fourth row that agrees with a at 4 of the 5 it fills (0.8) joins a only.
        w = cc.weights({"a": "MKVLAWYRTE", "b": "MKVLAWYRGG", "c": "MKVGG-----"})
        assert w == {"a": 0.5, "b": 0.5, "c": 1.0}
        w = cc.weights({"a": "MKVLAWYRTE", "b": "MKVLAWYRGG", "d": "-----WYRTA"})
        assert w["a"] == pytest.approx(1 / 3) and w["b"] == 0.5 and w["d"] == 0.5
        assert cc.weights({}) == {}


class TestJointAndPmi:
    def cols(self):
        # 40 species: A with K in 20, D with R in 20: the two columns covary perfectly
        a = {f"s{i}": ("A" if i < 20 else "D") for i in range(40)}
        b = {f"s{i}": ("K" if i < 20 else "R") for i in range(40)}
        return a, b, {s: 1.0 for s in a}

    def test_joint_by_hand(self):
        a, b, w = self.cols()
        p, n, neff = cc.joint(a, b, w, theta=0.05)
        assert (n, neff) == (40, 40.0) and p.sum() == pytest.approx(1.0)
        assert p[cc.IDX["A"], cc.IDX["K"]] == pytest.approx(0.95 * 0.5 + 0.05 / 400)
        assert p[cc.IDX["A"], cc.IDX["R"]] == pytest.approx(0.05 / 400)

    def test_pmi_and_loss_by_hand(self):
        a, b, w = self.cols()
        p, _, _ = cc.joint(a, b, w, theta=0.05)
        pa = 0.95 * 0.5 + 0.05 / 20
        assert cc.pmi(p, "A", "K") == pytest.approx(np.log((0.95 * 0.5 + 0.05 / 400) / pa ** 2))
        # host A with K, donor D: the chimera pairs D with K, never seen together
        want = np.log((0.475 + 0.000125) / pa ** 2) - np.log(0.000125 / pa ** 2)
        assert cc.loss(p, "A", "D", "K") == pytest.approx(want) and want > 8
        assert cc.loss(p, "A", "A", "K") == 0.0

    def test_skipped_species_gaps_and_weights(self):
        a, b, w = self.cols()
        a["s0"] = "-"; b["s1"] = "X"; w["s2"] = 0.5
        p, n, neff = cc.joint(a, b, w, skip=frozenset({"s3"}))
        assert n == 37 and neff == pytest.approx(36.5)

    def test_independent_columns_lose_nothing(self):
        a = {f"s{i}": "AD"[i % 2] for i in range(40)}
        b = {f"s{i}": "KR"[(i // 2) % 2] for i in range(40)}
        p, _, _ = cc.joint(a, b, {s: 1.0 for s in a})
        assert cc.loss(p, "A", "D", "K") == pytest.approx(0.0, abs=1e-9)


class TestScorePair:
    def setup(self):
        # one mitochondrial gene of 3 columns (0 contact, 1 its control, 2 unused), one nuclear gene of 1 column
        mt = {f"s{i}": ("A" if i < 20 else "D") + "AD"[i % 2] + "G" for i in range(40)}
        nuc = {f"s{i}": ("K" if i < 20 else "R") for i in range(40)}
        mt["host"] = "AAG"; nuc["host"] = "K"; mt["donor"] = "DDG"; nuc["donor"] = "R"
        pr = pd.DataFrame([{"mt_gene": "X", "mt_col": 0, "ctl_col": 1, "nuc_gene": "N", "nuc_col": 0}])
        return pr, {"X": mt}, {"N": nuc}, {s: 1.0 for s in mt}

    def test_contact_loses_and_control_does_not(self):
        pr, mt, nuc, w = self.setup()
        host, donor = {"X": {0: "A", 1: "A", 2: "G"}}, {"X": {0: "D", 1: "D", 2: "G"}}
        s = cc.score_pair(pr, mt, nuc, w, "host", "donor", host, donor, frozenset({"host", "donor"}))
        assert (s["pairs_scored"], s["pairs_changed"], s["pairs_unscored"]) == (1, 1, 0)
        assert s["coev_contact"] > 8 and abs(s["coev_control"]) < 1e-9
        assert s["coev_excess"] == pytest.approx(s["coev_contact"] - s["coev_control"]) == s["coev_excess_mean"]

    def test_host_and_donor_are_left_out_of_the_counts(self):
        pr, mt, nuc, w = self.setup()
        host, donor = {"X": {0: "A", 1: "A"}}, {"X": {0: "D", 1: "D"}}
        nuc["N"]["donor"] = "K"                          # were the donor counted, D would have been seen with K
        left_out = cc.score_pair(pr, mt, nuc, w, "host", "donor", host, donor, frozenset({"host", "donor"}))
        counted = cc.score_pair(pr, mt, nuc, w, "host", "donor", host, donor, frozenset())
        assert left_out["coev_contact"] > counted["coev_contact"] + 1

    def test_unknown_residues_and_thin_pairs_are_unscored_not_zero(self):
        pr, mt, nuc, w = self.setup()
        host = {"X": {0: "A", 1: "A"}}
        s = cc.score_pair(pr, mt, nuc, w, "host", "donor", host, {"X": {0: "?", 1: "D"}}, frozenset())
        assert s["pairs_unscored"] == 1 and np.isnan(s["coev_excess"]) and np.isnan(s["coev_excess_mean"])
        few = {k: v for k, v in list(mt["X"].items())[:10]}
        s = cc.score_pair(pr, {"X": {**few, "host": "AAG"}}, nuc, w, "host", "donor", host, {"X": {0: "D", 1: "D"}}, frozenset())
        assert s["pairs_unscored"] == 1 and np.isnan(s["coev_excess"])


class TestColumns:
    def test_a_protein_on_the_columns_of_the_human_row(self):
        human_row = "M-KVLAW--YRTE"
        got = cc.on_columns(human_row, "MKVLAYRSE", partial=False)        # W deleted, T -> S
        assert got == {0: "M", 2: "K", 3: "V", 4: "L", 5: "A", 6: "-", 9: "Y", 10: "R", 11: "S", 12: "E"}
        idx = ce.index_on_columns(human_row, "MKVLAYRSE")
        assert idx == {0: 0, 2: 1, 3: 2, 4: 3, 5: 4, 9: 5, 10: 6, 11: 7, 12: 8}


class TestEsmBookkeeping:
    def setup(self):
        host_seq = {"X": "MKWLA"}; host_idx = {"X": {c: c for c in range(5)}}
        donor = {"X": {0: "M", 1: "R", 2: "W", 3: "I", 4: "?"}}            # K->R at column 1, L->I at column 3
        pr = pd.DataFrame([{"mt_gene": "X", "mt_col": 1, "ctl_col": 3, "nuc_gene": "N", "nuc_col": 0}])
        return host_seq, host_idx, donor, pr

    def lp(self, **aa):
        v = np.full(20, -5.0)
        for a, x in aa.items():
            v[cc.IDX[a]] = x
        return v

    def test_tasks_mask_only_differing_positions_and_skip_long_pairs(self):
        host_seq, host_idx, donor, pr = self.setup()
        t, too_long = ce.tasks(host_seq, host_idx, {"d": donor}, pr, {"N": "AAAA"})
        assert too_long == 0 and t[("X", ce.ALONE)] == ("MKWLA", 0, [1, 3])
        assert t[("X", "N")] == ("AAAA" + ce.LINKER + "MKWLA", 4 + 25, [1, 3])
        t, too_long = ce.tasks(host_seq, host_idx, {"d": donor}, pr, {"N": "A" * 1000})
        assert too_long == 1 and ("X", "N") not in t

    def test_scores_by_hand(self):
        host_seq, host_idx, donor, pr = self.setup()
        logp = {("X", ce.ALONE, 1): self.lp(K=-1.0, R=-3.0), ("X", ce.ALONE, 3): self.lp(L=-1.0, I=-2.0),
                ("X", "N", 1): self.lp(K=-0.5, R=-4.5), ("X", "N", 3): self.lp(L=-1.0, I=-2.5)}
        s = ce.score(pr, {"X": {1}}, {"X": {3}}, host_seq, host_idx, donor, logp)
        assert (s["esm_contact"], s["esm_control"], s["esm_all"]) == (2.0, 1.0, 3.0)
        assert s["positions_unknown"] == 1 and s["n_contact"] == 1 and s["esm_contact_mean"] == 2.0
        # with the partner: contact LLR -4 against -2 alone (effect 2), control -1.5 against -1 (effect 0.5)
        assert s["esm_partner_effect"] == 2.0 and s["esm_partner_effect_control"] == 0.5 and s["esm_partner_excess"] == 1.5

    def test_a_partner_that_was_not_run_is_not_assessed(self):
        host_seq, host_idx, donor, pr = self.setup()
        logp = {("X", ce.ALONE, 1): self.lp(K=-1.0, R=-3.0), ("X", ce.ALONE, 3): self.lp(L=-1.0, I=-2.0)}
        s = ce.score(pr, {"X": {1}}, {"X": {3}}, host_seq, host_idx, donor, logp)
        assert s["partner_residues_assessed"] == 0 and np.isnan(s["esm_partner_excess"]) and s["esm_contact"] == 2.0


class TestReadout:
    def test_missing_tables_are_said_and_left_out(self, tmp_path, capsys):
        key = ["host", "donor", "complex"]
        meas = cybrid.read_measurements()
        pts = meas[key].drop_duplicates()
        pts = pts[pts.complex != "all"]
        allp = pd.DataFrame([(h, d, c) for h, d in set(zip(meas.host, meas.donor)) for c in cybrid.STRUCTURES], columns=key)
        rng = np.random.default_rng(0)
        allp.assign(d_contact=rng.integers(0, 50, len(allp)), d_all=rng.integers(0, 90, len(allp)), f_all=rng.random(len(allp))
                    ).to_csv(tmp_path / "scores.tsv", sep="\t", index=False)
        with pytest.raises(SystemExit):
            cybrid.readout(tmp_path)
        allp.assign(coev_excess=rng.random(len(allp)), coev_excess_mean=rng.random(len(allp))).to_csv(tmp_path / "coev_scores.tsv", sep="\t", index=False)
        s = cybrid.readout(tmp_path)
        assert "NOT ASSESSED" in capsys.readouterr().out and "esm_scores.tsv" not in set(s.score)
        assert {"coev_excess", "d_all"} <= set(s.score[s.test == "T1"]) and not s.score.str.startswith("esm").any()
        assert set(s.test) == {"T1", "T2", "T2b", "T3", "T4"}
