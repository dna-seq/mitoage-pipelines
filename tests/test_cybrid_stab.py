"""H28: the bookkeeping of StaB-ddG's scores (the chains kept, the breaks, where a substitution lands in each of the
three descriptions, the sets, the collection, the rule's verdict) by hand; and, where the model and its weights are
installed (`uv sync --group stab`, `mitoage download stabddg`), that an empty set gives exactly 0 and a seed repeats."""
import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.mitoport import cybrid, cybrid_stab as cs

P = "MKVLAWYRTEQNDSGHIFPCMKVLAWYRTEQNDSGHIFPC"          # 40 residues
Q = "GSTNDEQKRHGSTNDEQKRHGSTNDEQKRHGSTNDEQKRH"


def chain(cid, seq, origin, num=None):
    """A straight chain along x from `origin`, 3.8 A a residue; the four backbone atoms 1 A apart."""
    n = len(seq)
    ca = np.array(origin, dtype=float) + np.c_[3.8 * np.arange(n), np.zeros(n), np.zeros(n)]
    bb = np.stack([ca + [-1, 0, 0], ca, ca + [1, 0, 0], ca + [1, 1, 0]], axis=1)
    return {"id": cid, "seq": seq, "num": list(num or range(1, n + 1)), "bb": bb, "heavy": bb.reshape(-1, 3)}


class TestSystem:
    chains = [chain("A", P, (0, 0, 0)), chain("B", Q, (0, 4, 0)), chain("C", Q, (0, 60, 0)), chain("D", P[::-1], (0, -4, 0)),
              chain("E", P, (0, 8, 0))]

    def test_the_mt_chains_every_copy_and_what_touches_them(self):
        kept, mt = cs.system(self.chains, {"X": P}, {"Y": P[::-1]})
        assert [c["id"] for c in kept] == ["A", "B", "E"]          # C is far; D is another complex's mtDNA-encoded protein
        assert sorted(mt) == [0, 2] and mt[0][0] == "X" and mt[0][1] == {i: i for i in range(40)}

    def test_a_gene_without_a_chain_is_refused(self):
        with pytest.raises(SystemExit):
            cs.system(self.chains, {"X": P, "Z": Q[::-1] + "WWWWWWWWWWWWWWWWWWWWWWWW"}, {})


class TestStructure:
    def test_breaks_become_masked_places_up_to_the_cap(self):
        c = chain("A", "MKVLAW", (0, 0, 0), num=[1, 2, 5, 6, 6, 200])
        seq, bb, where = cs.gapped(c)
        assert seq == "MK--VLA" + "-" * cs.GAP_CAP + "W"
        assert where == [0, 1, 4, 5, 6, 7 + cs.GAP_CAP]
        assert np.isnan(bb[2]).all() and np.isnan(bb[3]).all() and np.isfinite(bb[4]).all()
        assert (bb[4] == c["bb"][2]).all() and len(bb) == len(seq)

    def test_a_substitution_lands_on_its_residue_in_each_description(self):
        chains = [chain("A", "MKVLAW", (0, 0, 0), num=[1, 2, 5, 6, 7, 8]), chain("B", "GSTND", (0, 4, 0)), chain("C", "MKVLAW", (0, 8, 0))]
        whole, at = cs.structure("x", chains, [0, 1, 2])
        mt, at_mt = cs.structure("x_mt", chains, [0, 2])
        nuc, at_nuc = cs.structure("x_nuc", chains, [1])
        assert whole["seq"] == "MK--VLAW" + "GSTND" + "MKVLAW" and mt["seq"] == "MK--VLAW" + "MKVLAW" and nuc["seq"] == "GSTND"
        assert whole["num_of_chains"] == 3 and whole["masked_list"] == ["c0", "c1", "c2"] and whole["visible_list"] == []
        for d, a in ((whole, at), (mt, at_mt)):
            assert d["seq"][a[(0, 2)]] == "V" and d["seq"][a[(2, 5)]] == "W"
        assert (0, 2) not in at_nuc and nuc["seq"][at_nuc[(1, 4)]] == "D"
        assert np.array(whole["coords_chain_c1"]["CA_chain_c1"]).shape == (5, 3)
        subs = {(0, 2): "I", (2, 5): "F"}
        enc = cs.encode(whole["seq"], {at[k]: aa for k, aa in subs.items()})
        assert "".join(cs.ALPHABET[i] for i in enc) == "MKXXILAW" + "GSTND" + "MKVLAF"
        assert (cs.encode(nuc["seq"], {at_nuc[k]: aa for k, aa in subs.items() if k in at_nuc}) == cs.encode(nuc["seq"], {})).all()

    def test_a_substitution_at_a_break_is_refused(self):
        with pytest.raises(ValueError):
            cs.encode("MK--V", {2: "A"})


class TestSets:
    def test_sites_are_carried_to_the_host_by_alignment(self):
        human = {"X": "MKVLAWYRTE"}
        host = {"X": "MKLAWYRTEG"}                      # V3 deleted in the host: later positions sit one lower
        sites = pd.DataFrame({"gene": ["X"] * 4, "human_pos": [2, 3, 6, 9], "human_aa": list("KVWT"),
                              "group": ["contact", "contact", "control", "other"]})
        assert cs.site_groups(sites, human, host) == {"X": {1: "contact", 4: "control"}}
        sites.loc[0, "human_aa"] = "R"
        with pytest.raises(SystemExit):
            cs.site_groups(sites, human, host)

    def test_the_sets_and_every_copy(self):
        host = {"X": "MKVLAWYRTE"}
        chains = [chain("A", "KVLAWYRT", (0, 0, 0)), chain("B", "GSTND", (0, 4, 0)), chain("C", "KVLAWYRT", (0, 8, 0))]
        to_protein = {i: i + 1 for i in range(8)}
        mt_chain = {0: ("X", to_protein), 2: ("X", to_protein)}
        donor = {"X": list("MRVI-WYKTD")}               # K2R contact; L4I control; A5 deleted; R8K neither; E10D unresolved
        groups = {"X": {1: "contact", 3: "control", 4: "contact"}}
        sets, unmodelled = cs.plan(chains, mt_chain, host, donor, groups)
        assert sets["all"] == {(c, i): aa for c in (0, 2) for i, aa in ((0, "R"), (2, "I"), (6, "K"))}
        assert sets["contact"] == {(0, 0): "R", (2, 0): "R"} and sets["control"] == {(0, 2): "I", (2, 2): "I"}
        assert sets["none"] == {} and unmodelled == 4      # the deletion and the unresolved change, in each copy
        assert sets["gene:X"] == sets["all"] and set(sets) == {*cs.SETS, "gene:X"}

    def test_a_protein_alone_holds_its_own_substitutions_only(self):
        chains = [chain("A", "KVLAWYRT", (0, 0, 0)), chain("B", "GSTNDEQK", (0, 4, 0))]
        mt_chain = {0: ("X", {i: i + 1 for i in range(8)}), 1: ("Y", {i: i for i in range(8)})}
        sets, _ = cs.plan(chains, mt_chain, {"X": "MKVLAWYRTE", "Y": "GSTNDEQK"}, {"X": list("MRVLAWYRTE"), "Y": list("GSTNDEQR")}, {"X": {}, "Y": {}})
        assert sets["gene:X"] == {(0, 0): "R"} and sets["gene:Y"] == {(1, 7): "R"} and len(sets["all"]) == 2


class TestCollect:
    def table(self, tmp_path, host="homo-sapiens", donor="pan-paniscus", cx="IV", none=0.0, shift=0.0, fold_none=0.0):
        rows = []
        for w in cs.WEIGHTS:
            for seed in cs.SEEDS:
                for k, s in enumerate([*cs.SETS, *(("gene:COX1", "gene:COX2") if w == "stabddg" else ())]):
                    v = none if s == "none" else 1.0 + k + seed + shift + (10 if w != "stabddg" else 0)
                    rows.append({"host": host, "donor": donor, "complex": cx, "pdb": "5Z62", "chains": 9, "mt_chains": 3, "places": 900,
                                 "unmodelled": 2, "device": "cpu", "samples": 20, "weights": w, "seed": seed, "set": s,
                                 "substituted": max(3 - k, 0), "ddG": v, "fold": fold_none if s == "none" else 3 * v,
                                 "fold_one": v, "fold_two": v})
        d = tmp_path / "stab"; d.mkdir(exist_ok=True)
        pd.DataFrame(rows).to_csv(d / f"{host}.{donor}.{cx}.tsv", sep="\t", index=False)

    def test_nothing_to_collect_is_refused(self, tmp_path):
        with pytest.raises(SystemExit):
            cs.collect(tmp_path)

    def test_means_spread_and_what_is_missing(self, tmp_path, capsys):
        self.table(tmp_path)
        res = cs.collect(tmp_path)
        assert "NOT ASSESSED: 43 of 44" in capsys.readouterr().out
        r = res.iloc[0]
        assert len(res) == 1 and r.seeds == 3 and r.n_all == 3 and r.n_contact == 2 and r.n_control == 1
        assert r.stab_all == pytest.approx(2.0) and r.stab_contact == pytest.approx(3.0) and r.stab_control == pytest.approx(4.0)
        assert r.mpnn_all == pytest.approx(12.0) and r.stabfold_control == pytest.approx(14.0)
        assert r.stab_all_sd == pytest.approx(1.0)
        by = pd.read_csv(tmp_path / "stab_by_seed.tsv", sep="\t")
        assert list(by.seed) == [0, 1, 2] and list(by.stab_contact) == [2.0, 3.0, 4.0] and "stab_none" not in by.columns
        assert r.stab_fold_all == pytest.approx(6.0) and r.mpnn_fold_all == pytest.approx(36.0) and list(by.stab_fold_all) == [3.0, 6.0, 9.0]
        assert not [c for c in res.columns if "gene" in c or "COX" in c]
        g = pd.read_csv(tmp_path / "stab_genes.tsv", sep="\t")
        assert list(g.gene) == ["COX1", "COX2"] and list(g.seeds) == [3, 3]
        assert list(g.stab_bind) == [6.0, 7.0] and list(g.stab_fold) == [18.0, 21.0] and list(g.stab_bind_sd) == [1.0, 1.0]

    def test_an_empty_set_that_is_not_zero_is_refused(self, tmp_path):
        self.table(tmp_path, none=1e-6)
        with pytest.raises(SystemExit):
            cs.collect(tmp_path)
        self.table(tmp_path, fold_none=1e-6)
        with pytest.raises(SystemExit):
            cs.collect(tmp_path)


class TestVerdict:
    def summary(self, t1, t2, t3):
        rows = [{"test": "T1", "score": "s", "value": t1}, {"test": "T1", "score": "d_all", "value": -0.36},
                {"test": "T2", "score": "s", "value": t2}, {"test": "T2", "score": "d_all", "value": 0.74}]
        rows += [{"test": "T3", "score": "s", "value": v} for v in t3] + [{"test": "T3", "score": "d_all", "value": 1.0}] * 6
        return pd.DataFrame(rows)

    def test_all_three_count_two_are_partial(self):
        v = cybrid.verdict(self.summary(-0.5, 0.8, [1, 1, 1, 1, 1, 0]), ("s",)).iloc[0]
        assert v.verdict == "counts" and v.T3_named == 5 and v.T3_pairs == 6
        assert cybrid.verdict(self.summary(-0.5, 0.8, [1, 1, 1, 1, 0, 0]), ("s",)).iloc[0].verdict == "partial"
        assert cybrid.verdict(self.summary(-0.36, 0.8, [1] * 6), ("s",)).iloc[0].verdict == "partial"      # equal does not beat
        assert cybrid.verdict(self.summary(-0.5, 0.74, [1] * 6), ("s",)).iloc[0].verdict == "partial"
        assert cybrid.verdict(self.summary(-0.1, 0.5, [1] * 6), ("s",)).iloc[0].verdict == "does not count"

    def test_a_missing_readout_is_not_assessed(self):
        assert cybrid.verdict(self.summary(np.nan, 0.9, [1] * 6), ("s",)).iloc[0].verdict == "not assessed"
        assert cybrid.verdict(self.summary(-0.5, 0.9, [1, 1, 1, 1, 1, np.nan]), ("s",)).iloc[0].verdict == "not assessed"
        assert cybrid.verdict(self.summary(-0.5, 0.9, []), ("s",)).iloc[0].verdict == "not assessed"


def test_the_prediction_is_the_mean_of_the_draws_with_the_sign_turned():
    torch = pytest.importorskip("torch")

    class Fixed:
        device = "cpu"
        draws = iter([-5.0, -1.0, -2.0, -9.0, -3.0, -2.0])       # whole, one, two; twice. The model's sign: negative is worse

        def folding_ddG(self, d, s):
            return torch.tensor([next(self.draws)])

    s = np.zeros(3, dtype=np.int64)
    v = cs.predict(Fixed(), ({}, {}, {}), (s, s, s), samples=2)
    assert v == pytest.approx({"fold": 7.0, "fold_one": 2.0, "fold_two": 2.0, "ddG": 3.0})


class Counting:
    """Stands in for the model, in the model's own sign (negative is worse): for the whole system minus twice the
    places changed, for a partner alone minus the places changed. A substitution in an mtDNA-encoded chain is seen in
    the system and in partner one, so it counts 2 in `fold` and 1 in the binding prediction."""
    device = "cpu"

    def __init__(self):
        self.calls = []

    def folding_ddG(self, d, s):
        import torch
        n = int((torch.from_numpy(cs.encode(d["seq"], {})) != s[0]).sum())
        self.calls.append((d["name"], n))
        return torch.tensor([-(1.0 if d["name"].endswith(("_mt", "_nuc")) else 2.0) * n])


def test_one_task_on_the_real_structure_with_a_counting_model(tmp_path):
    pytest.importorskip("torch")
    from mitoage_pipelines.mitoport.complexes import COMPLEXES
    cif, gb, sites = paths.PDB / "5Z62.cif", paths.CYBRID_GB / "records.gb", COMPLEXES["IV"].work / "sites.tsv"
    if not (cif.is_file() and gb.is_file() and sites.is_file()):
        pytest.skip("needs `mitoage download pdb 5Z62`, `mitoage download cybrid` and `mitoage mitoport burial`")
    models = {}
    res = cs.run("homo-sapiens", "pan-paniscus", "IV", gb, paths.PDB, lambda c: COMPLEXES[c].work, tmp_path, tmp_path,
                 seeds=(0, 1), samples=2, model_for=lambda name: models.setdefault(name, Counting()))
    genes = sorted(s for s in set(res.set) if s.startswith("gene:"))
    assert genes and set(genes) <= {"gene:COX1", "gene:COX2", "gene:COX3"}
    assert len(res) == 2 * (len(cs.WEIGHTS) * len(cs.SETS) + len(genes)) and set(res.weights) == set(cs.WEIGHTS)
    assert set(res.set[res.weights != "stabddg"]) == set(cs.SETS)
    one = res[(res.weights == "stabddg") & (res.seed == 0)].set_index("set")
    # every substitution sits in an mtDNA-encoded chain: seen in the system and in partner one, never in partner two
    assert (one.ddG == one.substituted).all() and (one.fold == 2 * one.substituted).all() and (one.fold_two == 0).all()
    assert one.ddG["none"] == 0 and one.substituted["all"] > 0
    assert one.substituted["contact"] + one.substituted["control"] <= one.substituted["all"]
    assert one.substituted[genes].sum() == one.substituted["all"]
    assert all(n == 0 for m in models.values() for name, n in m.calls if name.endswith("_nuc"))
    assert res.mt_chains.iloc[0] == 3 and res.chains.iloc[0] > 3
    again = pd.read_csv(tmp_path / "stab" / "homo-sapiens.pan-paniscus.IV.tsv", sep="\t")
    assert len(again) == len(res)
    scores = pd.read_csv(paths.CYBRID_OUT / "scores.tsv", sep="\t") if (paths.CYBRID_OUT / "scores.tsv").is_file() else None
    if scores is not None:      # H16 counts differences over the whole proteins; the structure cannot hold more
        d_all = scores[(scores.host == "homo-sapiens") & (scores.donor == "pan-paniscus") & (scores.complex == "IV")].d_all.iloc[0]
        assert one.substituted["all"] + res.unmodelled.iloc[0] >= d_all >= one.substituted["all"]


def real_model():
    pytest.importorskip("torch"); pytest.importorskip("stabddg")
    w = paths.STABDDG / cs.WEIGHTS["stabddg"][0]
    if not w.is_file():
        pytest.skip("needs `mitoage download stabddg`")
    assert cs.sha256(w) == cs.WEIGHTS["stabddg"][1]
    return cs.load_model(w, "cpu")


def test_the_model_gives_zero_for_no_substitution_and_repeats_with_its_seed():
    import math
    model = real_model()
    import torch
    rng = np.random.default_rng(0)

    def helix(cid, seq, shift):       # an ideal alpha helix: CA at radius 2.3 A, 100 degrees and 1.5 A a residue
        n = len(seq); t = np.radians(100) * np.arange(n)
        ca = np.c_[2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)] + shift
        bb = np.stack([ca + rng.normal(0, 0.4, (n, 3)) + o for o in ([-1.2, 0, -0.6], [0, 0, 0], [1.2, 0, 0.6], [1.6, 1.0, 0.9])], axis=1)
        return {"id": cid, "seq": seq, "num": list(range(1, n + 1)), "bb": bb, "heavy": bb.reshape(-1, 3)}

    chains = [helix("A", P, np.zeros(3)), helix("B", Q, np.array([9.0, 0, 0]))]
    parts = [cs.structure("t", chains, [0, 1]), cs.structure("t_a", chains, [0]), cs.structure("t_b", chains, [1])]
    structs = tuple(d for d, _ in parts)

    def score(subs, seed):
        torch.manual_seed(seed)
        return cs.predict(model, structs, tuple(cs.encode(d["seq"], {at[k]: aa for k, aa in subs.items() if k in at}) for d, at in parts), samples=3)

    assert score({}, 0) == {"ddG": 0.0, "fold": 0.0, "fold_one": 0.0, "fold_two": 0.0}
    a = score({(0, 5): "D", (0, 9): "G"}, 0)
    assert math.isfinite(a["ddG"]) and a["ddG"] != 0.0 and a["fold"] != 0.0 and a["fold_two"] == 0.0
    assert score({(0, 5): "D", (0, 9): "G"}, 0) == a
    assert score({(0, 5): "D", (0, 9): "G"}, 1) != a
    # the three terms asked in the model's own order give the model's own binding prediction, draw for draw
    torch.manual_seed(0)
    seqs = [torch.from_numpy(cs.encode(d["seq"], {at[k]: aa for k, aa in {(0, 5): "D", (0, 9): "G"}.items() if k in at}))[None, :] for d, at in parts]
    with torch.no_grad():
        own = -float(np.mean([float(model(*structs, *seqs).squeeze()) for _ in range(3)]))
    assert own == pytest.approx(a["ddG"], abs=1e-5)
