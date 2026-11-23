"""H34, part B: the per-position reading of ProteinMPNN against the repository's own function on a real structure,
what the reading sees, the model's input, the run's bookkeeping with a stand-in reader on the real structure, the
collection's arithmetic and the rule's readout."""
import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.mitoport import cybrid
from mitoage_pipelines.mitoport import cybrid_stab as cs
from mitoage_pipelines.mitoport import variant_mpnn as vm
from mitoage_pipelines.mitoport import variant_stab as vs

CIF, GB = paths.PDB / "5Z62.cif", paths.CYBRID_GB / "records.gb"


def have_tool() -> bool:
    return (paths.LIGANDMPNN / "model_utils.py").is_file() and all((paths.LIGANDMPNN / "model_params" / f).is_file() for f in vm.CHECKPOINTS)


@pytest.fixture(scope="module")
def two_chains():
    """Complex IV's COX2 (its first 120 residues) and the smallest nuclear chain of H29's system around it, from the real structure."""
    if not (CIF.is_file() and GB.is_file()):
        pytest.skip("needs `mitoage download pdb 5Z62` and `mitoage download cybrid`")
    from mitoage_pipelines.mitoport.complexes import COMPLEXES
    hp = cybrid.species_proteins(cybrid.read_proteins(GB), cybrid.HUMAN)
    genes = tuple(COMPLEXES["IV"].mt)
    chains, mt = cs.system(cs.read_chains(CIF), {g: hp[g][0] for g in genes}, {g: p for g, (p, _, _) in hp.items() if g not in genes})
    cox2 = next(k for k, (g, _) in mt.items() if g == "COX2")
    small = min((k for k in range(len(chains)) if k not in mt), key=lambda k: len(chains[k]["seq"]))
    # COX2's first 120 residues: the repository's own function costs the fourth power of the length
    chains[cox2] = {**chains[cox2], **{k: chains[cox2][k][:120] for k in ("seq", "num", "bb")}}
    return chains, [cox2, small]


@pytest.fixture(scope="module")
def models():
    pytest.importorskip("torch")
    if not have_tool():
        pytest.skip("needs the LigandMPNN repository with its weights in data/input/ligandmpnn")
    return {arm: vm.load(paths.LIGANDMPNN, arm) for arm in vm.ARMS}


@pytest.mark.parametrize("arm", list(vm.ARMS))
def test_the_reading_is_the_repositorys_own_on_a_real_structure(two_chains, models, arm):
    import torch
    chains, which = two_chains
    feats, _ = vm.features(chains, which, vm.ARMS[arm][2])
    L = feats["S"].shape[1]
    noise = vm.order_noise(L, 7)
    with torch.no_grad():
        tool = models[arm].single_aa_score({**feats, "randn": noise[None], "batch_size": 1}, use_sequence=False)
    mine = vm.read(models[arm], feats, list(range(L)), noise, batch=5)
    last = (tool["decoding_order"][0, torch.arange(L), -1] == torch.arange(L)).numpy()
    assert last.mean() > 0.98             # the repository puts the scored place last unless its draw is near 0
    assert np.abs(tool["log_probs"][0].numpy()[last] - mine[last]).max() < 1e-4
    some = [3, 50, L - 1]                 # the positions read alone, in another batch, give the same numbers
    assert np.abs(vm.read(models[arm], feats, some, noise, batch=2) - mine[some]).max() < 1e-5


def test_the_scored_place_sees_every_other_residue_and_not_its_own(two_chains, models):
    import torch
    chains, which = two_chains
    model = models["memb"]
    feats, _ = vm.features(chains, which, 1)
    L = feats["S"].shape[1]
    noise = vm.order_noise(L, 0)
    p = 40
    base = vm.read(model, feats, [p], noise)
    own = {**feats, "S": feats["S"].clone()}
    own["S"][0, p] = (own["S"][0, p] + 1) % 20
    assert np.abs(vm.read(model, own, [p], noise) - base).max() < 1e-6
    _, _, E_idx = model.encode(feats)
    q = int(next(j for j in E_idx[0, p].tolist() if j != p))            # a neighbour in the graph
    other = {**feats, "S": feats["S"].clone()}
    other["S"][0, q] = (other["S"][0, q] + 1) % 20
    assert np.abs(vm.read(model, other, [p], noise) - base).max() > 1e-3
    assert np.abs(vm.read(model, feats, [p], vm.order_noise(L, 1)) - base).max() > 1e-6       # the order of the others matters
    assert np.array_equal(vm.read(model, feats, [p], vm.order_noise(L, 0)), base)
    assert torch.equal(vm.order_noise(L, 0), vm.order_noise(L, 0)) and not torch.equal(vm.order_noise(L, 0), vm.order_noise(L, 1))


def test_the_label_changes_the_membrane_model_only(two_chains, models):
    chains, which = two_chains
    f1, _ = vm.features(chains, which, 1)
    f0, _ = vm.features(chains, which, 0)
    noise = vm.order_noise(f1["S"].shape[1], 0)
    assert np.abs(vm.read(models["memb"], f1, [10, 60], noise) - vm.read(models["memb"], f0, [10, 60], noise)).max() > 1e-3


def test_the_models_input():
    pytest.importorskip("torch")
    nan = np.full(3, np.nan)
    bb = lambda n: np.arange(n * 12, dtype=float).reshape(n, 4, 3)
    a = {"seq": "MKZ", "num": [5, 6, 9], "bb": bb(3)}
    a["bb"][1, 3] = nan                     # residue 6 lacks its O
    b = {"seq": "GA", "num": [100, 100], "bb": bb(2)}     # an insertion code: two residues numbered 100
    feats, at = vm.features([a, {"seq": "W", "num": [1], "bb": bb(1)}, b], [0, 2], 1)
    assert at == {(0, 0): 0, (0, 1): 1, (0, 2): 2, (2, 0): 3, (2, 1): 4}
    assert feats["mask"][0].tolist() == [1, 0, 1, 1, 1]
    assert feats["X"][0, 1, 3].tolist() == [0.0, 0.0, 0.0] and feats["X"][0, 1, 2].tolist() == a["bb"][1, 2].tolist()
    assert feats["R_idx"][0].tolist() == [5, 6, 9, 100, 101] and feats["chain_labels"][0].tolist() == [0, 0, 0, 1, 1]
    assert feats["S"][0].tolist() == [cs.ALPHABET.index(c) for c in "MKXGA"]
    assert feats["membrane_per_residue_labels"][0].tolist() == [1] * 5
    assert vm.features([a], [0], 0)[0]["membrane_per_residue_labels"][0].tolist() == [0] * 3
    assert "membrane_per_residue_labels" not in vm.features([a], [0], None)[0]
    assert vm.renumber([1, 2, 2, 2, 3, 7, 7]) == [1, 2, 3, 4, 5, 9, 10]


def test_a_file_that_is_not_the_pinned_one_is_refused(tmp_path):
    pytest.importorskip("torch")
    (tmp_path / "model_utils.py").write_text("import torch\n")
    with pytest.raises(SystemExit):
        vm.tool_module(tmp_path)
    if not have_tool():
        pytest.skip("needs the LigandMPNN repository")
    import shutil
    shutil.copy(paths.LIGANDMPNN / "model_utils.py", tmp_path / "model_utils.py")
    (tmp_path / "model_params").mkdir()
    for f in vm.CHECKPOINTS:
        (tmp_path / "model_params" / f).write_bytes(b"not the weights")
    with pytest.raises(SystemExit):
        vm.load(tmp_path, "memb")


class Letters:
    """Stands in for the reading: the log-probability of letter a at input place i is -(1000 i + a), plus the seed's
    first draw, so each number says which place and letter it was taken from."""

    def __call__(self, model, feats, positions, noise, batch):
        self.S = feats["S"][0].numpy()
        return np.array([[-(1000.0 * p + a) + float(noise[0]) for a in range(21)] for p in positions])


def test_one_complex_on_the_real_structure_with_a_stand_in_reader(tmp_path):
    pytest.importorskip("torch")
    if not ((paths.PDB / "8H9S.cif").is_file() and GB.is_file() and (paths.VARIANTS / "clinvar_mt.json").is_file()):
        pytest.skip("needs `mitoage download pdb 8H9S`, `mitoage download cybrid` and `mitoage download variants`")
    reader = Letters()
    res = vm.run("V", GB, paths.PDB, paths.VARIANTS, tmp_path, tmp_path, seeds=(0, 1), model_for=lambda arm: None, reader=reader)
    v, _ = vs.variants(GB, paths.VARIANTS / "clinvar_mt.json", paths.VARIANTS / "gnomad.chrM.tsv")
    mine = v[v.gene.isin(["ATP6", "ATP8"])]      # complex V: some of its variants cannot be placed
    assert set(res.id) == set(mine.id) and (res.protein_pos > 1).all()
    sc = res[res.status == "scored"]
    assert len(sc) and sc.groupby("id").size().eq(len(vm.ARMS) * 2 * sc.groupby("id").copies.first()).all()
    a = np.array([cs.ALPHABET.index(x) for x in sc.ref_aa]); b = np.array([cs.ALPHABET.index(x) for x in sc.alt_aa])
    noise0 = np.array([float(vm.order_noise(len(reader.S), s)[0]) for s in sc.seed])
    place = (-(sc.lp_ref.to_numpy() - noise0) - a) / 1000.0
    assert np.allclose(place, np.round(place)) and (reader.S[np.round(place).astype(int)] == a).all()       # read where the reference residue is
    assert np.allclose(sc.lp_alt.to_numpy() - noise0, -(1000 * place + b))
    assert np.allclose(sc.score, sc.lp_ref - sc.lp_alt) and np.allclose(sc.score, b - a)
    left = res[res.status != "scored"]
    assert len(left) and left.score.isna().all() and len(left) == mine.id.nunique() - sc.id.nunique()


def test_the_score_is_the_mean_over_copies_then_seeds():
    row = dict(gene="CYTB", protein_pos=5, ref_aa="L", alt_aa="P", label="harmful", pos=1, ref="T", alt="C", af_hom=np.nan, complex="III", status="scored")
    long = pd.DataFrame([{**row, "id": "x", "arm": arm, "seed": s, "chain": c, "copies": 2, "score": float(10 * s + c + k)}
                         for k, arm in enumerate(vm.ARMS) for s in (0, 1, 2) for c in (0, 1)]
                        + [{**row, "id": "y", "status": "not assessed: the position is not resolved", "arm": "", "seed": -1, "chain": -1, "copies": 0, "score": np.nan}])
    r = vm.scores(long).set_index("id")
    # seed s: copies 10s + k and 10s + 1 + k, mean 10s + k + 0.5; over seeds 0, 1, 2: 10 + k + 0.5
    assert [r.loc["x", arm] for arm in vm.ARMS] == pytest.approx([10.5, 11.5, 12.5, 13.5])
    assert r.loc["x", "memb_sd"] == pytest.approx(10.0) and r.loc["x", "seeds"] == 3
    assert np.isnan(r.loc["y", "memb"]) and r.loc["y", "status"].startswith("not assessed")
    with pytest.raises(SystemExit):              # a change of the start codon in the tables: they predate the correction
        vm.scores(long.assign(protein_pos=lambda t: np.where(t.id == "y", 1, t.protein_pos)))


class TestReadout:
    def table(self, tmp_path, harmful, harmless, n_h=1, soluble=None):
        rows = []
        for k, x in enumerate(harmful):
            rows.append(dict(id=f"h{k}", gene="ND1", complex="I" if k % 2 else "IV", label="harmful", status="scored", memb=x, memb0=x, pmpnn=x,
                             smpnn=x if soluble is None else soluble[0][k], blosum=-2.0, n_aa=n_h, contact=k % 2 == 0))
        for k, x in enumerate(harmless):
            rows.append(dict(id=f"b{k}", gene="ND1", complex="I" if k % 2 else "IV", label="harmless", status="scored", memb=x, memb0=x, pmpnn=x,
                             smpnn=x if soluble is None else soluble[1][k], blosum=float(k % 5 - 2), n_aa=1 + k % 6, contact=k % 2 == 0))
        rows.append(dict(id="u", gene="ND1", complex="I", label="harmful", status="not assessed: the position is not resolved", blosum=0.0, n_aa=2, contact=False))
        pd.DataFrame(rows).to_csv(tmp_path / "mpnn_scores.tsv", sep="\t", index=False)

    def rule(self, tmp_path):
        r = vm.readout(tmp_path)
        return r[r.what == "rule"].reads.iloc[0].split(":")[0], r

    def test_a_score_that_separates_and_beats_the_baseline_counts(self, tmp_path):
        rng = np.random.default_rng(1)
        self.table(tmp_path, 5 + rng.normal(0, 1, 30), rng.normal(0, 1, 60), n_h=3)
        verdict, r = self.rule(tmp_path)
        assert verdict == "counts"
        assert r[(r.what == "area") & (r.score == "memb")].value.iloc[0] > 0.95
        assert r.iloc[0].n_harmful == 31 and "1 not assessed" in r.iloc[0].reads
        assert set(r[r.what.str.startswith("area (reported")].group) == {"all", "contact", "not contact", "complex I", "complex IV"}

    def test_a_good_score_the_baseline_matches_does_not_count(self, tmp_path):
        rng = np.random.default_rng(1)
        self.table(tmp_path, 5 + rng.normal(0, 1, 30), rng.normal(0, 1, 60), n_h=0)        # the site's variability alone separates them
        verdict, r = self.rule(tmp_path)
        assert verdict == "does not count"
        assert r[(r.what == "area") & (r.score == vm.BASELINE)].value.iloc[0] == 1.0

    def test_a_weak_score_does_not_count_even_above_the_baseline(self, tmp_path):
        rng = np.random.default_rng(2)
        self.table(tmp_path, 0.3 + rng.normal(0, 1, 40), rng.normal(0, 1, 80), n_h=7)      # the baseline points the wrong way
        verdict, r = self.rule(tmp_path)
        assert r[(r.what == "area") & (r.score == "memb")].value.iloc[0] < vm.AUC_MIN and verdict == "does not count"
        assert float(r[r.what == "difference"].reads.iloc[0].split("interval ")[1].split(" to")[0]) > 0

    def test_the_difference_to_the_soluble_weights_is_read_the_right_way(self, tmp_path):
        rng = np.random.default_rng(3)
        h, b = 5 + rng.normal(0, 1, 30), rng.normal(0, 1, 60)
        self.table(tmp_path, h, b, n_h=3, soluble=(rng.normal(0, 1, 30), rng.normal(0, 1, 60)))
        r = vm.readout(tmp_path)
        d = r[r.score == "memb - smpnn"].iloc[0]
        assert d.value > 0.3 and float(d.reads.split("interval ")[1].split(" to")[0]) > 0

    def test_one_kind_only_is_not_assessed(self, tmp_path):
        self.table(tmp_path, [], [0.1, 0.2, 0.3])
        pd.read_csv(tmp_path / "mpnn_scores.tsv", sep="\t").query("id != 'u'").to_csv(tmp_path / "mpnn_scores.tsv", sep="\t", index=False)
        assert self.rule(tmp_path)[0] == "not assessed"
