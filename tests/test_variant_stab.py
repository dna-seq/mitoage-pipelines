"""H29: which variants are kept and what they change (by hand, on a small made-up genome), how they are dealt out
and placed on the chains, the scoring's arithmetic with a stand-in for the model, and the rule's readout."""
import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.mitoport import variant_stab as vs

#        1   4   7   10  13  16  19  22  25  28
SEQ = "CCATGAAACTTTGGTAAGGCTACCGTTTCATGG"
# plus strand, positions 3 to 17: ATG AAA CTT TGG TAA -> M K L W *
# minus strand, positions 20 to 31: reverse complement of TACCGTTTCATG = CAT GAA ACG GTA... read as ATG AAA CGG TA.
CDS = [("ND1", 2, 17, 1), ("ND6", 19, 31, -1)]


class TestConsequence:
    def test_plus_strand(self):
        assert vs.consequence(SEQ, CDS, 6, "A", "C") == [("missense", "ND1", 2, "K", "Q")]       # AAA to CAA
        assert vs.consequence(SEQ, CDS, 8, "A", "G") == [("synonymous", "ND1", 2, "K", "K")]     # AAA to AAG
        assert vs.consequence(SEQ, CDS, 14, "G", "A") == [("synonymous", "ND1", 4, "W", "W")]    # TGA is tryptophan in this code
        assert vs.consequence(SEQ, CDS, 13, "G", "A") == [("nonsense", "ND1", 4, "W", "*")]      # TGG to TAG
        assert vs.consequence(SEQ, CDS, 16, "A", "C") == [("stop lost", "ND1", 5, "*", "S")]

    def test_a_change_in_the_first_codon_is_a_change_of_the_start(self):
        # the first codon of the plus-strand gene is positions 3 to 5 (ATG); of the minus-strand gene 29 to 31
        assert vs.consequence(SEQ, CDS, 4, "T", "C") == [("start", "ND1", 1, "M", "T")]          # ATG to ACG
        assert vs.consequence(SEQ, CDS, 5, "G", "A") == [("start", "ND1", 1, "M", "M")]          # ATG to ATA, still a start
        assert [c[:3] for c in vs.consequence(SEQ, CDS, 30, "A", "G")] == [("start", "ND6", 1)]
        assert vs.consequence(SEQ, CDS, 6, "A", "C")[0][0] == "missense"                         # the second codon is not

    def test_minus_strand(self):
        # the gene reads the reverse complement of positions 20 to 31: ATG AAA CGG TAG; position 31 is its first base
        assert vs.consequence(SEQ, CDS, 27, "T", "G") == [("missense", "ND6", 2, "K", "T")]      # codon AAA, second base to C

    def test_outside_the_genes_and_a_wrong_reference_base(self):
        assert vs.consequence(SEQ, CDS, 1, "C", "T") == [] and vs.consequence(SEQ, CDS, 19, "G", "A") == []
        with pytest.raises(ValueError):
            vs.consequence(SEQ, CDS, 6, "G", "C")

    def test_two_genes_that_overlap_each_give_their_change(self):
        both = vs.consequence(SEQ, [("ATP8", 2, 17, 1), ("ATP6", 0, 15, 1)], 6, "A", "C")     # in the second codon of each
        assert [(k, g) for k, g, *_ in both] == [("missense", "ATP8"), ("missense", "ATP6")] and both[0][2:] != both[1][2:]


def test_the_collection_leaves_out_the_changes_of_a_start():
    long = pd.DataFrame({"id": ["a", "a", "b", "c", "c"], "protein_pos": [1, 1, 2, 1, 1], "fold": [1.0, 2.0, 3.0, 4.0, 5.0]})
    kept, n = vs.without_start_changes(long)
    assert list(kept.id) == ["b"] and list(kept.fold) == [3.0] and n == 2
    same, none = vs.without_start_changes(long[long.protein_pos > 1])
    assert list(same.id) == ["b"] and none == 0


def test_what_counts_as_pathogenic():
    assert vs.is_pathogenic("Pathogenic") and vs.is_pathogenic("Likely pathogenic; drug response") and vs.is_pathogenic("Pathogenic/Likely pathogenic")
    assert not vs.is_pathogenic("Conflicting classifications of pathogenicity") and not vs.is_pathogenic("Benign") and not vs.is_pathogenic(None)
    assert not vs.is_pathogenic("Uncertain significance")


def test_the_variants_kept_and_what_is_left_out():
    cv = [dict(uid=1, title="", type="single nucleotide variant", spdi="NC:5:A:C", classification="Pathogenic", review=vs.PANEL),
          dict(uid=2, title="", type="single nucleotide variant", spdi="NC:7:A:G", classification="Pathogenic", review=vs.PANEL),          # synonymous
          dict(uid=3, title="", type="single nucleotide variant", spdi="NC:9:T:A", classification="Likely pathogenic", review="criteria provided, single submitter"),
          dict(uid=4, title="", type="Deletion", spdi="", classification="Pathogenic", review=vs.PANEL),
          dict(uid=5, title="", type="single nucleotide variant", spdi="NC:10:T:G", classification="Benign", review=vs.PANEL),
          dict(uid=6, title="", type="single nucleotide variant", spdi="", classification="Pathogenic", review=vs.PANEL)]
    g = pd.DataFrame({"position": [10, 12, 6, 27, 9, 26], "ref": list("TTATCT"), "alt": list("ACCGAG"), "filters": ["PASS", "PASS", "PASS", "npg", "PASS", "PASS"],
                      "AF_hom": [0.01, 0.0005, 0.5, 0.2, 0.02, 0.001]})
    # 12T>C is missense and rare; 27T>G does not pass gnomAD's filter
    v, counts = vs.classify(SEQ, CDS, cv, g)
    # 6A>C is the panel's; 10T>A is pathogenic for a single submitter: neither is harmful-by-the-panel twice nor harmless
    assert list(v.id) == ["m.6A>C:ND1", "m.9C>A:ND1", "m.26T>G:ND6"]
    assert list(v.label) == ["harmful", "harmless", "harmless"]
    assert counts == {"harmless left out: pathogenic in ClinVar": 2, "harmful without a sequence change given": 1}
    assert v.set_index("id").loc["m.9C>A:ND1", ["protein_pos", "ref_aa", "alt_aa", "af_hom"]].tolist() == [3, "L", "I", 0.02]


def test_a_complexs_variants_are_dealt_out_in_turn():
    v = pd.DataFrame({"id": [f"v{i}" for i in range(30)], "gene": ["CYTB"] * 25 + ["COX1"] * 5})
    parts = [vs.share(v, "III", k) for k in range(vs.CHUNKS["III"])]
    assert sorted(sum((list(p.id) for p in parts), [])) == sorted(v.id[:25]) and abs(len(parts[0]) - len(parts[1])) <= 1
    assert len(vs.tasks()) == sum(vs.CHUNKS.values()) and vs.tasks()[0] == ("I", 0)


def test_a_variant_is_placed_in_every_copy_or_said_why_not():
    chains = [{"seq": "KVLAW"}, {"seq": "GSTND"}, {"seq": "KVLAW"}, {"seq": "KVIAW"}]
    to_protein = {i: i + 1 for i in range(5)}
    row = lambda pos, ref, alt, gene="X": pd.Series({"gene": gene, "protein_pos": pos, "ref_aa": ref, "alt_aa": alt})
    two = {0: ("X", to_protein), 2: ("X", to_protein), 1: ("Y", {i: i for i in range(5)})}
    assert vs.place(row(4, "L", "F"), chains, two) == ({(0, 2): "F", (2, 2): "F"}, "scored")
    assert vs.place(row(9, "L", "F"), chains, two) == ({}, "not assessed: the position is not resolved")
    subs, status = vs.place(row(4, "L", "F"), chains, {0: ("X", to_protein), 3: ("X", to_protein)})
    assert subs == {} and status.startswith("not assessed: the structure holds I")


def test_the_scores_are_means_of_the_draws_with_the_sign_turned():
    torch = pytest.importorskip("torch")

    class Fixed:
        device = "cpu"
        draws = iter([[-4.0, -2.0], [-1.0, 0.0], [-6.0, -2.0], [-3.0, 0.0]])        # whole, one; twice. Two variants

        def folding_ddG(self, d, s):
            assert s.shape == (2, 3)
            return torch.tensor(next(self.draws))

    s = np.zeros((2, 3), dtype=np.int64)
    fold, ddg = vs.score_batch(Fixed(), {}, {}, s, s, samples=2)
    assert fold.tolist() == pytest.approx([5.0, 2.0]) and ddg.tolist() == pytest.approx([3.0, 2.0])


class TestReadout:
    def table(self, tmp_path, harmful, harmless, blosum_h=-2.0, n_h=1):
        rows = []
        for k, x in enumerate(harmful):
            rows.append(dict(id=f"h{k}", gene="ND1", label="harmful", status="scored", stab_fold=x, stab_bind=x / 2, mpnn_fold=x, stabfold_fold=x,
                             blosum=blosum_h, n_aa=n_h, contact=k % 2 == 0))
        for k, x in enumerate(harmless):
            rows.append(dict(id=f"b{k}", gene="ND1", label="harmless", status="scored", stab_fold=x, stab_bind=x / 2, mpnn_fold=x, stabfold_fold=x,
                             blosum=float(k % 5 - 2), n_aa=1 + k % 6, contact=k % 2 == 0))
        rows.append(dict(id="u", gene="ND1", label="harmful", status="not assessed: the position is not resolved", blosum=0.0, n_aa=2, contact=False))
        pd.DataFrame(rows).to_csv(tmp_path / "variant_scores.tsv", sep="\t", index=False)

    def rule(self, tmp_path):
        r = vs.readout(tmp_path)
        return r[r.what == "rule"].reads.iloc[0].split(":")[0], r

    def test_a_score_that_separates_and_beats_both_baselines_counts(self, tmp_path):
        rng = np.random.default_rng(1)
        self.table(tmp_path, 5 + rng.normal(0, 1, 30), rng.normal(0, 1, 60), blosum_h=0.0, n_h=3)
        verdict, r = self.rule(tmp_path)
        assert verdict == "counts"
        assert r[(r.what == "area") & (r.score == "stab_fold")].value.iloc[0] > 0.95
        first = r.iloc[0]
        assert first.n_harmful == 31 and first.n_harmless == 60 and "1 not assessed" in first.reads
        assert set(r[r.what.str.startswith("area (reported")].group) == {"ND1", "contact", "not contact"}

    def test_a_good_score_that_a_baseline_matches_does_not_count(self, tmp_path):
        rng = np.random.default_rng(1)
        self.table(tmp_path, 5 + rng.normal(0, 1, 30), rng.normal(0, 1, 60), blosum_h=-9.0)      # BLOSUM62 alone separates them
        assert self.rule(tmp_path)[0] == "does not count"

    def test_a_weak_score_does_not_count_even_above_the_baselines(self, tmp_path):
        rng = np.random.default_rng(2)
        self.table(tmp_path, 0.3 + rng.normal(0, 1, 40), rng.normal(0, 1, 80), blosum_h=3.0, n_h=7)      # both baselines point the wrong way
        verdict, r = self.rule(tmp_path)
        assert r[(r.what == "area") & (r.score == "stab_fold")].value.iloc[0] < vs.AUC_MIN and verdict == "does not count"
        assert (r[r.what == "difference"].reads.str.extract(r"interval (-?[\d.]+) to")[0].astype(float) > 0).all()

    def test_one_kind_only_is_not_assessed(self, tmp_path):
        self.table(tmp_path, [], [0.1, 0.2, 0.3])
        pd.read_csv(tmp_path / "variant_scores.tsv", sep="\t").query("id != 'u'").to_csv(tmp_path / "variant_scores.tsv", sep="\t", index=False)
        assert self.rule(tmp_path)[0] == "not assessed"


def test_the_difference_interval():
    rng = np.random.default_rng(0)
    h = np.r_[np.ones(40, bool), np.zeros(80, bool)]
    good = np.r_[3 + rng.normal(0, 1, 40), rng.normal(0, 1, 80)]
    noise = rng.normal(0, 1, 120)
    a, b, lo, hi = vs.compare(h, good, noise, n=300)
    assert a > 0.9 and 0.3 < b < 0.7 and lo > 0
    a, b, lo, hi = vs.compare(h, good, good, n=300)
    assert a == b and lo == 0 == hi
    with pytest.raises(ValueError):
        vs.compare(np.ones(5, bool), good[:5], noise[:5], n=3)
    from mitoage_pipelines.mitoport import cybrid
    tied = np.round(good)         # many ties: the area from ranks is the pairwise one
    assert vs.compare(h, tied, noise, n=5)[0] == pytest.approx(cybrid.auc(pd.Series(tied), pd.Series(h)))


class Counting:
    """Stands in for the model, in its own sign: for the whole system minus twice the places changed, for the
    mtDNA-encoded chains alone minus the places changed."""
    device = "cpu"

    def folding_ddG(self, d, s):
        import torch
        from mitoage_pipelines.mitoport import cybrid_stab as cs
        n = (torch.from_numpy(cs.encode(d["seq"], {}))[None, :] != s).sum(dim=1).float()
        return -(1.0 if d["name"].endswith("_mt") else 2.0) * n


def test_one_task_on_the_real_structure_with_a_counting_model(tmp_path):
    pytest.importorskip("torch")
    cif, gb = paths.PDB / "5Z62.cif", paths.CYBRID_GB / "records.gb"
    if not (cif.is_file() and gb.is_file() and (paths.VARIANTS / "clinvar_mt.json").is_file()):
        pytest.skip("needs `mitoage download pdb 5Z62`, `mitoage download cybrid` and `mitoage download variants`")
    res = vs.run("IV", 0, gb, paths.PDB, paths.VARIANTS, tmp_path, tmp_path, samples=2, model_for=lambda name: Counting())
    v, _ = vs.variants(gb, paths.VARIANTS / "clinvar_mt.json", paths.VARIANTS / "gnomad.chrM.tsv")
    mine = vs.share(v, "IV", 0)
    assert set(res.id) == set(mine.id) and set(res.gene) <= {"COX1", "COX2", "COX3"}
    sc = res[res.status == "scored"]
    assert len(sc) and (sc.copies >= 1).all() and (sc.fold == 2 * sc.copies).all() and (sc.ddG == sc.copies).all()
    assert sc.groupby("id").size().eq(sum(len(s) for s in vs.SEEDS.values())).all()
    assert res[res.status != "scored"].fold.isna().all() and len(res[res.status != "scored"]) == mine.id.nunique() - sc.id.nunique()
