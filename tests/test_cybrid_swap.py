"""H36, step 2: the swap's bookkeeping (where an ortholog's residues land, which species the null may draw, what a
restoration is, the candidate rule, the sums and the correlation) on frames small enough to count by hand and on the
real alignments and structures where they are installed; and, with the model's code (`uv sync --group stab`), that the
batched reading equals the repository's own `folding_ddG`."""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from mitoage_pipelines import paths
from mitoage_pipelines.mitoport import cybrid, cybrid_stab as cs, cybrid_swap as sw, donor_orthologs as do
from mitoage_pipelines.mitoport.complexes import COMPLEXES

HUMAN = "MKVLAWYRTEQNDSGHIFPC"


def chain(cid, seq, origin, num=None):
    n = len(seq)
    ca = np.array(origin, dtype=float) + np.c_[3.8 * np.arange(n), np.zeros(n), np.zeros(n)]
    bb = np.stack([ca + [-1, 0, 0], ca, ca + [1, 0, 0], ca + [1, 1, 0]], axis=1)
    return {"id": cid, "seq": seq, "num": list(num or range(1, n + 1)), "bb": bb, "heavy": bb.reshape(-1, 3), "n_heavy": [4] * n}


class TestOrthologs:
    def test_a_chain_is_placed_on_the_human_row_through_its_gaps(self):
        row = "--MKV--LAWYRTEQNDSGHIFPC-"
        cols = sw.chain_columns("KVLAWYRTE", row)             # the structure holds residues 2 to 10 of the protein
        assert [row[c] for c in cols] == list("KVLAWYRTE") and cols[0] == 3 and cols[2] == 7     # LA across the gap
        cols = sw.chain_columns("GGKVLAWYRTEQNDSGHIFPCGG", row)
        assert cols[0] is None and cols[-1] is None and [row[c] for c in cols[2:-2]] == list("KVLAWYRTEQNDSGHIFPC")

    def test_the_residues_that_change_are_the_species_and_the_rest_keep_the_human_one(self):
        seq = "KVLAWYRTE"
        cols = sw.chain_columns(seq, "--MKV--LAWYRTEQ")
        row = "--MKI--LAGYRXEQ"                                 # V to I (index 1); W to G (index 4); T unknown X (index 7)
        subs, placed, lost = sw.ortholog_subs(seq, cols, row)
        assert subs == {1: "I", 4: "G"} and placed == 8 and lost == 1 and placed + lost == len(seq)
        subs, placed, lost = sw.ortholog_subs(seq, cols, "--MK---LAWYRTEQ")      # a gap where the species lacks the residue
        assert subs == {} and lost == 1                                          # and the ortholog identical to human gives none
        assert sw.ortholog_subs(seq, cols, "--MKV--LAWYRTEQ") == ({}, 9, 0)

    def test_an_unmapped_residue_and_an_unknown_in_the_structure_are_not_placed(self):
        subs, placed, lost = sw.ortholog_subs("GKXAW", [None, 3, 4, 5, 6], "--MKVAWL")
        assert subs == {} and placed == 3 and lost == 2            # G has no column, X is no amino acid, only K A W are read
        assert sw.ortholog_subs("KAW", [3, 5, 6], "--MKVSWL")[0] == {1: "S"}

    def test_the_real_alignments_hold_the_structures_residues(self):
        f = paths.MITOPORT / "partners" / "aln" / "COX4I1.afa"
        if not (f.is_file() and (paths.PDB / "5Z62.cif").is_file()):
            pytest.skip("needs the vertebrate alignments of `mitoport iface-align` and 5Z62")
        rows = sw.read_alignment(f)
        ct = pd.read_csv(COMPLEXES["IV"].work / "chains.tsv", sep="\t", dtype={"chain": str})
        chains = [c for c in cs.read_chains(paths.PDB / "5Z62.cif") if c["id"] == ct.chain[ct.gene == "COX4I1"].iloc[0]]
        assert len(chains) == 1
        c = chains[0]
        cols = sw.chain_columns(c["seq"], rows["homo-sapiens"])
        assert all(col is not None for col in cols)
        subs, placed, lost = sw.ortholog_subs(c["seq"], cols, rows["homo-sapiens"])
        assert subs == {} and lost == 0 and placed == len(c["seq"])               # the human protein's own residues
        for s in ("gorilla-gorilla", "mus-musculus"):
            sub, p, l = sw.ortholog_subs(c["seq"], cols, rows[s])
            assert p + l == len(c["seq"]) and all(rows[s][cols[i]] == aa != c["seq"][i] for i, aa in sub.items())


    def test_a_chain_that_is_the_start_of_the_protein_keeps_to_the_start(self):
        f = paths.MITOPORT / "partners" / "aln" / "UQCRFS1.afa"
        if not (f.is_file() and COMPLEXES["III"].cif.is_file() and (COMPLEXES["III"].work / "chains.tsv").is_file()):
            pytest.skip("needs the vertebrate alignments of `mitoport iface-align` and 5XTE")
        ct = pd.read_csv(COMPLEXES["III"].work / "chains.tsv", sep="\t", dtype={"chain": str})
        ids = set(ct.chain[ct.gene == "UQCRFS1"])
        pre = [c for c in cs.read_chains(COMPLEXES["III"].cif) if c["id"] in ids and c["seq"].startswith("MLSVASRSG")]
        assert pre, "the cleaved presequence of the Rieske protein is in the structure"
        h = sw.read_alignment(f)["homo-sapiens"]
        ungapped = [c for c, ch in enumerate(h) if ch != "-"]
        # a chain of 57 residues that are the first 57 of the protein: each goes to its own column, the last too (an alignment
        # with end gaps that cost sent the last residue to the final glycine of the row)
        assert sw.chain_columns(pre[0]["seq"], h) == ungapped[:len(pre[0]["seq"])]

    def test_a_residue_placed_against_another_residue_is_not_placed(self):
        row = "MKVLAWYRTEQNDSGHIFPC"
        cols = sw.chain_columns("KVLAFYRTE", row)                 # the structure holds F where the human row has W
        assert cols[4] is None and [row[c] for i, c in enumerate(cols) if i != 4] == list("KVLAYRTE")


class TestContacts:
    def test_a_residue_touches_at_4_5_angstrom_and_not_beyond(self):
        mt = chain("M", "A", (0, 0, 0))                  # heavy atoms at (-1,0,0), (0,0,0), (1,0,0), (1,1,0)
        near = chain("N", "AAA", (0, 0, 0))
        near["heavy"] = np.array([[1000.0 + k, 0.0, 0.0] for k in range(12)])           # every atom far away ...
        near["heavy"][4] = [0.0, -4.5, 0.0]                                              # ... but residue 1's first: exactly 4.5 from (0,0,0)
        near["heavy"][8] = [0.0, -4.6, 0.0]                                              # and residue 2's: 4.6
        assert list(sw.contact_flags([mt, near], [0], [1])[1]) == [False, True, False]

    def test_every_copy_adds_its_columns(self):
        cols = {1: [4, 5, None, 7], 2: [4, 5, 6, 7]}
        flags = {1: np.array([True, False, True, False]), 2: np.array([False, False, True, True])}
        assert sw.contact_columns(cols, flags) == [4, 6, 7]              # a residue without a column adds none

    def test_a_gap_or_an_unknown_is_not_a_difference(self):
        human = "MKVLAW"
        assert sw.contact_differences("MRVL-W", human, [1, 2, 4]) == 1
        assert sw.contact_differences("MXVLAF", human, [1, 5]) == 1 and sw.contact_differences(human, human, [0, 1, 2]) == 0
        assert sw.contact_differences("MRVLAW", "M-VLAW", [1]) == 0          # nor is a column the human row lacks

    def test_the_real_structure_gives_contacts_in_complex_iv(self):
        cif = paths.PDB / "5Z62.cif"
        if not cif.is_file():
            pytest.skip("needs `mitoage download pdb 5Z62`")
        chains = cs.read_chains(cif)
        assert all(sum(c["n_heavy"]) == len(c["heavy"]) and len(c["n_heavy"]) == len(c["seq"]) for c in chains)
        mt = [i for i, c in enumerate(chains) if c["id"] in ("A", "B", "C")]
        other = [i for i in range(len(chains)) if i not in mt]
        flags = sw.contact_flags(chains, mt, other)
        assert any(f.any() for f in flags.values()) and not all(f.all() for f in flags.values())
        # the residues counted are the ones within 4.5 A, recomputed by brute force on one chain
        i = [k for k in other if flags[k].any()][0]
        d = np.linalg.norm(chains[i]["heavy"][:, None, :] - np.concatenate([chains[k]["heavy"] for k in mt])[None], axis=2).min(axis=1)
        edges = np.concatenate([[0], np.cumsum(chains[i]["n_heavy"])])
        assert list(flags[i]) == [bool((d[a:b] <= cs.CONTACT_A).any()) for a, b in zip(edges[:-1], edges[1:])]


class TestNull:
    rows = {"homo-sapiens": "MKVLAWYRTE", "donor": "RKVLAWYRTA",          # donor differs at 2 of the 4 contact columns below
            "m_zero": "MKVLAWYRTE", "m_one": "RKVLAWYRTE", "m_two": "RKVLAWYRTA", "m_two_b": "MKVLAWYQTA", "m_three": "RKVLCWYRTA",
            "m_gap": "R-VLAWYRT-", "bird_two": "RKVLAWYQTA", "m_four": "RQVLAWFRTA", "m_unk": "XKVLAWYRTX"}
    mammals = {s for s in rows if not s.startswith("bird")}

    def test_the_window_is_a_quarter_either_way_of_the_donors_differences(self):
        columns = [0, 6, 7, 9]
        human = self.rows["homo-sapiens"]
        got = {s: sw.contact_differences(r, human, columns) for s, r in self.rows.items()}
        assert got["donor"] == 2 and got["m_zero"] == 0 and got["m_one"] == 1 and got["m_gap"] == 1 and got["m_unk"] == 0
        assert got["m_two"] == got["m_two_b"] == got["m_three"] == 2 and got["m_four"] == 3
        cand, n = sw.null_candidates(self.rows, self.mammals, "donor", columns)
        assert n == 2 and cand == ["m_three", "m_two", "m_two_b"]                  # 1.5 <= k <= 2.5: the species with 2, not those with 1 or 3
        assert sw.null_candidates(self.rows, self.mammals, "donor", columns, window=(1, 3))[0] == ["m_gap", "m_one"]    # 0.5 <= k <= 1.5

    def test_the_exact_bounds_of_three_quarters_and_five_quarters(self):
        human = "A" * 8
        rows = {"homo-sapiens": human, "donor": "C" * 4 + "A" * 4}                               # n = 4 differences
        for k in range(8):
            rows[f"m{k}"] = "C" * k + "A" * (8 - k)
        cand, n = sw.null_candidates(rows, set(rows), "donor", list(range(8)))
        assert n == 4 and cand == ["m3", "m4", "m5"]                                            # 3 <= k <= 5, the donor itself left out

    def test_a_donor_without_a_difference_draws_only_species_without_one(self):
        cand, n = sw.null_candidates(self.rows, self.mammals, "m_zero", [0, 6, 7, 9])
        assert n == 0 and cand == ["m_unk"]                      # the host is not one; m_unk holds unknowns there, which are no difference
        rows = {**self.rows, "m_zero_b": "MKVLCWYRTE"}
        assert sw.null_candidates(rows, self.mammals | {"m_zero_b"}, "m_zero", [0, 6, 7, 9])[0] == ["m_unk", "m_zero_b"]

    def test_the_donor_the_host_and_non_mammals_are_never_candidates(self):
        cand, _ = sw.null_candidates(self.rows, self.mammals, "donor", [0, 6, 7, 9], window=(0, 100))
        assert "donor" not in cand and "homo-sapiens" not in cand and "bird_two" not in cand and "m_two" in cand

    def test_a_gap_in_the_species_is_no_difference_at_the_contact_columns(self):
        assert sw.contact_differences(self.rows["m_gap"], self.rows["homo-sapiens"], [1, 9]) == 0

    def test_draws_come_from_the_generator_seeded_with_the_draw_and_repeat(self):
        cand = [f"s{i}" for i in range(40)]
        got = [sw.draw_species(cand, k) for k in range(50)]
        assert got == [sw.draw_species(cand, k) for k in range(50)]
        assert got == [cand[int(np.random.default_rng((0, k)).integers(40))] for k in range(50)]
        assert len(set(got)) > 10 and sw.draw_species(cand, 3, seed=1) == cand[int(np.random.default_rng((1, 3)).integers(40))]
        assert sw.draw_species(["only"], 7) == "only"


class TestSizes:
    def test_the_largest_are_three_or_a_quarter_rounded_up_above_twelve(self):
        assert [sw.top_size(n) for n in (1, 3, 12, 13, 16, 17, 27)] == [3, 3, 3, 4, 4, 5, 7]


def frame(donors, complex_, gene, kinds):
    """The tasks' scores: `kinds` is {(donor, gene): {"B0": [per seed], "B1": [per seed], "null": [per draw: per seed]}}."""
    cases, scores = [], []
    for (d, g), v in kinds.items():
        cases.append({"host": "homo-sapiens", "donor": d, "complex": complex_, "gene": g, "status": "scored", "substituted": 1,
                      "contact_differences": 1, "candidates": 5})
        for seed, x in enumerate(v["B0"]):
            scores.append({"donor": d, "complex": complex_, "gene": g, "kind": "B0", "draw": -1, "seed": seed, "ddG": x})
        for seed, x in enumerate(v["B1"]):
            scores.append({"donor": d, "complex": complex_, "gene": g, "kind": "B1", "draw": -1, "seed": seed, "ddG": x})
        for k, per in enumerate(v.get("null", [])):
            for seed, x in enumerate(per):
                scores.append({"donor": d, "complex": complex_, "gene": g, "kind": "null", "draw": k, "seed": seed, "ddG": x})
    return pd.DataFrame(cases), pd.DataFrame(scores)


class TestRestoration:
    def test_restoration_is_b0_minus_b1_on_the_means_and_the_null_is_b0_minus_each_draws_mean(self):
        null = [[2.0, 4.0], [0.0, 2.0], [5.0, 5.0], [1.0, 1.0], [0.0, 0.0]]          # draw means 3, 1, 5, 1, 0
        cases, scores = frame(["d"], "IV", "A", {("d", "A"): {"B0": [10.0, 12.0], "B1": [6.0, 6.0], "null": null}})
        r = sw.restoration_table(cases, scores).iloc[0]
        assert r.B0 == 11.0 and r.B1 == 6.0 and r.restoration == 5.0
        drawn = np.array([11 - 3, 11 - 1, 11 - 5, 11 - 1, 11 - 0], float)
        assert r.null_mean == drawn.mean() and r.null_p95 == np.percentile(drawn, 95)
        assert r.restoration_seed_sd == pytest.approx(np.std([4.0, 6.0], ddof=1)) and r.seeds == 2 and r.null_draws == 5
        assert r.null_p == (1 + np.sum(drawn >= 5.0)) / 6 and not r.exceeds_null          # 5 is below the null's 95th percentile (10.8)
        assert np.percentile(drawn, 95) != np.percentile(drawn, 90)

    def test_equal_to_the_95th_percentile_does_not_exceed_it_and_ties_count_in_the_share(self):
        cases, scores = frame(["d"], "IV", "A", {("d", "A"): {"B0": [10.0], "B1": [4.0], "null": [[4.0]] * 4}})      # every draw restores 6, as the real one
        r = sw.restoration_table(cases, scores).iloc[0]
        assert r.restoration == 6.0 and r.null_p95 == 6.0 and not r.exceeds_null and r.null_p == 1.0

    def test_a_restoration_above_the_null_percentile_exceeds_and_a_missing_null_never_does(self):
        null = [[9.0], [8.0], [7.0], [9.0]]                                          # restorations 1, 2, 3, 1 against B0 10
        cases, scores = frame(["d"], "IV", "A", {("d", "A"): {"B0": [10.0], "B1": [4.0], "null": null}, ("d", "B"): {"B0": [10.0], "B1": [4.0]}})
        r = sw.restoration_table(cases, scores).set_index("gene")
        assert r.loc["A", "restoration"] == 6.0 and bool(r.loc["A", "exceeds_null"])
        assert r.loc["B", "null_draws"] == 0 and np.isnan(r.loc["B", "null_p95"]) and not r.loc["B", "exceeds_null"]

    def test_a_negative_restoration_is_below_a_null_of_zeros_and_never_among_the_largest(self):
        cases, scores = frame(["d"], "IV", "A", {("d", g): {"B0": [1.0], "B1": [b1], "null": [[1.0]] * 3} for g, b1 in (("A", 3.0), ("B", 1.0), ("C", 0.5))})
        r = sw.restoration_table(cases, scores).set_index("gene")
        assert list(r.restoration) == [0.5, 0.0, -2.0] and list(r.index) == ["C", "B", "A"]
        assert list(r.in_top) == [True, False, False]                                  # 0 is not among the largest, though only one is above it
        assert not r.loc["B", "exceeds_null"]                                          # 0 against a null of 0

    def test_rank_and_the_number_of_genes_ranked(self):
        spec = {("d", g): {"B0": [10.0], "B1": [10.0 - v], "null": [[10.0]]} for g, v in zip("ABCDE", (5, 4, 4, 3, 1))}
        r = sw.restoration_table(*frame(["d"], "IV", "A", spec)).set_index("gene")
        assert list(r["rank"]) == [1, 2, 2, 4, 5] and set(r.genes_ranked) == {5}
        assert list(r.in_top) == [True, True, True, False, False]                    # three largest; the tie shares rank 2, the next is 4th

    def test_more_than_twelve_genes_keep_a_quarter(self):
        spec = {("d", f"G{i:02d}"): {"B0": [20.0], "B1": [20.0 - (16 - i)], "null": [[20.0]]} for i in range(16)}
        r = sw.restoration_table(*frame(["d"], "I", "A", spec))
        assert int(r.in_top.sum()) == 4 and list(r[r.in_top].gene) == ["G00", "G01", "G02", "G03"]       # a quarter of 16

    def test_b0_and_b1_must_both_be_there_for_every_seed(self):
        cases, scores = frame(["d"], "IV", "A", {("d", "A"): {"B0": [1.0, 2.0], "B1": [1.0]}})
        with pytest.raises(SystemExit):
            sw.restoration_table(cases, scores)

    def test_a_gene_that_was_not_assessed_has_no_row(self):
        cases, scores = frame(["d"], "IV", "A", {("d", "A"): {"B0": [1.0], "B1": [0.0], "null": [[1.0]]}})
        cases = pd.concat([cases, pd.DataFrame([{"host": "homo-sapiens", "donor": "d", "complex": "IV", "gene": "Z", "status": "not assessed: no alignment of the gene"}])])
        assert list(sw.restoration_table(cases, scores).gene) == ["A"]


class TestCandidates:
    def spec(self, per_donor):
        """per_donor: {donor: {gene: (restoration, null draws' restorations)}}, with B0 = 10, one seed."""
        out = {}
        for d, genes in per_donor.items():
            for g, (rest, nulls) in genes.items():
                out[(d, g)] = {"B0": [10.0], "B1": [10.0 - rest], "null": [[10.0 - n] for n in nulls]}
        return out

    def run(self, per_donor, cases_extra=()):
        cases, scores = frame(list(per_donor), "IV", "A", self.spec(per_donor))
        if cases_extra:
            cases = pd.concat([cases, pd.DataFrame(list(cases_extra))], ignore_index=True)
        res = sw.restoration_table(cases, scores)
        return res, sw.candidate_table(cases, res).set_index("gene")

    def test_both_conditions_in_more_than_half_of_the_donors_with_an_ortholog(self):
        low = [0.1, 0.2, 0.3, 0.4]
        genes = lambda a, b: {"A": (a, low), "B": (b, low)}
        # three donors; A is the largest and above its null in two of three, B in one
        res, c = self.run({"d1": genes(5.0, 0.05), "d2": genes(5.0, 0.05), "d3": genes(0.05, 5.0)})
        assert c.loc["A", "both"] == 2 and bool(c.loc["A", "candidate"]) and c.loc["B", "both"] == 1 and not c.loc["B", "candidate"]
        assert c.loc["A", "donors_with_ortholog"] == 3 and c.loc["A", "in_top"] == 3

    def test_largest_but_not_above_the_null_is_not_a_candidate_and_neither_is_above_the_null_but_not_largest(self):
        wide = [0.0, 3.0, 6.0, 9.0]
        res, c = self.run({"d1": {"A": (5.0, wide)}, "d2": {"A": (5.0, wide)}})
        assert c.loc["A", "in_top"] == 2 and c.loc["A", "exceeds_null"] == 0 and not c.loc["A", "candidate"]
        genes = {g: (v, [0.0] * 4) for g, v in zip("ABCDE", (9, 8, 7, 6, 1))}                  # E is above its null of 0 but is the 5th
        res, c = self.run({"d1": genes, "d2": genes})
        assert c.loc["E", "exceeds_null"] == 2 and c.loc["E", "in_top"] == 0 and not c.loc["E", "candidate"] and bool(c.loc["A", "candidate"])

    def test_exactly_half_is_not_more_than_half_and_a_donor_without_an_ortholog_is_not_counted(self):
        low = [0.0] * 4
        res, c = self.run({"d1": {"A": (5.0, low)}, "d2": {"A": (0.0, low)}})
        assert c.loc["A", "both"] == 1 and c.loc["A", "donors_with_ortholog"] == 2 and not c.loc["A", "candidate"]
        extra = [{"host": "homo-sapiens", "donor": "d3", "complex": "IV", "gene": "A", "status": "not assessed: no ortholog of the donor in the alignment"}]
        res, c = self.run({"d1": {"A": (5.0, low)}, "d2": {"A": (0.0, low)}}, extra)
        assert c.loc["A", "donors"] == 3 and c.loc["A", "donors_with_ortholog"] == 2 and not c.loc["A", "candidate"]
        extra2 = extra + [{**extra[0], "donor": "d4"}]
        res, c = self.run({"d1": {"A": (5.0, low)}}, extra2)
        assert c.loc["A", "donors"] == 3 and c.loc["A", "both"] == 1 and bool(c.loc["A", "candidate"])     # one of the one donor that has an ortholog
        res, c = self.run({"d1": {"A": (5.0, low)}, "d2": {"A": (5.0, low)}}, extra)
        assert c.loc["A", "both"] == 2 and bool(c.loc["A", "candidate"]) and "d3: not assessed: no ortholog" in c.loc["A", "not_assessed"]

    def test_a_donor_whose_null_is_not_assessed_has_an_ortholog_and_is_not_above_it(self):
        cases, scores = frame(["d"], "IV", "A", {("d1", "A"): {"B0": [10.0], "B1": [5.0], "null": [[10.0]] * 3}, ("d2", "A"): {"B0": [10.0], "B1": [5.0]}})
        res = sw.restoration_table(cases, scores)
        c = sw.candidate_table(cases, res).iloc[0]
        assert c.donors_with_ortholog == 2 and c.both == 1 and not c.candidate

    def test_a_gene_nowhere_assessed_is_not_a_candidate(self):
        cases = pd.DataFrame([{"host": "homo-sapiens", "donor": "d1", "complex": "IV", "gene": "Z", "status": "not assessed: no alignment of the gene"}])
        c = sw.candidate_table(cases, sw.restoration_table(cases, pd.DataFrame(columns=["donor", "complex", "gene", "kind", "draw", "seed", "ddG"]))).iloc[0]
        assert c.donors_with_ortholog == 0 and not c.candidate and "no alignment" in c.not_assessed


class TestSums:
    def test_sums_over_the_genes_scored_for_every_donor_and_the_correlation(self):
        spec = {}
        vals = {"d1": {"A": 4.0, "B": 1.0, "C": 3.0}, "d2": {"A": 2.0, "B": 2.0}, "d3": {"A": 1.0, "B": 0.5, "C": 9.0}}
        for d, genes in vals.items():
            for g, v in genes.items():
                spec[(d, g)] = {"B0": [10.0], "B1": [10.0 - v], "null": [[10.0]]}
        cases, scores = frame(list(vals), "I", "A", spec)
        cases = pd.concat([cases, pd.DataFrame([{"host": "homo-sapiens", "donor": "d2", "complex": "I", "gene": "C", "status": "not assessed: no ortholog of the donor in the alignment"}])], ignore_index=True)
        res = sw.restoration_table(cases, scores)
        meas = pd.DataFrame({"host": "homo-sapiens", "donor": ["d1", "d2", "d3", "d1"], "complex": ["I", "I", "I", "IV"], "assay": "x", "outcome": "measured",
                             "pct_low": [59.0, 57.0, 55.0, np.nan], "pct_high": [59.0, 57.0, 55.0, np.nan], "significant": ["yes", "yes", "yes", "no"],
                             "source_doi": "x", "location": "x", "note": ""})
        meas["pct"] = (meas.pct_low + meas.pct_high) / 2
        sums, rho = sw.activity_table(res, cases, meas)
        s = sums.set_index("donor")
        assert list(s.genes_summed) == [2, 2, 2] and list(s.restoration_sum) == [5.0, 4.0, 1.5]       # C is left out for all: d2 has none
        assert list(s.pct) == [59.0, 57.0, 55.0] and list(s.genes_scored) == [3, 2, 3]
        assert rho.set_index("complex").loc["I", "n"] == 3 and rho.set_index("complex").loc["I", "spearman"] == pytest.approx(1.0)
        assert rho.set_index("complex").loc["pooled", "spearman"] == pytest.approx(stats.spearmanr([59, 57, 55], [5.0, 4.0, 1.5])[0])

    def test_fewer_than_three_cases_or_a_constant_gives_no_correlation(self):
        cases, scores = frame(["d"], "IV", "A", {("d1", "A"): {"B0": [10.0], "B1": [8.0], "null": [[10.0]]}, ("d2", "A"): {"B0": [10.0], "B1": [7.0], "null": [[10.0]]}})
        meas = pd.DataFrame({"host": "homo-sapiens", "donor": ["d1", "d2"], "complex": "IV", "assay": "x", "outcome": "measured", "pct_low": [50.0, 60.0],
                             "pct_high": [50.0, 60.0], "significant": "yes", "source_doi": "x", "location": "x", "note": ""})
        meas["pct"] = meas.pct_low
        _, rho = sw.activity_table(sw.restoration_table(cases, scores), cases, meas)
        assert rho.n.tolist() == [2, 2] and rho.spearman.isna().all()


class TestGeneChains:
    def test_nuclear_chains_are_grouped_by_gene_and_the_mitochondrial_ones_left_out(self):
        chains = [{"id": "A"}, {"id": "B"}, {"id": "C"}, {"id": "D"}, {"id": "E"}]
        got = sw.gene_chains(chains, {0: ("ND1", {}), 4: ("ND1", {})}, {"A": "ND1", "B": "U1", "C": "U2", "D": "U1", "E": "ND1"})
        assert got == {"U1": [1, 3], "U2": [2]}
        with pytest.raises(SystemExit):
            sw.gene_chains(chains, {0: ("ND1", {})}, {"A": "ND1", "B": "U1", "C": "U2", "E": "ND1"})          # D has no gene


class TestTasks:
    def test_a_task_for_each_ape_donor_with_an_outcome_and_each_nuclear_gene_of_the_complex(self):
        t = sw.tasks()
        assert {d for d, _, _ in t} == {"pan-troglodytes", "pan-paniscus", "gorilla-gorilla"}          # the orangutans have only a whole-cybrid row
        for c in COMPLEXES:
            assert sorted(g for d, cx, g in t if cx == c and d == "gorilla-gorilla") == sorted(COMPLEXES[c].nuclear)
        assert len(t) == 3 * sum(len(c.nuclear) for c in COMPLEXES.values()) and t == sorted(t, key=lambda x: (x[1], x[0], COMPLEXES[x[1]].nuclear.index(x[2])))

    def test_a_donor_or_host_without_a_stated_outcome_is_not_a_task(self, tmp_path):
        f = tmp_path / "m.tsv"
        m = cybrid.read_measurements().copy()
        m = pd.concat([m, m[m.donor == "gorilla-gorilla"].assign(host="mus-musculus", donor="mus-spretus")])
        m = pd.concat([m, m[(m.donor == "pan-troglodytes") & (m["complex"] == "IV")].assign(donor="pongo-abelii")])      # a donor with no ortholog
        m.drop(columns="pct").to_csv(f, sep="\t", index=False)
        assert {d for d, _, _ in sw.tasks(f)} == {"pan-troglodytes", "pan-paniscus", "gorilla-gorilla"}
        # the human host's row for a pair removed, another host's row for the same pair added: no task for it
        m2 = m[~((m.donor == "pan-paniscus") & (m["complex"] == "V"))]
        m2 = pd.concat([m2, m[(m.donor == "pan-paniscus") & (m["complex"] == "V")].assign(host="mus-musculus")])
        m2.drop(columns="pct").to_csv(f, sep="\t", index=False)
        assert ("pan-paniscus", "V", COMPLEXES["V"].nuclear[0]) not in sw.tasks(f) and ("pan-paniscus", "IV", "COX4I1") in sw.tasks(f)


class TestCollect:
    def test_nothing_to_collect_is_refused_and_a_missing_task_is_said(self, tmp_path, capsys):
        with pytest.raises(SystemExit):
            sw.collect(tmp_path)
        (tmp_path / "swap").mkdir()
        d, c, g = sw.tasks()[0]
        case = {"host": "homo-sapiens", "donor": d, "complex": c, "gene": g, "status": "scored", "substituted": 1, "contact_differences": 1, "candidates": 3}
        pd.DataFrame([case]).to_csv(tmp_path / "swap" / f"{d}.{c}.{g}.cases.tsv", sep="\t", index=False)
        rows = [{"host": "homo-sapiens", "donor": d, "complex": c, "gene": g, "kind": k, "draw": dr, "species": "s", "seed": s, "substituted": 1, "fold": 0.0,
                 "fold_one": 0.0, "fold_two": 0.0, "ddG": v} for s in (0, 1, 2) for k, dr, v in (("B0", -1, 4.0), ("B1", -1, 1.0), ("null", 0, 3.0), ("null", 1, 2.0), ("null", 2, 3.5))]
        pd.DataFrame(rows).to_csv(tmp_path / "swap" / f"{d}.{c}.{g}.tsv", sep="\t", index=False)
        res = sw.collect(tmp_path)
        assert f"MISSING: {len(sw.tasks()) - 1} of {len(sw.tasks())}" in capsys.readouterr().out
        assert len(res) == 1 and res.iloc[0].restoration == 3.0
        r = res.iloc[0]
        assert r.null_draws == 3 and r.null_p95 == np.percentile([1.0, 2.0, 0.5], 95)          # the rows of kind "null" are read as they are
        for name in ("scores", "cases", "restoration", "candidates", "sums", "spearman"):
            assert (tmp_path / f"swap_{name}.tsv").is_file()
        assert len(pd.read_csv(tmp_path / "swap_scores.tsv", sep="\t")) == 15
        assert not (tmp_path / "swap_unassessed.tsv").exists()           # the human host's tables are the first run's six


# --- the model and the structures ------------------------------------------------------------------------------------


def small_model():
    torch = pytest.importorskip("torch"); pytest.importorskip("stabddg")
    from stabddg.model import StaBddG
    from stabddg.mpnn_utils import ProteinMPNN
    torch.manual_seed(0)
    pm = ProteinMPNN(node_features=128, edge_features=128, hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3, k_neighbors=24,
                     dropout=0.0, augment_eps=0.0)
    return StaBddG(pmpnn=pm, use_antithetic_variates=True, noise_level=cs.NOISE, device="cpu").eval()


def helix_domain():
    rng = np.random.default_rng(1)
    def helix(cid, seq, shift):
        n = len(seq); t = np.radians(100) * np.arange(n)
        ca = np.c_[2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)] + shift
        bb = np.stack([ca + rng.normal(0, 0.3, (n, 3)) + o for o in ([-1.2, 0, -0.6], [0, 0, 0], [1.2, 0, 0.6], [1.6, 1.0, 0.9])], axis=1)
        return {"id": cid, "seq": seq, "num": list(range(1, n + 1)), "bb": bb, "heavy": bb.reshape(-1, 3), "n_heavy": [4] * n}
    chains = [helix("A", "MKVLAWYRTEQNDSGHIFPC" * 2, np.zeros(3)), helix("B", "GSTNDEQKRHGSTNDEQKRH", np.array([9.0, 0, 0]))]
    return cs.structure("t", chains, [0, 1])


class TestModel:
    def test_the_batched_reading_equals_the_repositorys_own_folding_ddG_sample_for_sample(self):
        import torch
        model = small_model()
        domain, at = helix_domain()
        subs = [{}, {(0, 5): "D"}, {(0, 5): "D", (1, 3): "G", (0, 30): "W"}, {(1, 7): "A"}, {(0, 2): "P", (0, 9): "K"}]
        seqs = [cs.encode(domain["seq"], {at[k]: aa for k, aa in s.items()}) for s in subs]
        samples = 3
        got = sw.ddg_rows(model, domain, seqs, seed=4, samples=samples, batch=2)
        torch.manual_seed(4)
        T = torch.from_numpy(np.stack(seqs))
        with torch.no_grad():
            own = -np.mean([model.folding_ddG(domain, T).reshape(-1).double().numpy() for _ in range(samples)], axis=0)
        assert got[0] == 0.0 and np.abs(got[1:]).min() > 0
        assert got == pytest.approx(own, abs=1e-4)

    def test_a_row_depends_on_its_inputs_and_the_seed_alone(self):
        model = small_model()
        domain, at = helix_domain()
        a = cs.encode(domain["seq"], {at[(0, 5)]: "D"}); b = cs.encode(domain["seq"], {at[(1, 3)]: "G"})
        one = sw.ddg_rows(model, domain, [a], 2, samples=2)[0]
        both = sw.ddg_rows(model, domain, [b, a], 2, samples=2, batch=1)
        assert both[1] == pytest.approx(one, abs=1e-4)                         # the same row with another neighbour in the batch
        assert sw.ddg_rows(model, domain, [a], 2, samples=2)[0] == one          # and again, the same bytes
        assert sw.ddg_rows(model, domain, [a], 3, samples=2)[0] != one
        assert sw.ddg_rows(model, domain, [cs.encode(domain["seq"], {})], 2, samples=2).tolist() == [0.0]

    def test_binding_is_the_complex_term_minus_the_two_partners(self):
        assert sw.binding(np.array([5.0, 1.0]), 2.0, np.array([1.0, 0.0])).tolist() == [2.0, -1.0]


class Counting:
    """Stands in for the model, in its own sign (negative is worse): every substituted place costs 2 in the system and 1 in a
    partner, so the binding prediction is the number of substitutions of both partners together."""
    device = "cpu"
    use_antithetic_variates = True

    def __init__(self):
        self.seen = []

    def _get_decoding_order(self, chain_M):
        import torch
        return torch.zeros_like(chain_M)

    def _get_backbone_noise(self, X):
        import torch
        return torch.zeros_like(X)

    def folding_dG(self, domain, seqs, decoding_order=None, backbone_noise=None):
        import torch
        wt = torch.from_numpy(cs.encode(domain["seq"], {}))
        n = (seqs != wt[None, :]).sum(dim=1).double()
        self.seen.append((domain["name"], len(seqs)))
        return -(1.0 if domain["name"].endswith(("_mt", "_nuc")) else 2.0) * n


def real_inputs(complex_="IV", gene="COX4I1"):
    pytest.importorskip("torch"); pytest.importorskip("stabddg")
    cx = COMPLEXES[complex_]
    gb, aln, sp = paths.CYBRID_GB / "records.gb", paths.MITOPORT / "partners" / "aln", paths.BUILD / "species.parquet"
    if not (cx.cif.is_file() and gb.is_file() and (aln / f"{gene}.afa").is_file() and sp.is_file() and (cx.work / "chains.tsv").is_file()):
        pytest.skip("needs the structures, the cybrid records, the vertebrate alignments and the build")
    return gb, aln, sp


class TestTask:
    def test_one_task_on_the_real_structure_with_a_counting_model(self, tmp_path):
        gb, aln, sp = real_inputs("IV", "COX5B")
        models = {}
        scores, case = sw.run("gorilla-gorilla", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, seeds=(0, 1), samples=2, draws=6,
                              model_for=lambda name: models.setdefault(name, Counting()))
        assert case["status"] == "scored" and case["copies"] == 1 and case["candidates"] >= sw.MIN_CANDIDATES and case["null"] == "assessed"
        assert case["placed"] + case["not_placed"] == 98 and case["substituted"] > 0 and case["contact_differences"] <= case["substituted"]
        for seed in (0, 1):
            s = scores[scores.seed == seed].set_index(["kind", "draw"])
            nmt = s.loc[("B0", -1), "fold_one"]
            assert nmt == case["mt_substituted"] > 0
            assert s.loc[("B0", -1), "ddG"] == nmt and s.loc[("B0", -1), "substituted"] == 0 and s.loc[("B0", -1), "fold_two"] == 0
            # each substitution of the subunit counts once more in the binding prediction, and twice in the system
            assert s.loc[("B1", -1), "ddG"] == nmt + case["substituted"] and s.loc[("B1", -1), "fold_two"] == case["substituted"]
            assert s.loc[("B1", -1), "fold"] == 2 * (nmt + case["substituted"])
            nulls = s.loc["null"]
            assert len(nulls) == 6 and (nulls.ddG == nmt + nulls.substituted).all() and (nulls.fold == 2 * (nmt + nulls.substituted)).all()
        again = pd.read_csv(tmp_path / "swap" / "gorilla-gorilla.IV.COX5B.tsv", sep="\t")
        assert len(again) == len(scores) == 2 * (2 + 6) and list(again.columns) == sw.SCORE_COLUMNS
        assert pd.read_csv(tmp_path / "swap" / "gorilla-gorilla.IV.COX5B.cases.tsv", sep="\t").status.iloc[0] == "scored"
        # the three structures a task reads: the system with the subunit, the mtDNA-encoded chains, the subunit
        assert {n for m in models.values() for n, _ in m.seen} == {"5Z62_COX5B", "5Z62_mt", "5Z62_COX5B_nuc"}

    def test_the_species_drawn_are_mammals_within_the_window_of_the_donors_differences(self, tmp_path):
        gb, aln, sp = real_inputs("IV", "COX5B")
        scores, case = sw.run("pan-paniscus", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, seeds=(0,), samples=1, draws=12,
                              model_for=lambda name: Counting())
        rows = sw.read_alignment(aln / "COX5B.afa")
        mammals = set(pd.read_parquet(sp, columns=["slug", "class"]).query("`class` == 'Mammalia'").slug)
        ct = pd.read_csv(COMPLEXES["IV"].work / "chains.tsv", sep="\t", dtype={"chain": str})
        chains = cs.read_chains(COMPLEXES["IV"].cif)
        by_id = {c["id"]: i for i, c in enumerate(chains)}
        h = rows["homo-sapiens"]
        u = by_id[ct.chain[ct.gene == "COX5B"].iloc[0]]
        columns = sw.contact_columns({u: sw.chain_columns(chains[u]["seq"], h)}, sw.contact_flags(chains, [by_id[c] for c in "ABC"], [u]))
        n = sum(1 for c in columns if rows["pan-paniscus"][c] in sw.THREE and h[c] in sw.THREE and rows["pan-paniscus"][c] != h[c])
        assert n == case["contact_differences"] > 0
        drawn = set(scores[scores.kind == "null"].species)
        assert drawn and drawn <= mammals - {"pan-paniscus", "homo-sapiens"}
        for sp_ in drawn:
            k = sum(1 for c in columns if rows[sp_][c] in sw.THREE and h[c] in sw.THREE and rows[sp_][c] != h[c])
            assert 0.75 * n <= k <= 1.25 * n
        by_species = scores[scores.kind == "null"].groupby("species").substituted.nunique()
        assert (by_species == 1).all()                                                  # a species is read as one set, however often it is drawn

    def test_every_copy_of_a_subunit_carries_the_ortholog(self, tmp_path):
        gb, aln, sp = real_inputs("III", "UQCRFS1")
        scores, case = sw.run("pan-troglodytes", "III", "UQCRFS1", gb, paths.PDB, aln, sp, tmp_path, tmp_path, seeds=(0,), samples=1, draws=3,
                              model_for=lambda name: Counting())
        rows = sw.read_alignment(aln / "UQCRFS1.afa")
        ct = pd.read_csv(COMPLEXES["III"].work / "chains.tsv", sep="\t", dtype={"chain": str})
        hp = cybrid.species_proteins(cybrid.read_proteins(gb), cybrid.HUMAN)
        genes = tuple(COMPLEXES["III"].mt)
        kept, mt_chain = cs.system(cs.read_chains(COMPLEXES["III"].cif), {g: hp[g][0] for g in genes}, {g: p for g, (p, _, _) in hp.items() if g not in genes})
        ids = set(ct.chain[ct.gene == "UQCRFS1"])
        copies = [c for i, c in enumerate(kept) if i not in mt_chain and c["id"] in ids]
        assert case["copies"] == len(copies) == 2                                 # of the four chains of the Rieske protein, the two that touch cytochrome b
        per_copy = [sw.ortholog_subs(c["seq"], sw.chain_columns(c["seq"], rows["homo-sapiens"]), rows["pan-troglodytes"])[0] for c in copies]
        assert case["substituted"] == sum(len(x) for x in per_copy) > 0
        b1 = scores[scores.kind == "B1"].iloc[0]
        assert b1.fold_two == case["substituted"] and b1.ddG == case["mt_substituted"] + case["substituted"]

    def test_an_ortholog_identical_to_human_changes_nothing_and_runs_no_null(self, tmp_path):
        gb, aln, sp = real_inputs("IV", "COX6B1")
        scores, case = sw.run("gorilla-gorilla", "IV", "COX6B1", gb, paths.PDB, aln, sp, tmp_path, tmp_path, seeds=(0,), samples=1, draws=5,
                              model_for=lambda name: Counting())
        assert case["status"] == "scored" and case["substituted"] == 0 and case["null"].startswith("not run")
        assert set(scores.kind) == {"B0", "B1"} and (scores[scores.kind == "B1"].ddG.values == scores[scores.kind == "B0"].ddG.values).all()
        res = sw.restoration_table(pd.DataFrame([case]), scores)
        assert res.restoration.iloc[0] == 0 and res.null_draws.iloc[0] == 0 and not res.in_top.iloc[0] and not res.exceeds_null.iloc[0]

    def test_a_gene_without_an_alignment_a_donor_or_a_chain_is_not_assessed_and_said(self, tmp_path):
        gb, aln, sp = real_inputs("IV", "COX5B")
        m = lambda name: Counting()
        _, c = sw.run("gorilla-gorilla", "IV", "COX6B2", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=m)
        assert c["status"] == "not assessed: no alignment of the gene"
        _, c = sw.run("pongo-abelii", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=m)
        assert c["status"] == "not assessed: no ortholog of the donor in the alignment"
        s, c = sw.run("gorilla-gorilla", "IV", "COX5B", gb, paths.PDB, tmp_path / "none", sp, tmp_path, tmp_path, model_for=m)
        assert c["status"] == "not assessed: no alignment of the gene" and s.empty
        assert pd.read_csv(tmp_path / "swap" / "gorilla-gorilla.IV.COX5B.tsv", sep="\t").empty
        if COMPLEXES["III"].cif.is_file() and (COMPLEXES["III"].work / "chains.tsv").is_file():
            # a gene with an alignment whose chain is not in this complex's structure
            _, c = sw.run("gorilla-gorilla", "III", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=m)
            assert c["status"] == "not assessed: no chain of the gene touches the mitochondrial chains"


FIXTURES = __import__("pathlib").Path(__file__).parent / "data" / "donor_orthologs"
REAL = do.parse_fasta((FIXTURES / "proteins.faa").read_text())          # RefSeq proteins of NDUFS2 and ATP5MC1 (tests/test_donor_orthologs.py)
MOUSE_NDUFS2, RAT_NDUFS2 = REAL["NP_694704.1"], REAL["NP_001011907.1"]


class TestProteinReading:
    """H36 (donor orthologs from RefSeq): the donor's protein is aligned to the structure's chain itself."""

    def test_the_residue_set_against_each_residue_of_the_chain_is_the_donors_through_a_presequence_and_a_tail(self):
        chain = MOUSE_NDUFS2[40:140]                                        # a part of the protein, as a structure holds it
        assert sw.place_protein(chain, MOUSE_NDUFS2) == chain                # the whole protein around it: end gaps are free
        assert sw.place_protein(chain, "MAAGGSTPLR" + MOUSE_NDUFS2 + "GGGG") == chain
        row = sw.place_protein(chain, RAT_NDUFS2)
        assert len(row) == len(chain) and row == RAT_NDUFS2[40:140]         # one length, no insertion: the rat's own residues
        assert row != chain                                                  # and the rat's NDUFS2 does differ from the mouse's there

    def test_the_substitutions_are_where_the_donor_differs_and_every_other_residue_is_placed(self):
        chain = MOUSE_NDUFS2[40:140]
        subs, placed, lost, row = sw.protein_subs(chain, RAT_NDUFS2)
        assert subs == {i: RAT_NDUFS2[40 + i] for i in range(100) if RAT_NDUFS2[40 + i] != chain[i]} and 0 < len(subs) < 20
        assert placed == 100 and lost == 0 and row == RAT_NDUFS2[40:140]
        assert sw.protein_subs(chain, MOUSE_NDUFS2)[:3] == ({}, 100, 0)       # the same protein changes nothing

    def test_a_residue_set_against_a_gap_or_an_unknown_or_in_a_chain_that_is_not_an_amino_acid_is_not_placed(self):
        chain = MOUSE_NDUFS2[40:140]
        gapped = chain[:50] + chain[52:]                                      # the donor lacks residues 50 and 51
        subs, placed, lost, row = sw.protein_subs(chain, gapped)
        assert row[50:52] == "--" and subs == {} and placed == 98 and lost == 2
        subs, placed, lost, _ = sw.protein_subs(chain, chain[:30] + "X" + chain[31:])
        assert subs == {} and placed == 99 and lost == 1                       # an unknown in the donor
        subs, placed, lost, _ = sw.protein_subs(chain[:30] + "X" + chain[31:], chain)
        assert subs == {} and placed == 99 and lost == 1                       # an unknown in the structure
        subs, placed, lost, _ = sw.protein_subs(chain, chain[:30] + "a" + chain[31:])
        assert subs == {30: "A"} if chain[30] != "A" else subs == {}           # the donor's letters are read in upper case
        assert placed + lost == 100

    def test_a_chain_that_is_the_start_or_the_end_of_the_protein_keeps_to_that_end(self):
        # with end gaps that cost, the first residues of the chain were set against a shifted stretch (MAALRAVQCL for the rat's MAALRAVCSL)
        assert sw.place_protein(MOUSE_NDUFS2[:10], RAT_NDUFS2) == RAT_NDUFS2[:10]
        assert sw.place_protein(MOUSE_NDUFS2[-10:], RAT_NDUFS2) == RAT_NDUFS2[-10:]
        assert sw.protein_subs(MOUSE_NDUFS2[:10], RAT_NDUFS2)[3] == RAT_NDUFS2[:10]

    def test_a_donor_that_holds_only_a_part_of_the_chain_leaves_the_rest_not_placed(self):
        chain = MOUSE_NDUFS2[40:140]
        subs, placed, lost, row = sw.protein_subs(chain, MOUSE_NDUFS2[60:120])
        assert row == "-" * 20 + MOUSE_NDUFS2[60:120] + "-" * 20 and subs == {} and placed == 60 and lost == 40

    def test_the_real_alignments_and_the_proteins_give_the_same_residues(self):
        f = paths.MITOPORT / "partners" / "aln" / "COX4I1.afa"
        if not (f.is_file() and (paths.PDB / "5Z62.cif").is_file()):
            pytest.skip("needs the vertebrate alignments of `mitoport iface-align` and 5Z62")
        rows = sw.read_alignment(f)
        ct = pd.read_csv(COMPLEXES["IV"].work / "chains.tsv", sep="\t", dtype={"chain": str})
        chain = [c for c in cs.read_chains(paths.PDB / "5Z62.cif") if c["id"] == ct.chain[ct.gene == "COX4I1"].iloc[0]][0]
        cols = sw.chain_columns(chain["seq"], rows["homo-sapiens"])
        for sp in ("mus-musculus", "gorilla-gorilla", "rattus-norvegicus"):
            by_alignment = sw.ortholog_subs(chain["seq"], cols, rows[sp])
            by_protein = sw.protein_subs(chain["seq"], rows[sp].replace("-", ""))
            assert by_protein[0] == by_alignment[0] and by_protein[1] == by_alignment[1] and by_protein[2] == by_alignment[2]
            assert by_alignment[0]                                            # each differs from human somewhere

    def test_contact_differences_count_a_position_once_across_copies_and_not_a_gap_an_unknown_or_an_unplaced_residue(self):
        seqs = {0: "KVLAWYRTE", 1: "KVLAWYRTE"}
        rows = {0: "KILAWYRTE", 1: "KILGWY-TX"}                                # copy 0: V2->I; copy 1: V2->I, A4->G, a gap and an unknown
        ids = {0: list(range(9)), 1: list(range(9))}
        flags = {0: np.array([1, 1, 0, 0, 0, 0, 0, 0, 0], bool), 1: np.array([1, 1, 0, 1, 0, 0, 1, 0, 1], bool)}
        assert sw.protein_contact_differences(seqs, rows, ids, flags) == 2     # position 1 (in both copies) and 3; the gap (6) and the unknown (8) are not
        assert sw.protein_contact_differences(seqs, rows, ids, {0: flags[0] & False, 1: flags[1] & False}) == 0
        ids[1] = [None] * 9
        assert sw.protein_contact_differences(seqs, rows, ids, flags) == 1     # a residue without an identity (no column) is not counted
        assert sw.protein_contact_differences({0: "KVX"}, {0: "KIX"}, {0: [0, 1, 2]}, {0: np.array([1, 1, 1], bool)}) == 1
        assert sw.protein_contact_differences({0: "KVX"}, {0: "KIA"}, {0: [0, 1, 2]}, {0: np.array([1, 1, 1], bool)}) == 1      # an unknown of the chain is no difference

    def test_the_donors_number_of_differences_can_be_given_for_a_donor_the_alignment_does_not_hold(self):
        rows = {"homo-sapiens": "ACDEFGHIKL", "a": "ACDEFGHIKL", "b": "VCDEFGHIKL", "c": "VCDEFGHIKM", "d": "VCDEFGHIKM", "e": "VCDEFGHIMM"}
        mammals = set(rows)
        columns = list(range(10))
        cand, n = sw.null_candidates(rows, mammals, "not-in-the-alignment", columns, n=1)
        assert n == 1 and cand == ["b"]                                       # 0.75 <= k <= 1.25 of 1: exactly one difference
        cand, n = sw.null_candidates(rows, mammals, "not-in-the-alignment", columns, n=2)
        assert n == 2 and cand == ["c", "d"]                                  # 1.5 <= k <= 2.5
        cand, n = sw.null_candidates(rows, mammals, "c", columns)             # without it, the donor's row gives it
        assert n == 2 and cand == ["d"]


class TestProteinTask:
    def test_a_protein_that_is_the_alignments_row_gives_the_alignment_readings_whole_task(self, tmp_path):
        gb, aln, sp = real_inputs("IV", "COX5B")
        rows = sw.read_alignment(aln / "COX5B.afa")
        m = lambda name: Counting()
        by_alignment, case_a = sw.run("gorilla-gorilla", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path / "a", tmp_path, seeds=(0,), samples=1, draws=6, model_for=m)
        source = do.Source({("gorilla-gorilla", "COX5B"): do.Ortholog(rows["gorilla-gorilla"].replace("-", ""), "XP_0.1", "1")}, {}, {"gorilla-gorilla": "x"})
        by_protein, case_p = sw.run("gorilla-gorilla", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path / "p", tmp_path, seeds=(0,), samples=1, draws=6,
                                    model_for=m, donor_source=source)
        for k in ("status", "copies", "contact_residues", "placed", "not_placed", "substituted", "contact_differences", "candidates", "null", "mt_substituted"):
            assert case_p[k] == case_a[k], k
        assert case_p["reading"] == "protein" and case_p["donor_protein"] == "XP_0.1" and "reading" not in case_a        # the alignment reading's table is as it was
        pd.testing.assert_frame_equal(by_protein, by_alignment)

    def test_every_copy_of_a_subunit_carries_the_donors_protein_residues(self, tmp_path):
        gb, aln, sp = real_inputs("III", "UQCRFS1")
        rows = sw.read_alignment(aln / "UQCRFS1.afa")
        m = lambda name: Counting()
        _, by_alignment = sw.run("pan-troglodytes", "III", "UQCRFS1", gb, paths.PDB, aln, sp, tmp_path / "a", tmp_path, seeds=(0,), samples=1, draws=3, model_for=m)
        source = do.Source({("pan-troglodytes", "UQCRFS1"): do.Ortholog(rows["pan-troglodytes"].replace("-", ""), "XP_0.1", "1")}, {}, {"pan-troglodytes": "x"})
        scores, case = sw.run("pan-troglodytes", "III", "UQCRFS1", gb, paths.PDB, aln, sp, tmp_path / "p", tmp_path, seeds=(0,), samples=1, draws=3,
                              model_for=m, donor_source=source)
        assert case["copies"] == by_alignment["copies"] == 2 and case["substituted"] == by_alignment["substituted"] > 1
        b1 = scores[scores.kind == "B1"].iloc[0]
        assert b1.fold_two == case["substituted"] and b1.ddG == case["mt_substituted"] + case["substituted"]

    def test_a_gene_without_an_alignment_is_scored_with_the_donors_protein_and_has_no_null(self, tmp_path):
        d = paths.DONOR_ORTHOLOGS
        if not (d / "mus-pahari.faa").is_file():
            pytest.skip("needs `mitoage download donor-orthologs`")
        gb, aln, sp = real_inputs("IV", "COX5B")
        source = do.read_source(d)
        # every gene a task can use has an alignment (they are the genes that touch the mtDNA-encoded chains); a folder without one stands for a gene without
        scores, case = sw.run("mus-pahari", "IV", "COX5B", gb, paths.PDB, tmp_path / "none", sp, tmp_path, tmp_path, seeds=(0,), samples=1, draws=5,
                              model_for=lambda name: Counting(), donor_source=source)
        assert case["status"] == "scored" and case["reading"] == "protein" and case["donor_protein"] == source.orthologs[("mus-pahari", "COX5B")].accession
        assert case["copies"] == 1 and case["substituted"] > 0 and case["contact_residues"] > 0 and case["candidates"] == 0 and case["contact_differences"] > 0
        assert case["null"].startswith("not assessed: no alignment of the gene") and set(scores.kind) == {"B0", "B1"}
        b1 = scores[scores.kind == "B1"].iloc[0]
        assert b1.fold_two == case["substituted"] and b1.ddG == case["mt_substituted"] + case["substituted"]
        # the same gene with its alignment: the contact residues are those of the alignment's columns and the null is assessed
        _, with_aln = sw.run("mus-pahari", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path / "w", tmp_path, seeds=(0,), samples=1, draws=5,
                             model_for=lambda name: Counting(), donor_source=source)
        assert with_aln["substituted"] == case["substituted"] and with_aln["contact_differences"] == case["contact_differences"]
        assert with_aln["contact_residues"] == case["contact_residues"] and with_aln["candidates"] >= sw.MIN_CANDIDATES and with_aln["null"] == "assessed"

    def test_a_donor_the_table_covers_is_never_read_from_the_alignment_for_a_gene_it_lacks(self, tmp_path):
        gb, aln, sp = real_inputs("IV", "COX5B")
        m = lambda name: Counting()
        source = do.Source({}, {("gorilla-gorilla", "COX5B"): "no ortholog"}, {"gorilla-gorilla": "GCF_x", "mus-spretus": ""})
        s, c = sw.run("gorilla-gorilla", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=m, donor_source=source)
        assert c["status"] == "not assessed: no ortholog of the donor in RefSeq (no ortholog)" and s.empty       # the alignment holds a gorilla row, which is not used
        s, c = sw.run("mus-spretus", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=m, donor_source=source)
        assert c["status"].startswith("not assessed: no ortholog source") and s.empty
        _, c = sw.run("gorilla-gorilla", "IV", "COX6B2", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=m, donor_source=do.Source({}, {}, {"mus-musculus": "x"}))
        assert c["status"] == "not assessed: no alignment of the gene"           # a donor the table does not cover keeps the alignment reading


# --- the mouse host (H36, 2026-11-24) --------------------------------------------------------------------------------

MOUSE = sw.MOUSE


class TestMouseTasks:
    def test_the_donors_with_a_refseq_source_and_an_outcome_have_a_task_for_each_gene_and_no_complex_v(self):
        t = sw.tasks(host=MOUSE)
        assert {d for d, _, _ in t} == {"rattus-norvegicus", "mus-caroli", "mus-pahari"} and {c for _, c, _ in t} == {"I", "III", "IV"}
        for d in {"rattus-norvegicus", "mus-caroli", "mus-pahari"}:
            for c in ("I", "III", "IV"):
                assert [g for dd, cx, g in t if dd == d and cx == c] == COMPLEXES[c].nuclear
        assert len(t) == 3 * sum(len(COMPLEXES[c].nuclear) for c in ("I", "III", "IV"))
        assert t == sorted(t, key=lambda x: (x[1], x[0], COMPLEXES[x[1]].nuclear.index(x[2])))

    def test_a_donor_without_a_source_has_no_task_and_is_listed_as_not_assessed(self):
        t = sw.tasks(host=MOUSE)
        for d in ("mus-spretus", "mus-terricolor", "otomys-irroratus", "cricetulus-griseus", "mesocricetus-auratus"):
            assert d not in {x for x, _, _ in t}
        u = sw.unassessed_donors(host=MOUSE)
        assert set(zip(u.donor, u["complex"])) == {(d, c) for d in ("mus-spretus", "mus-terricolor", "otomys-irroratus") for c in ("I", "III", "IV")}
        assert (u.host == MOUSE).all() and u.status.str.startswith("not assessed: no ortholog source").all()      # the hamsters have only whole-cybrid rows

    def test_a_donor_the_species_table_does_not_know_is_refused(self, tmp_path):
        m = cybrid.read_measurements().copy()
        m = pd.concat([m, m[(m.host == MOUSE) & (m.donor == "mus-pahari")].assign(donor="mus-unknown")])
        f = tmp_path / "m.tsv"
        m.drop(columns="pct").to_csv(f, sep="\t", index=False)
        with pytest.raises(SystemExit):
            sw.tasks(f, MOUSE)

    def test_the_human_host_keeps_its_donors(self):
        assert {d for d, _, _ in sw.tasks()} == {"pan-troglodytes", "pan-paniscus", "gorilla-gorilla"} and len(sw.tasks()) == 216


class TestChainsBySequence:
    A = "MKVLAWYRTEQNDSGHIFPCMKVLAWYRTEQNDSGHIFPC"
    B = "GSTNDEQKRHGSTNDEQKRHAAAAKKKKLLLLMMMMVVVV"

    def test_a_chain_is_the_gene_whose_protein_it_copies_and_every_copy_is_listed(self):
        chains = [{"id": "m", "seq": "X" * 30}, {"id": "a1", "seq": self.A}, {"id": "b", "seq": self.B}, {"id": "a2", "seq": self.A}]
        got = sw.gene_chains_by_sequence(chains, {0: ("ND1", {})}, {"GA": self.A, "GB": self.B})
        assert got == {"GA": [1, 3], "GB": [2]}

    def test_the_mitochondrial_chains_another_isoform_and_a_short_chain_are_no_gene(self):
        other = self.A[:20] + "W" * 20                                         # half of the residues of GA changed: another isoform
        chains = [{"id": "m", "seq": self.A}, {"id": "iso", "seq": other}, {"id": "short", "seq": self.A[:15]}, {"id": "none", "seq": "P" * 40}]
        assert sw.gene_chains_by_sequence(chains, {0: ("ND1", {})}, {"GA": self.A}) == {}
        assert sw.gene_chains_by_sequence([], {}, {"GA": self.A}) == {}
        near = self.A[:37] + "WWW"                                             # 3 of 40 changed (92.5% identical): still the gene's
        assert sw.gene_chains_by_sequence([{"id": "n", "seq": near}], {}, {"GA": self.A}) == {"GA": [0]}

    def test_the_real_mouse_structures_name_their_subunits_and_leave_the_other_isoforms_out(self):
        pdb, d = paths.PDB / "7O37.cif", paths.DONOR_ORTHOLOGS
        if not (pdb.is_file() and (d / "mus-musculus.faa").is_file() and (paths.CYBRID_GB / "records.gb").is_file()):
            pytest.skip("needs 7O37, `mitoage download donor-orthologs` and the cybrid records")
        source = do.read_source(d)
        hp = cybrid.species_proteins(cybrid.read_proteins(paths.CYBRID_GB / "records.gb"), MOUSE)
        genes = tuple(COMPLEXES["IV"].mt)
        chains, mt_chain = cs.system(cs.read_chains(pdb), {g: hp[g][0] for g in genes}, {g: p for g, (p, _, _) in hp.items() if g not in genes})
        nuclear = {g: source.get(MOUSE, g)[0].protein for g in COMPLEXES["IV"].nuclear}
        by = sw.gene_chains_by_sequence(chains, mt_chain, nuclear)
        assert set(by) == {"COX4I1", "COX5A", "COX5B", "COX6B1", "COX6C", "COX7B", "COX7C"}          # not COX6A1, COX7A2, COX8A: the structure holds other isoforms
        assert all(len(v) == 1 and v[0] not in mt_chain for v in by.values())
        for g, (i,) in by.items():                                                # the chains are the mouse's RefSeq proteins, residue for residue
            assert sw.protein_subs(chains[i]["seq"], nuclear[g])[:3] == ({}, len(chains[i]["seq"]), 0)


def mouse_inputs(complex_="IV", gene="COX5B"):
    pytest.importorskip("torch"); pytest.importorskip("stabddg")
    pid = sw.STRUCTURES[MOUSE][complex_]
    gb, aln, sp, d = paths.CYBRID_GB / "records.gb", paths.MITOPORT / "partners" / "aln", paths.BUILD / "species.parquet", paths.DONOR_ORTHOLOGS
    if not ((paths.PDB / f"{pid}.cif").is_file() and gb.is_file() and (aln / f"{gene}.afa").is_file() and sp.is_file() and (d / "mus-pahari.faa").is_file()):
        pytest.skip("needs the mouse structures, the cybrid records, the vertebrate alignments, the build and `mitoage download donor-orthologs`")
    return gb, aln, sp, do.read_source(d)


class TestMouseTask:
    def run(self, tmp_path, donor="rattus-norvegicus", complex_="IV", gene="COX5B", source=None, draws=8, **kw):
        gb, aln, sp, src = mouse_inputs(complex_, gene)
        models = {}
        out = sw.run(donor, complex_, gene, gb, paths.PDB, aln, sp, tmp_path, tmp_path, seeds=(0,), samples=1, draws=draws,
                     model_for=lambda name: models.setdefault(name, Counting()), donor_source=source or src, host=MOUSE, **kw)
        return (*out, models, (gb, aln, sp, src))

    def test_the_mouse_complex_is_the_system_its_own_subunit_is_b0_and_the_donors_is_b1(self, tmp_path):
        scores, case, models, (gb, aln, sp, src) = self.run(tmp_path)
        assert case["host"] == MOUSE and case["status"] == "scored" and case["pdb"] == "7O37" and case["host_differs"] == 0
        assert case["host_protein"] == src.get(MOUSE, "COX5B")[0].accession and case["donor_protein"] == src.get("rattus-norvegicus", "COX5B")[0].accession
        assert {n for m in models.values() for n, _ in m.seen} == {"7O37_COX5B", "7O37_mt", "7O37_COX5B_nuc"}       # never the human structure
        assert (scores.host == MOUSE).all() and scores[scores.kind == "B0"].species.iloc[0] == MOUSE and scores[scores.kind == "B1"].species.iloc[0] == "rattus-norvegicus"
        chain = [c for c in cs.read_chains(paths.PDB / "7O37.cif") if c["id"] == "f"][0]
        subs = sw.protein_subs(chain["seq"], src.get("rattus-norvegicus", "COX5B")[0].protein)[0]
        assert case["substituted"] == len(subs) > 0
        b0, b1 = scores[scores.kind == "B0"].iloc[0], scores[scores.kind == "B1"].iloc[0]
        assert b0.substituted == 0 and b0.ddG == case["mt_substituted"] > 0           # B0: the mouse's subunit, the rat's mitochondrial changes
        assert b1.fold_two == case["substituted"] and b1.ddG == case["mt_substituted"] + case["substituted"]
        assert pd.read_csv(tmp_path / "swap" / "rattus-norvegicus.IV.COX5B.cases.tsv", sep="\t").host.iloc[0] == MOUSE

    def test_the_null_draws_mammals_that_differ_from_the_mouse_as_often_as_the_donor_does(self, tmp_path):
        scores, case, _, (gb, aln, sp, src) = self.run(tmp_path, draws=12)
        rows = sw.read_alignment(aln / "COX5B.afa")
        mammals = set(pd.read_parquet(sp, columns=["slug", "class"]).query("`class` == 'Mammalia'").slug)
        hp = cybrid.species_proteins(cybrid.read_proteins(gb), MOUSE)
        genes = tuple(COMPLEXES["IV"].mt)
        chains, mt_chain = cs.system(cs.read_chains(paths.PDB / "7O37.cif"), {g: hp[g][0] for g in genes}, {g: p for g, (p, _, _) in hp.items() if g not in genes})
        u = [i for i, c in enumerate(chains) if c["id"] == "f"][0]
        mt = sorted(mt_chain)
        m = rows[MOUSE]
        cols = sw.chain_columns(chains[u]["seq"], m)
        columns = sw.contact_columns({u: cols}, sw.contact_flags(chains, mt, [u]))
        assert case["contact_residues"] == len(columns) > 0
        drawn = set(scores[scores.kind == "null"].species)
        assert drawn and drawn <= mammals - {MOUSE, "rattus-norvegicus"}
        n = case["contact_differences"]
        for s in drawn:
            k = sum(1 for c in columns if rows[s][c] in sw.THREE and m[c] in sw.THREE and rows[s][c] != m[c])
            assert 0.75 * n <= k <= 1.25 * n
        assert case["candidates"] == len(sw.null_candidates(rows, mammals, "rattus-norvegicus", columns, n=n, host=MOUSE)[0]) >= sw.MIN_CANDIDATES
        by_human = len(sw.null_candidates(rows, mammals, "rattus-norvegicus", columns, n=n)[0])           # measured from the human row it is another set
        assert by_human != case["candidates"] or by_human == 0

    def test_a_mouse_host_needs_the_pinned_orthologs(self, tmp_path):
        gb, aln, sp, _ = mouse_inputs()
        with pytest.raises(SystemExit):
            sw.run("rattus-norvegicus", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=lambda n: Counting(), host=MOUSE)

    def test_a_complex_without_a_mouse_structure_is_not_assessed_and_the_human_one_is_not_used(self, tmp_path):
        gb, aln, sp, src = mouse_inputs()
        m = {}
        s, c = sw.run("rattus-norvegicus", "V", "ATP5PF", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=lambda n: m.setdefault(n, Counting()), donor_source=src, host=MOUSE)
        assert c["status"].startswith("not assessed: no solved structure of the mus-musculus complex V") and s.empty and not m

    def test_a_gene_without_a_chain_of_the_host_in_the_structure_is_not_assessed_and_says_why(self, tmp_path):
        s, c, models, _ = self.run(tmp_path, gene="COX7A2")          # the structure holds another isoform (COX7A2L) there
        assert c["status"].startswith("not assessed: no chain of 7O37") and "COX7A2" in c["status"] and s.empty and not models

    def test_a_missing_ortholog_of_the_host_or_the_donor_is_not_assessed_and_a_donor_the_table_lacks_is_not_read_from_the_alignment(self, tmp_path):
        gb, aln, sp, src = mouse_inputs()
        m = lambda n: Counting()
        no_host = do.Source({k: v for k, v in src.orthologs.items() if k != (MOUSE, "COX5B")}, {**src.status, (MOUSE, "COX5B"): "no ortholog"}, src.species)
        s, c = sw.run("rattus-norvegicus", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=m, donor_source=no_host, host=MOUSE)
        assert c["status"] == "not assessed: no ortholog of the host in RefSeq (no ortholog)" and s.empty and c["host"] == MOUSE
        s, c = sw.run("rattus-norvegicus", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=m,
                      donor_source=do.Source({}, {}, {"rattus-norvegicus": "x"}), host=MOUSE)
        assert c["status"] == "not assessed: no pinned source for the host"
        s, c = sw.run("mus-spretus", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=m, donor_source=src, host=MOUSE)
        assert c["status"].startswith("not assessed: no ortholog source") and s.empty
        gone = do.Source(src.orthologs, src.status, {k: v for k, v in src.species.items() if k != "rattus-norvegicus"})
        s, c = sw.run("rattus-norvegicus", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=m, donor_source=gone, host=MOUSE)
        assert c["status"] == "not assessed: the pinned table does not cover the donor" and s.empty         # the alignment holds a rat row; it is not used
        no_donor = do.Source(src.orthologs, {**src.status, ("mus-caroli", "COX5B"): "no ortholog"}, src.species)
        no_donor.orthologs = {k: v for k, v in src.orthologs.items() if k != ("mus-caroli", "COX5B")}
        _, c = sw.run("mus-caroli", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=m, donor_source=no_donor, host=MOUSE)
        assert c["status"] == "not assessed: no ortholog of the donor in RefSeq (no ortholog)"

    def test_a_chain_that_is_not_the_hosts_protein_is_not_assessed_because_b0_would_not_be_the_hosts_own(self, tmp_path):
        gb, aln, sp, src = mouse_inputs()
        host = src.get(MOUSE, "COX5B")[0]
        p = host.protein
        changed = p[:40] + ("W" if p[40] != "W" else "F") + p[41:60] + ("W" if p[60] != "W" else "F") + p[61:]
        other = do.Source({**src.orthologs, (MOUSE, "COX5B"): do.Ortholog(changed, host.accession, host.gene_id)}, src.status, src.species)
        s, c = sw.run("rattus-norvegicus", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path, tmp_path, model_for=lambda n: Counting(), donor_source=other, host=MOUSE)
        assert c["status"] == "not assessed: the chains of 7O37 differ from the host's RefSeq protein at 2 placed residues, so the chain as solved is not the host's own subunit"
        assert s.empty

    def test_a_donor_identical_to_the_host_or_with_a_gap_changes_nothing_and_runs_no_null(self, tmp_path):
        gb, aln, sp, src = mouse_inputs()
        host = src.get(MOUSE, "COX5B")[0]
        for name, protein, lost in (("same", host.protein, 0), ("gap", host.protein[:80] + host.protein[83:], 3)):
            other = do.Source({**src.orthologs, ("rattus-norvegicus", "COX5B"): do.Ortholog(protein, "XP_0.1", "1")}, src.status, src.species)
            s, c = sw.run("rattus-norvegicus", "IV", "COX5B", gb, paths.PDB, aln, sp, tmp_path / name, tmp_path, seeds=(0,), samples=1, draws=5,
                          model_for=lambda n: Counting(), donor_source=other, host=MOUSE)
            assert c["status"] == "scored" and c["substituted"] == 0 and c["null"].startswith("not run") and c["not_placed"] == lost
            assert set(s.kind) == {"B0", "B1"} and (s[s.kind == "B1"].ddG.values == s[s.kind == "B0"].ddG.values).all()
            r = sw.restoration_table(pd.DataFrame([c]), s)
            assert r.restoration.iloc[0] == 0 and not r.in_top.iloc[0] and not r.exceeds_null.iloc[0]

    def test_every_copy_of_the_complex_iii_subunit_carries_the_donors_protein_in_the_mouse_dimer(self, tmp_path):
        gb, aln, sp, src = mouse_inputs("III", "UQCRFS1")
        s, c = sw.run("mus-pahari", "III", "UQCRFS1", gb, paths.PDB, aln, sp, tmp_path, tmp_path, seeds=(0,), samples=1, draws=3,
                      model_for=lambda n: Counting(), donor_source=src, host=MOUSE)
        chains = cs.read_chains(paths.PDB / "7O3H.cif")
        copies = [c_ for c_ in chains if c_["id"] in ("E", "P")]
        prot = src.get("mus-pahari", "UQCRFS1")[0].protein
        assert c["status"] == "scored" and c["copies"] == 2 and c["pdb"] == "7O3H"
        assert c["substituted"] == sum(len(sw.protein_subs(x["seq"], prot)[0]) for x in copies)
        assert s[s.kind == "B1"].iloc[0].fold_two == c["substituted"] > 0


class TestMouseCollect:
    def test_the_mouse_tasks_are_collected_with_the_donors_that_have_no_task_listed(self, tmp_path, capsys):
        (tmp_path / "swap").mkdir()
        d, c, g = sw.tasks(host=MOUSE)[0]
        case = {"host": MOUSE, "donor": d, "complex": c, "gene": g, "status": "scored", "substituted": 1, "contact_differences": 1, "candidates": 3}
        pd.DataFrame([case]).to_csv(tmp_path / "swap" / f"{d}.{c}.{g}.cases.tsv", sep="\t", index=False)
        rows = [{"host": MOUSE, "donor": d, "complex": c, "gene": g, "kind": k, "draw": dr, "species": "s", "seed": s, "substituted": 1, "fold": 0.0,
                 "fold_one": 0.0, "fold_two": 0.0, "ddG": v} for s in (0, 1, 2) for k, dr, v in (("B0", -1, 4.0), ("B1", -1, 1.0), ("null", 0, 3.0), ("null", 1, 2.0), ("null", 2, 3.5))]
        pd.DataFrame(rows).to_csv(tmp_path / "swap" / f"{d}.{c}.{g}.tsv", sep="\t", index=False)
        res = sw.collect(tmp_path, host=MOUSE)
        n = len(sw.tasks(host=MOUSE))
        assert f"MISSING: {n - 1} of {n}" in capsys.readouterr().out and len(res) == 1 and res.iloc[0].restoration == 3.0
        u = pd.read_csv(tmp_path / "swap_unassessed.tsv", sep="\t")
        assert len(u) == 9 and set(u.donor) == {"mus-spretus", "mus-terricolor", "otomys-irroratus"}

    def test_the_activity_of_the_mouse_host_is_the_mouse_rows(self):
        cases = pd.DataFrame([{"host": MOUSE, "donor": d, "complex": "I", "gene": "A", "status": "scored", "substituted": 1, "contact_differences": 1, "candidates": 5}
                              for d in ("d1", "d2", "d3")])
        scores = pd.DataFrame([{"donor": d, "complex": "I", "gene": "A", "kind": k, "draw": -1, "seed": 0, "ddG": v}
                               for d, b1 in (("d1", 1.0), ("d2", 2.0), ("d3", 3.0)) for k, v in (("B0", 4.0), ("B1", b1))])
        res = sw.restoration_table(cases, scores)
        base = {"assay": "x", "outcome": "measured", "significant": "yes", "source_doi": "x", "location": "x", "note": "", "complex": "I"}
        meas = pd.DataFrame([{**base, "host": MOUSE, "donor": d, "pct_low": p, "pct_high": p, "pct": p} for d, p in (("d1", 30.0), ("d2", 50.0), ("d3", 70.0))]
                            + [{**base, "host": "homo-sapiens", "donor": d, "pct_low": p, "pct_high": p, "pct": p} for d, p in (("d1", 90.0), ("d2", 50.0), ("d3", 10.0))])
        sums, rho = sw.activity_table(res, cases, meas, MOUSE)
        assert list(sums.pct) == [30.0, 50.0, 70.0] and rho.set_index("complex").loc["I", "spearman"] == pytest.approx(-1.0)
        assert list(sw.activity_table(res, cases, meas)[0].pct) == [90.0, 50.0, 10.0]          # the default host reads the human rows


class TestReducedReading:
    base = {"assay": "x", "outcome": "measured", "source_doi": "x", "location": "x", "note": "", "complex": "IV", "host": MOUSE}

    def frames(self, restorations):
        """restorations: {donor: {gene: (restoration, null draws' restorations)}}, B0 = 10."""
        cases, scores = frame(list(restorations), "IV", "A", {(d, g): {"B0": [10.0], "B1": [10.0 - r], "null": [[10.0 - n] for n in nulls]}
                                                             for d, genes in restorations.items() for g, (r, nulls) in genes.items()})
        cases["host"] = MOUSE
        return cases, sw.restoration_table(cases, scores)

    def meas(self, **status):
        rows = []
        for d, (sig, pct) in status.items():
            rows.append({**self.base, "donor": d, "significant": sig, "pct_low": pct, "pct_high": pct, "pct": pct})
        rows.append({**self.base, "donor": "dis", "significant": "yes", "pct_low": 50.0, "pct_high": 50.0, "pct": 50.0})
        rows.append({**self.base, "donor": "dis", "significant": "no", "pct_low": 95.0, "pct_high": 95.0, "pct": 95.0})
        rows.append({**self.base, "donor": "human-host-row", "host": "homo-sapiens", "significant": "yes", "pct_low": 1.0, "pct_high": 1.0, "pct": 1.0})
        return pd.DataFrame(rows)

    def test_the_rule_is_read_on_the_reduced_pairs_alone_and_the_not_reduced_are_counted_apart(self):
        low = [0.0] * 4
        cases, res = self.frames({"r1": {"A": (5.0, low), "B": (0.0, low)}, "r2": {"A": (5.0, low), "B": (5.0, low)},
                                  "n1": {"A": (5.0, low), "B": (0.0, low)}, "dis": {"A": (5.0, low), "B": (5.0, low)}})
        t = sw.reduced_table(res, cases, self.meas(r1=("yes", 50.0), r2=("yes", 40.0), n1=("no", 90.0)), MOUSE).set_index("gene")
        a, b = t.loc["A"], t.loc["B"]
        assert (a.reduced_donors, a.reduced_both, bool(a.candidate_reduced), a.not_reduced_donors, a.not_reduced_both) == (2, 2, True, 1, 1)
        assert (b.reduced_donors, b.reduced_both, bool(b.candidate_reduced), b.not_reduced_donors, b.not_reduced_both) == (2, 1, False, 1, 0)   # one of two is not more than half
        assert a.best_reduced_restoration == 5.0                               # the donor the papers disagree on ("dis") is in neither count

    def test_no_reduced_pair_is_never_a_candidate_and_another_hosts_rows_are_not_read(self):
        cases, res = self.frames({"n1": {"A": (5.0, [0.0] * 4)}, "h": {"A": (5.0, [0.0] * 4)}})
        t = sw.reduced_table(res, cases, self.meas(n1=("no", 90.0)), MOUSE).iloc[0]
        assert t.reduced_donors == 0 and not t.candidate_reduced and np.isnan(t.best_reduced_restoration) and t.not_reduced_donors == 1
