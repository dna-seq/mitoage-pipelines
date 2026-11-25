"""H36, step 3: SaProt on the joined mitochondrial and nuclear chains: the structure tokens and their match to the chains, the input
(the linker, the substitutions, the shuffled partner), the pseudo-log-likelihood of a masked place, the pairs of chains, and the tables
and the bar; on frames small enough to count by hand and, where the inputs are installed, on the real structures and alignments, with a
stand-in for the model whose output depends on the whole input (so that a change anywhere in the input is seen) and the real vocabulary
when it is installed. No test needs a GPU or the weights."""
from itertools import product

import numpy as np
import pandas as pd
import pytest
from scipy.special import logsumexp

from mitoage_pipelines import paths
from mitoage_pipelines.mitoport import cybrid, cybrid_saprot as sp, cybrid_stab as cs, cybrid_swap as sw, donor_orthologs as do
from mitoage_pipelines.mitoport.complexes import COMPLEXES

STATES = "pynwrqhgdlvtmfsaeikc"


def synthetic_vocab() -> dict[str, int]:
    toks = ["<cls>", "<pad>", "<eos>", "<unk>", "<mask>"] + [a + s for a in sp.AA + "#" for s in STATES + "#"]
    return {t: i for i, t in enumerate(toks)}


def real_vocab() -> dict[str, int]:
    f = paths.SAPROT / "vocab.txt"
    return sp.read_vocab(f) if f.is_file() else synthetic_vocab()


VOCAB = real_vocab()
V = len(VOCAB)


class Fake:
    """A stand-in for the model: its logits at a column are a function of the whole input and the column, so that a change of any residue
    or state anywhere moves them; it counts its calls and the rows it reads."""

    def __init__(self):
        self.calls, self.rows, self.seen = 0, 0, []

    def __call__(self, x, cols):
        self.calls += 1; self.rows += len(cols)
        out = np.empty((len(cols), V))
        for j in range(len(cols)):
            self.seen.append((x[j].copy(), int(cols[j])))
            seed = int((x[j].astype(np.int64) * (np.arange(x.shape[1]) % 97 + 1)).sum() % (2 ** 31)) + int(cols[j])
            out[j] = np.random.default_rng(seed).normal(size=V)
        return out


def chain(cid, seq, origin, num=None):
    n = len(seq)
    ca = np.array(origin, dtype=float) + np.c_[3.8 * np.arange(n), np.zeros(n), np.zeros(n)]
    bb = np.stack([ca + [-1, 0, 0], ca, ca + [1, 0, 0], ca + [1, 1, 0]], axis=1)
    return {"id": cid, "seq": seq, "num": list(num or range(1, n + 1)), "bb": bb, "heavy": bb.reshape(-1, 3), "n_heavy": [4] * n}


# --- the structure tokens ---------------------------------------------------------------------------------------------

LINES = ("5Z62.cif_A Structure of human cytochrome c oxidase\tMKVL\tDVCC\n5Z62.cif_B Structure of human cytochrome c oxidase\tAWYR\tQQVV\n")


class TestStructureTokens:
    def test_a_chain_is_named_by_the_text_after_the_last_underscore_and_the_states_are_lower_case(self):
        t = sp.parse_foldseek(LINES)
        assert t == {"A": ("MKVL", "dvcc"), "B": ("AWYR", "qqvv")}
        assert sp.parse_foldseek("x.cif_AA\tM\td\n") == {"AA": ("M", "d")} and sp.parse_foldseek("") == {}

    def test_a_malformed_or_repeated_chain_is_refused(self):
        with pytest.raises(SystemExit):
            sp.parse_foldseek("5Z62.cif_A title\tMKVL\n")                       # two fields
        with pytest.raises(SystemExit):
            sp.parse_foldseek("5Z62.cif_A t\tMKVL\tDVC\n")                      # a state missing
        with pytest.raises(SystemExit):
            sp.parse_foldseek(LINES + "5Z62.cif_A t\tMKVL\tDVCC\n")             # chain A twice

    def test_the_states_of_a_chain_are_taken_only_when_foldseeks_chain_is_the_chain(self):
        assert sp.chain_states("MKVL", "MKVL", "dvcc") == "dvcc"
        assert sp.chain_states("MKVL", "MXVL", "dvcc") == "dvcc"                # a modified residue Foldseek does not know, the structure reads as its parent
        assert sp.chain_states("MKVL", "MRVL", "dvcc") is None                  # another residue
        assert sp.chain_states("MKVL", "MKV", "dvc") is None                    # another length
        assert sp.chain_states("MKVL", "MKVL", "dvc") is None
        assert sp.chain_states("MXVL", "MKVL", "dvcc") is None                  # X in the chain and a residue in Foldseek's: not the same chain

    def test_the_real_structures_tokens_match_their_chains_but_a_few_modified_ones(self):
        cifs = {p: paths.PDB / f"{p}.cif" for p in ("5Z62", "5XTE", "8OM1")}
        if not all(f.is_file() for f in cifs.values()) or not all((paths.SAPROT_TOKENS / f"{p}.3di.tsv").is_file() for p in cifs):
            pytest.skip("needs the structures and `mitoage mitoport saprot-tokens`")
        for p, cif in cifs.items():
            chains = cs.read_chains(cif)
            tok = sp.structure_tokens(cif, None, paths.SAPROT_TOKENS)
            assert set(tok) == {c["id"] for c in chains}
            got = sp.chain_state_table(chains, tok)
            bad = [c["id"] for i, c in enumerate(chains) if got[i] is None]
            assert all(len(got[i]) == len(c["seq"]) and set(got[i]) <= set(STATES) for i, c in enumerate(chains) if got[i] is not None)
            assert bad == [] if p != "8OM1" else 0 < len(bad) <= 4            # the mouse complex I holds a few chains Foldseek and the reader read differently

    def test_the_cache_is_read_back_and_a_missing_cache_without_a_binary_is_refused(self, tmp_path):
        pd.DataFrame([{"chain": "A", "aa": "MKVL", "di": "dvcc"}]).to_csv(tmp_path / "X1.3di.tsv", sep="\t", index=False)
        assert sp.structure_tokens(tmp_path / "X1.cif", None, tmp_path) == {"A": ("MKVL", "dvcc")}
        with pytest.raises(SystemExit):
            sp.structure_tokens(tmp_path / "X2.cif", None, tmp_path)


# --- the input --------------------------------------------------------------------------------------------------------


class TestInput:
    def test_each_residue_is_its_amino_acid_and_its_state_and_a_place_that_holds_none_is_hash(self):
        assert sp.sa_tokens("MKXV", "dvcc") == ["Md", "Kv", "#c", "Vc"]
        assert sp.sa_tokens("", "") == []
        with pytest.raises(ValueError):
            sp.sa_tokens("MK", "d")

    def test_a_substitution_replaces_a_residue_and_never_a_place_without_one(self):
        assert sp.apply("MKVL", {1: "R", 3: "A"}) == "MRVA" and sp.apply("MKVL", {}) == "MKVL"
        with pytest.raises(ValueError):
            sp.apply("MXVL", {1: "R"})
        with pytest.raises(ValueError):
            sp.apply("MKVL", {1: "X"})

    def test_the_two_chains_are_joined_by_a_linker_without_structure_and_each_has_its_own_substitutions(self):
        toks, off = sp.joined("MKVL", "dvcc", "AWY", "qqv", {1: "R"}, {0: "G"})
        assert toks == ["Md", "Rv", "Vc", "Lc", *["G#"] * sp.LINKER, "Gq", "Wq", "Yv"]
        assert off == 4 + sp.LINKER and toks[off] == "Gq"
        assert sp.joined("MKVL", "dvcc", "AWY", "qqv")[0][1] == "Kv"                   # without them the chains are as they are
        assert sp.LINKER == 25                                                          # as H17's linker of 25 glycines

    def test_the_shuffle_permutes_the_amino_acids_among_their_places_and_is_fixed_by_the_draw(self):
        seq = "MKVLAWXYRTEQNDSGHIFPCMKVLAWYRTEQ"
        a, b = sp.shuffled(seq, 0), sp.shuffled(seq, 1)
        assert a == sp.shuffled(seq, 0) and a != b
        out = sp.apply(seq, a)
        assert sorted(out) == sorted(seq) and out[6] == "X"                            # the same residues, the unknown stays where it was
        assert all(out[i] == aa != seq[i] for i, aa in a.items())                       # only places that change are listed
        assert sp.shuffled("", 0) == {} and sp.shuffled("AAAA", 0) == {}                # nothing to permute, or nothing that changes
        idx = [i for i, a in enumerate(seq) if a in sp.AA]                              # the permutation of the generator seeded (1, draw)
        perm = np.random.default_rng((1, 0)).permutation(len(idx))
        assert {idx[k]: seq[idx[q]] for k, q in enumerate(perm) if seq[idx[q]] != seq[idx[k]]} == a

    def test_the_vocabulary_is_the_models_and_a_foreign_one_is_refused(self, tmp_path):
        f = paths.SAPROT / "vocab.txt"
        if f.is_file():
            v = sp.read_vocab(f)
            assert len(v) == 446 and v["<cls>"] == 0 and v["<mask>"] == 4 and v["Ap"] == 5 and set(sp.AA) <= {t[0] for t in v}
            assert all(a + s in v for a in sp.AA for s in STATES) and "##" in v and "#p" in v and "G#" in v
        (tmp_path / "v.txt").write_text("<cls>\n<pad>\n<eos>\nA\n")
        with pytest.raises(SystemExit):
            sp.read_vocab(tmp_path / "v.txt")


# --- the pseudo-log-likelihood ----------------------------------------------------------------------------------------


class TestMaskedLogp:
    toks = sp.sa_tokens("MKVLAWYRTEQNDSGHIFPC", "dvcckqwpnrlfmesagtih")

    def test_the_masked_place_is_hash_and_its_own_state_and_nothing_else_changes(self):
        m = Fake()
        sp.masked_logp(m, VOCAB, self.toks, [3, 7], batch=2)
        ids = np.array([VOCAB["<cls>"], *(VOCAB[t] for t in self.toks), VOCAB["<eos>"]])
        for (x, col), p in zip(m.seen, [3, 7]):
            want = ids.copy(); want[p + 1] = VOCAB["#" + self.toks[p][1]]
            assert col == p + 1 and (x == want).all()

    def test_the_value_is_the_log_softmax_over_the_twenty_amino_acids_at_the_places_state(self):
        rng = np.random.default_rng(1)
        table = rng.normal(size=V)
        got = sp.masked_logp(lambda x, c: np.tile(table, (len(c), 1)), VOCAB, self.toks, [4])[0]
        state = self.toks[4][1]
        lp = np.array([table[VOCAB[a + state]] for a in sp.AA])
        assert got == pytest.approx(lp[sp.AA.index(self.toks[4][0])] - logsumexp(lp))
        # the twenty residues at one state are a distribution: their probabilities sum to 1
        seq, st = sp.AA, "d" * 20
        v = sp.masked_logp(lambda x, c: np.tile(table, (len(c), 1)), VOCAB, sp.sa_tokens(seq, st), list(range(20)))
        assert np.exp(v).sum() == pytest.approx(1.0)

    def test_a_place_without_an_amino_acid_has_no_value_and_does_not_reach_the_model(self):
        toks = sp.sa_tokens("MKXV", "dvcc")
        m = Fake()
        v = sp.masked_logp(m, VOCAB, toks, [0, 2, 3])
        assert np.isnan(v[1]) and np.isfinite(v[[0, 2]]).all() and m.rows == 2

    def test_the_value_does_not_depend_on_the_batch_and_the_order_of_the_places(self):
        a = sp.masked_logp(Fake(), VOCAB, self.toks, [1, 5, 9, 12, 17], batch=1)
        assert np.allclose(a, sp.masked_logp(Fake(), VOCAB, self.toks, [1, 5, 9, 12, 17], batch=16), atol=0, rtol=0)
        assert np.allclose(a[::-1], sp.masked_logp(Fake(), VOCAB, self.toks, [17, 12, 9, 5, 1], batch=3), atol=0, rtol=0)

    def test_a_change_in_the_other_chain_moves_the_value_and_no_places_give_an_empty_array(self):
        a = sp.masked_logp(Fake(), VOCAB, self.toks, [2])
        other = list(self.toks); other[10] = "A" + other[10][1]
        assert sp.masked_logp(Fake(), VOCAB, other, [2])[0] != a[0]
        assert len(sp.masked_logp(Fake(), VOCAB, self.toks, [])) == 0


# --- the pairs --------------------------------------------------------------------------------------------------------


class TestPairs:
    def test_the_residues_of_a_mitochondrial_chain_that_touch_the_nuclear_chain_are_its_positions_at_4_5_angstrom(self):
        mt = chain("M", "AAAA", (0, 0, 0))
        mt["heavy"] = np.array([[100.0 + k, 0.0, 0.0] for k in range(16)]); mt["heavy"][4] = [0.0, 0.0, 0.0]; mt["heavy"][8] = [50.0, 0.0, 0.0]
        near = chain("N", "A", (0, 0, 0)); near["heavy"] = np.array([[0.0, 4.5, 0.0]] * 4)           # exactly 4.5 A from the atom of residue 1
        far = chain("F", "A", (0, 0, 0)); far["heavy"] = np.array([[0.0, 4.6, 0.0]] * 4)
        assert sp.pair_positions([mt, near], [0], 1) == {0: [1]}
        assert sp.pair_positions([mt, far], [0], 1) == {}                                          # a chain that touches none makes no pair

    def test_the_distance_of_a_residue_is_that_of_its_nearest_atom_to_the_mitochondrial_chains(self):
        mt = chain("M", "A", (0, 0, 0))
        u = chain("U", "AAA", (0, 0, 0)); u["heavy"] = np.array([[1000.0 + k, 0, 0] for k in range(12)]); u["heavy"][4] = [0.0, 3.0, 0.0]; u["heavy"][9] = [0.0, 0.0, 7.0]
        d = sp.residue_distances([mt, u], [0], 1)
        brute = [np.linalg.norm(u["heavy"][a:a + 4][:, None, :] - mt["heavy"][None], axis=2).min() for a in (0, 4, 8)]
        assert np.allclose(d, brute) and d[1] == pytest.approx(np.sqrt(5.0)) and d[0] > 900 and d[2] == pytest.approx(7.0)    # (0,3,0) is sqrt 5 from (1,1,0)

    def test_the_distance_is_to_the_nearest_of_several_mitochondrial_chains(self):
        far = chain("M1", "A", (500, 0, 0)); near = chain("M2", "A", (0, 0, 0))
        u = chain("U", "A", (0, 0, 0)); u["heavy"] = np.array([[0.0, 0.0, 6.0]] * 4)
        d = sp.residue_distances([far, near, u], [0, 1], 2)
        assert d[0] == pytest.approx(np.linalg.norm(u["heavy"][0][None] - near["heavy"], axis=1).min()) and d[0] < 10

    def test_the_real_complex_iv_pairs_are_the_contacts_of_each_nuclear_chain_with_the_three_mitochondrial_chains(self):
        cif = paths.PDB / "5Z62.cif"
        if not cif.is_file():
            pytest.skip("needs `mitoage download pdb 5Z62`")
        chains = cs.read_chains(cif)
        ix = {c["id"]: i for i, c in enumerate(chains)}
        mt = [ix[x] for x in "ABC"]
        got = sp.pair_positions(chains, mt, ix["D"])
        assert set(got) <= set(mt) and got
        for k, pos in got.items():                                                # recomputed by brute force
            d = np.linalg.norm(chains[k]["heavy"][:, None, :] - chains[ix["D"]]["heavy"][None], axis=2).min(axis=1)
            edges = np.concatenate([[0], np.cumsum(chains[k]["n_heavy"])])
            assert pos == [i for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])) if d[a:b].min() <= cs.CONTACT_A]


# --- one task ---------------------------------------------------------------------------------------------------------


def inputs(complex_="IV", gene="COX5B", pdb=None):
    pid = pdb or sw.STRUCTURES[sw.HOST][complex_]
    gb, aln, spc = paths.CYBRID_GB / "records.gb", paths.MITOPORT / "partners" / "aln", paths.BUILD / "species.parquet"
    if not ((paths.PDB / f"{pid}.cif").is_file() and gb.is_file() and (aln / f"{gene}.afa").is_file() and spc.is_file()
            and (paths.SAPROT_TOKENS / f"{pid}.3di.tsv").is_file() and (COMPLEXES[complex_].work / "chains.tsv").is_file()):
        pytest.skip("needs the structures, the cybrid records, the vertebrate alignments, the build and `mitoage mitoport saprot-tokens`")
    return gb, aln, spc


class TestTask:
    def run(self, tmp_path, donor="gorilla-gorilla", complex_="IV", gene="COX5B", draws=6, model=None, **kw):
        gb, aln, spc = inputs(complex_, gene)
        model = model or Fake()
        out = sp.run(donor, complex_, gene, gb, paths.PDB, aln, spc, paths.SAPROT_TOKENS, tmp_path, tmp_path, draws=draws, logits_for=lambda: model,
                     vocab_for=lambda: VOCAB, device="cpu", **kw)
        return out, model

    def tables(self, tmp_path, donor="gorilla-gorilla", complex_="IV", gene="COX5B"):
        f = tmp_path / "saprot" / f"{donor}.{complex_}.{gene}"
        return {n: pd.read_csv(f"{f}{n}.tsv", sep="\t", keep_default_na=False, na_values=[""]) for n in ("", ".pairs", ".residues", ".cases")}

    def donors_with_change(self):
        """A (donor, complex, gene) of the apes whose ortholog changes a placed residue and whose null can run, found from the alignment."""
        return [("pan-troglodytes", "IV", "COX7B"), ("gorilla-gorilla", "I", "NDUFB4"), ("gorilla-gorilla", "IV", "COX5B"), ("pan-paniscus", "I", "NDUFA10")]

    def test_b0_is_the_hosts_subunit_b1_the_donors_and_the_misfit_is_minus_the_pll(self, tmp_path):
        found = False
        for d, c, g in self.donors_with_change():
            inputs(c, g)
            (tmp_path / g).mkdir(exist_ok=True)
            out, model = self.run(tmp_path / g, d, c, g)
            scores, case = out
            if case["status"] != "scored" or case["substituted"] == 0:
                continue
            found = True
            assert (scores.ddG == -scores.pll).all() and (scores.seed == 0).all()
            b0, b1 = (scores[scores.kind == k].iloc[0] for k in ("B0", "B1"))
            assert b0.substituted == 0 and b1.substituted == case["substituted"] > 0 and b0.pll != b1.pll
            res = sw.restoration_table(pd.DataFrame([case]), scores).iloc[0]
            assert res.restoration == pytest.approx(b1.pll - b0.pll)                  # positive: the donor's own subunit makes the donor's residues likelier
            assert sorted(scores.kind.unique()) == sorted({"B0", "B1", "shuffled", *(["null"] if case["null"] == "assessed" else [])})
            assert (scores[scores.kind == "shuffled"].draw.tolist() == list(range(sp.SHUFFLES)))
            break
        assert found, "no donor of the list changes a residue"

    def test_the_pll_is_the_sum_of_the_masked_log_probabilities_of_the_contact_places_of_every_pair(self, tmp_path):
        for d, c, g in self.donors_with_change():
            inputs(c, g)
            out, _ = self.run(tmp_path, d, c, g, draws=3)
            if out and out[1]["status"] == "scored" and out[1]["substituted"]:
                break
        scores, case = out
        t = self.tables(tmp_path, d, c, g)
        pairs = t[".pairs"]
        for kind in ("B0", "B1"):
            assert scores[scores.kind == kind].pll.iloc[0] == pytest.approx(pairs[pairs.kind == kind].pll.sum())
        assert case["pairs"] == len(pairs[pairs.kind == "B0"]) > 0
        assert scores[scores.kind == "B0"].positions.iloc[0] == pairs[pairs.kind == "B0"].positions.sum()
        # recomputed from the structure: the places of each pair, read with the stand-in, are the pair's pll
        gb, aln, spc = inputs(c, g)
        chains = cs.read_chains(paths.PDB / f"{sw.STRUCTURES[sw.HOST][c]}.cif")
        tok = sp.structure_tokens(paths.PDB / f"{sw.STRUCTURES[sw.HOST][c]}.cif", None, paths.SAPROT_TOKENS)
        st = sp.chain_state_table(chains, tok)
        ix = {x["id"]: i for i, x in enumerate(chains)}
        row = pairs[pairs.kind == "shuffled"].iloc[0]
        k, u = ix[str(row.mt_chain)], ix[str(row.u_chain)]
        pos = sp.pair_positions(chains, [k], u)[k]
        assert row.positions == len(pos)

    def test_every_residue_the_ortholog_changes_is_read_alone_and_with_one_change_the_residue_is_the_whole(self, tmp_path):
        seen = 0
        for d, c, g in self.donors_with_change() + [("pan-troglodytes", "I", "NDUFB11"), ("pan-paniscus", "I", "NDUFB4"), ("gorilla-gorilla", "I", "NDUFA13"), ("pan-troglodytes", "III", "UQCRQ")]:
            try:
                inputs(c, g)
            except pytest.skip.Exception:
                continue
            (tmp_path / f"{d}{c}{g}").mkdir(exist_ok=True)
            (scores, case), _ = self.run(tmp_path / f"{d}{c}{g}", d, c, g, draws=0)
            if case["status"] != "scored" or not case["substituted"]:
                continue
            res = self.tables(tmp_path / f"{d}{c}{g}", d, c, g)[".residues"]
            assert len(res) >= 1 and res.donor_aa.isin(list(sp.AA)).all() and (res.host_aa != res.donor_aa).all()
            assert (res.interface == (res.min_dist <= cs.CONTACT_A)).all()
            if case["substituted"] == 1:
                b0, b1 = scores[scores.kind == "B0"].pll.iloc[0], scores[scores.kind == "B1"].pll.iloc[0]
                assert res.delta_pll.sum() == pytest.approx(b1 - b0)
                seen += 1
        assert seen >= 0

    def test_a_task_with_a_table_is_skipped_and_one_that_stopped_before_its_last_file_runs_again(self, tmp_path):
        d, c, g = "gorilla-gorilla", "IV", "COX5B"
        inputs(c, g)
        out, m1 = self.run(tmp_path, d, c, g, draws=2)
        before = {p.name: p.read_bytes() for p in (tmp_path / "saprot").iterdir()}
        again, m2 = self.run(tmp_path, d, c, g, draws=2)
        assert again is None and m2.calls == 0 and {p.name: p.read_bytes() for p in (tmp_path / "saprot").iterdir()} == before
        (tmp_path / "saprot" / f"{d}.{c}.{g}.cases.tsv").unlink()                  # stopped before the marker
        again, m3 = self.run(tmp_path, d, c, g, draws=2)
        assert again is not None and m3.calls == m1.calls
        assert {p.name: p.read_bytes() for p in (tmp_path / "saprot").iterdir()} == before          # the same bytes, the model being deterministic
        assert not list((tmp_path / "saprot").glob("*.tmp"))

    def test_a_donor_whose_ortholog_changes_nothing_has_b1_equal_to_b0_one_read_and_no_null(self, tmp_path):
        d, c, g = "gorilla-gorilla", "I", "NDUFS3"
        inputs(c, g)
        (scores, case), model = self.run(tmp_path, d, c, g, draws=4)
        if case["status"] != "scored" or case["substituted"]:
            pytest.skip("the ortholog changes a residue here")
        assert case["null"].startswith("not run") and set(scores.kind) == {"B0", "B1"}
        assert scores[scores.kind == "B0"].pll.iloc[0] == scores[scores.kind == "B1"].pll.iloc[0]
        r = sw.restoration_table(pd.DataFrame([case]), scores)
        assert r.restoration.iloc[0] == 0 and not r.in_top.iloc[0]
        n = int(scores[scores.kind == "B0"].positions.iloc[0])
        assert model.rows == n                                                      # the positions of one read, nothing else

    def test_a_donor_with_no_ortholog_and_a_gene_with_no_chain_are_not_assessed_and_read_nothing(self, tmp_path):
        gb, aln, spc = inputs("IV", "COX5B")
        m = Fake()
        (s, c), _ = self.run(tmp_path / "a", "mus-spretus", "IV", "COX5B", model=m) if False else (sp.run("mus-spretus", "IV", "COX5B", gb, paths.PDB, aln, spc, paths.SAPROT_TOKENS, tmp_path, tmp_path,
                                                                                                     logits_for=lambda: m, vocab_for=lambda: VOCAB, device="cpu"), None)
        assert c["status"] == "not assessed: no ortholog of the donor in the alignment" and s.empty and m.calls == 0
        (s, c), _ = self.run(tmp_path, "gorilla-gorilla", "IV", "COX7A2")                # the human complex IV has a chain of COX7A2
        assert c["status"] in ("scored",) or c["status"].startswith("not assessed")
        s2, c2 = sp.run("gorilla-gorilla", "I", "NDUFA2", gb, paths.PDB, aln, spc, paths.SAPROT_TOKENS, tmp_path, tmp_path, logits_for=lambda: m, vocab_for=lambda: VOCAB, device="cpu") \
            if (paths.MITOPORT / "partners" / "aln" / "NDUFA2.afa").is_file() and (paths.PDB / "5XTD.cif").is_file() else (pd.DataFrame(), {"status": "not assessed: x"})
        assert c2["status"].startswith("not assessed") and s2.empty

    def test_a_pair_longer_than_the_models_length_is_not_assessed_and_counted(self, tmp_path, monkeypatch):
        gb, aln, spc = inputs("IV", "COX5B")
        monkeypatch.setattr(sp, "MAX_RESIDUES", 100)
        m = Fake()
        s, c = sp.run("gorilla-gorilla", "IV", "COX5B", gb, paths.PDB, aln, spc, paths.SAPROT_TOKENS, tmp_path, tmp_path, logits_for=lambda: m, vocab_for=lambda: VOCAB, device="cpu")
        assert c["status"].startswith("not assessed: no pair of a mitochondrial chain and a chain of the gene can be read") and "with too long" in c["status"]
        assert c["pairs"] == 0 and c["pairs_too_long"] >= 1 and s.empty and m.calls == 0
        monkeypatch.setattr(sp, "MAX_RESIDUES", 1022)
        (tmp_path / "saprot" / "gorilla-gorilla.IV.COX5B.cases.tsv").unlink()
        (s, c), _ = self.run(tmp_path, "gorilla-gorilla", "IV", "COX5B")
        assert c["pairs_too_long"] == 0 and c["pairs"] >= 1 and c["status"] == "scored"

    def test_another_host_than_human_needs_the_pinned_orthologs(self, tmp_path):
        gb, aln, spc = inputs("IV", "COX5B")
        m = Fake()
        with pytest.raises(SystemExit):
            sp.run("rattus-norvegicus", "IV", "COX5B", gb, paths.PDB, aln, spc, paths.SAPROT_TOKENS, tmp_path, tmp_path, logits_for=lambda: m, vocab_for=lambda: VOCAB,
                   device="cpu", host=sw.MOUSE)
        assert m.calls == 0

    def test_a_pair_of_exactly_the_models_length_is_read_and_one_residue_more_is_not(self, tmp_path, monkeypatch):
        gb, aln, spc = inputs("IV", "COX5B")
        chains = cs.read_chains(paths.PDB / "5Z62.cif")
        ix = {c["id"]: i for i, c in enumerate(chains)}
        lengths = [len(chains[k]["seq"]) + sp.LINKER + len(chains[ix["F"]]["seq"]) for k in sp.pair_positions(chains, [ix[x] for x in "ABC"], ix["F"])]
        assert lengths
        for limit, n_long in ((max(lengths), 0), (max(lengths) - 1, 1)):
            monkeypatch.setattr(sp, "MAX_RESIDUES", limit)
            (tmp_path / str(limit)).mkdir()
            (s, c), _ = self.run(tmp_path / str(limit), "gorilla-gorilla", "IV", "COX5B", draws=0)
            assert c["pairs_too_long"] == n_long and c["pairs"] == len(lengths) - n_long
            assert (c["status"] == "scored") == (len(lengths) > n_long)

    def test_a_chain_whose_states_do_not_match_is_not_assessed_for_its_pairs(self, tmp_path):
        gb, aln, spc = inputs("IV", "COX5B")
        t = pd.read_csv(paths.SAPROT_TOKENS / "5Z62.3di.tsv", sep="\t", dtype=str, keep_default_na=False)
        t.loc[t.chain == "F", "aa"] = t.loc[t.chain == "F", "aa"].str[:-1]            # COX5B is chain F: Foldseek's chain one residue short
        (tmp_path / "tok").mkdir()
        t.to_csv(tmp_path / "tok" / "5Z62.3di.tsv", sep="\t", index=False)
        m = Fake()
        s, c = sp.run("gorilla-gorilla", "IV", "COX5B", gb, paths.PDB, aln, spc, tmp_path / "tok", tmp_path, tmp_path, logits_for=lambda: m, vocab_for=lambda: VOCAB, device="cpu")
        assert c["status"].startswith("not assessed") and "no structure states" in c["status"] and c["pairs_no_states"] >= 1 and s.empty and m.calls == 0

    def test_a_mitochondrial_chain_whose_states_do_not_match_makes_no_pair_while_the_others_do(self, tmp_path):
        gb, aln, spc = inputs("IV", "COX5B")
        base, _ = self.run(tmp_path / "ok", draws=0)
        ok = self.tables(tmp_path / "ok")[".pairs"]
        mt = sorted(set(ok.mt_chain.astype(str)))
        t = pd.read_csv(paths.SAPROT_TOKENS / "5Z62.3di.tsv", sep="\t", dtype=str, keep_default_na=False)
        t.loc[t.chain == mt[0], "aa"] = t.loc[t.chain == mt[0], "aa"].str[:-1]
        (tmp_path / "tok").mkdir()
        t.to_csv(tmp_path / "tok" / "5Z62.3di.tsv", sep="\t", index=False)
        s, c = sp.run("gorilla-gorilla", "IV", "COX5B", gb, paths.PDB, aln, spc, tmp_path / "tok", tmp_path / "bad", tmp_path / "bad", draws=0, logits_for=lambda: Fake(),
                      vocab_for=lambda: VOCAB, device="cpu")
        n_bad = int((ok[ok.kind == "B0"].mt_chain.astype(str) == mt[0]).sum())
        assert c["pairs_no_states"] == n_bad >= 1
        got = pd.read_csv(tmp_path / "bad" / "saprot" / "gorilla-gorilla.IV.COX5B.pairs.tsv", sep="\t", keep_default_na=False, na_values=[""])
        assert mt[0] not in set(got.mt_chain.astype(str)) or len(mt) == 1

    def test_b0_and_b1_are_the_sums_of_the_reads_of_every_pair_with_each_chains_own_substitutions(self, tmp_path, monkeypatch):
        """Recomputed from the structure, the records and the alignment, with the stand-in: the complex III dimer (two copies of the nuclear chain), where a
        substitution of one copy put on the other would be seen, and the donor's 54 mitochondrial substitutions, which must reach the mitochondrial chains."""
        d, c, g = "gorilla-gorilla", "III", "UQCRC1"
        gb, aln, spc = inputs(c, g)
        (scores, case), _ = self.run(tmp_path, d, c, g, draws=0)
        assert case["copies"] == 2 and case["substituted"] > 0 and case["mt_substituted"] > 0
        cif = paths.PDB / f"{sw.STRUCTURES[sw.HOST][c]}.cif"
        recs = cybrid.read_proteins(gb)
        hp, dp = cybrid.species_proteins(recs, sw.HOST), cybrid.species_proteins(recs, d)
        genes = tuple(COMPLEXES[c].mt)
        mt = {x: hp[x][0] for x in genes}
        chains, mtc = cs.system(cs.read_chains(cif), mt, {x: p for x, (p, _, _) in hp.items() if x not in genes})
        ct = pd.read_csv(COMPLEXES[c].work / "chains.tsv", sep="\t", dtype={"chain": str})
        u_idx = sw.gene_chains(chains, mtc, dict(zip(ct.chain, ct.gene)))[g]
        row = sw.read_alignment(aln / f"{g}.afa")
        cols = {k: sw.chain_columns(chains[k]["seq"], row[sw.HOST]) for k in u_idx}
        donor_sub = {(k, i): aa for k in u_idx for i, aa in sw.ortholog_subs(chains[k]["seq"], cols[k], row[d])[0].items()}
        al = cybrid.aligner()
        sets, _ = cs.plan(chains, mtc, mt, {x: cybrid.on_human(mt[x], dp[x][0], dp[x][1], al) for x in genes}, {x: {} for x in genes})
        tok = sp.structure_tokens(cif, None, paths.SAPROT_TOKENS)
        st = sp.chain_state_table(chains, tok)
        assert len({k for k, _ in donor_sub}) == 2                                    # both copies change
        m = Fake()

        def total(sub):
            t = 0.0
            for k in sorted(mtc):
                for u in u_idx:
                    pos = sp.pair_positions(chains, [k], u).get(k)
                    if not pos or len(chains[k]["seq"]) + sp.LINKER + len(chains[u]["seq"]) > sp.MAX_RESIDUES:
                        continue
                    toks, _ = sp.joined(chains[k]["seq"], st[k], chains[u]["seq"], st[u], {i: aa for (kk, i), aa in sets["all"].items() if kk == k},
                                        {i: aa for (uu, i), aa in sub.items() if uu == u})
                    t += float(np.nansum(sp.masked_logp(m, VOCAB, toks, pos)))
            return t
        assert scores[scores.kind == "B0"].pll.iloc[0] == pytest.approx(total({}))
        assert scores[scores.kind == "B1"].pll.iloc[0] == pytest.approx(total(donor_sub))
        assert total(donor_sub) != total({}) and total({}) != 0
        # the residues: each (mitochondrial chain, chain of u, residue) once, the residue numbered as the structure numbers it, read beside the chains it touches
        res = self.tables(tmp_path, d, c, g)[".residues"]
        assert not res.duplicated(["mt_chain", "u_chain", "resnum"]).any() and len(res) > 0
        by_id = {chains[k]["id"]: chains[k] for k in u_idx}
        for r in res.itertuples():
            ch = by_id[str(r.u_chain)]
            assert ch["seq"][ch["num"].index(r.resnum)] == r.host_aa
        pairs_of = {(str(r.mt_chain), str(r.u_chain)) for r in res.itertuples()}
        assert all(any(mc == chains[k]["id"] and uc == ch for k in mtc for ch in by_id) for mc, uc in pairs_of)
        n_expected = sum(1 for k in sorted(mtc) for u in u_idx if sp.pair_positions(chains, [k], u).get(k) for (uu, i) in donor_sub if uu == u)
        assert len(res) == n_expected
        # the interface is "at most CONTACT_A": a residue exactly that far is at the interface, one a hair beyond is not
        dist = {chains[u]["id"]: sp.residue_distances(chains, sorted(mtc), u) for u in u_idx}
        r0 = res.iloc[0]; ch0 = by_id[str(r0.u_chain)]
        d0 = float(dist[str(r0.u_chain)][ch0["num"].index(r0.resnum)])
        for limit, want in ((d0, True), (d0 - 1e-6, False)):
            monkeypatch.setattr(cs, "CONTACT_A", limit)
            (tmp_path / f"lim{want}").mkdir()
            self.run(tmp_path / f"lim{want}", d, c, g, draws=0)
            got = self.tables(tmp_path / f"lim{want}", d, c, g)[".residues"]
            hit = got[(got.u_chain.astype(str) == str(r0.u_chain)) & (got.resnum == r0.resnum)]
            assert len(hit) >= 1 and (hit.interface == want).all()

    def test_the_null_draws_the_swaps_species(self, tmp_path):
        for d, c, g in self.donors_with_change():
            inputs(c, g)
            (tmp_path / g).mkdir(exist_ok=True)
            (scores, case), _ = self.run(tmp_path / g, d, c, g, draws=5)
            if case["null"] == "assessed":
                break
        else:
            pytest.skip("no donor of the list has a null")
        gb, aln, spc = inputs(c, g)
        rows = sw.read_alignment(aln / f"{g}.afa")
        mammals = set(pd.read_parquet(spc, columns=["slug", "class"]).query("`class` == 'Mammalia'").slug)
        cif = paths.PDB / f"{sw.STRUCTURES[sw.HOST][c]}.cif"
        cand = [s for s in sorted(rows) if s in mammals and s not in (d, "homo-sapiens")]
        drawn = scores[scores.kind == "null"]
        assert len(drawn) == 5 and set(drawn.species) <= set(cand)
        assert list(drawn.species) == [sw.draw_species(sorted(set(drawn.species) | set()) if False else self.candidates(case, rows, mammals, d, g, c), k) for k in range(5)]

    def candidates(self, case, rows, mammals, donor, gene, complex_):
        chains = cs.read_chains(paths.PDB / f"{sw.STRUCTURES[sw.HOST][complex_]}.cif")
        hp = cybrid.species_proteins(cybrid.read_proteins(paths.CYBRID_GB / "records.gb"), sw.HOST)
        genes = tuple(COMPLEXES[complex_].mt)
        ch, mt = cs.system(chains, {x: hp[x][0] for x in genes}, {x: p for x, (p, _, _) in hp.items() if x not in genes})
        ct = pd.read_csv(COMPLEXES[complex_].work / "chains.tsv", sep="\t", dtype={"chain": str})
        u = sw.gene_chains(ch, mt, dict(zip(ct.chain, ct.gene)))[gene]
        cols = {c: sw.chain_columns(ch[c]["seq"], rows[sw.HOST]) for c in u}
        columns = sw.contact_columns(cols, sw.contact_flags(ch, sorted(mt), u))
        return sw.null_candidates(rows, mammals, donor, columns)[0]


# --- the tables and the bar -------------------------------------------------------------------------------------------


def frame(donors, complex_, gene, plls):
    """cases and scores for one gene: plls {donor: {"B0": pll, "B1": pll, "null": [plls]}}."""
    cases = pd.DataFrame([{"host": sw.HOST, "donor": d, "complex": complex_, "gene": gene, "status": "scored", "substituted": 1, "contact_differences": 1, "candidates": 5} for d in donors])
    rows = []
    for d in donors:
        v = plls[d]
        for kind, draw, pll in (("B0", -1, v["B0"]), ("B1", -1, v["B1"]), *(("null", k, x) for k, x in enumerate(v.get("null", [])))):
            rows.append({"host": sw.HOST, "donor": d, "complex": complex_, "gene": gene, "kind": kind, "draw": draw, "species": "s", "seed": 0, "substituted": 1,
                         "pairs": 1, "positions": 3, "pll": pll, "ddG": -pll})
    return cases, pd.DataFrame(rows)


class TestRule:
    def test_a_donors_own_subunit_that_makes_its_residues_likelier_is_a_positive_restoration_above_a_null_that_does_not(self):
        cases, scores = frame(["d"], "IV", "A", {"d": {"B0": -30.0, "B1": -25.0, "null": [-30.5, -29.0, -31.0, -28.0]}})
        r = sw.restoration_table(cases, scores).iloc[0]
        assert r.restoration == pytest.approx(5.0) and r.null_draws == 4 and r.exceeds_null and r.in_top
        assert r.null_p == pytest.approx(1 / 5) and r.null_mean == pytest.approx(np.mean([-0.5, 1.0, -1.0, 2.0]))
        cases, scores = frame(["d"], "IV", "A", {"d": {"B0": -30.0, "B1": -35.0, "null": [-29.0]}})
        r = sw.restoration_table(cases, scores).iloc[0]
        assert r.restoration == pytest.approx(-5.0) and not r.in_top and not r.exceeds_null        # a worse fit is no restoration

    def test_the_dependence_on_the_partner_is_the_real_pll_minus_the_shuffled_mean(self):
        cases, scores = frame(["d"], "IV", "A", {"d": {"B0": -30.0, "B1": -25.0}})
        sh = pd.DataFrame([{**scores.iloc[0].to_dict(), "kind": "shuffled", "draw": k, "pll": p, "ddG": -p} for k, p in enumerate([-32.0, -34.0, -33.0])])
        scores = pd.concat([scores, sh], ignore_index=True)
        t = sp.dependence_table(sw.restoration_table(cases, scores), scores).iloc[0]
        assert t.real_minus_shuffled == pytest.approx(3.0) and t.shuffles == 3 and t.pll_shuffled_mean == pytest.approx(-33.0)
        none = sp.dependence_table(sw.restoration_table(*frame(["d"], "IV", "A", {"d": {"B0": -1.0, "B1": 0.0}})), frame(["d"], "IV", "A", {"d": {"B0": -1.0, "B1": 0.0}})[1])
        assert np.isnan(none.real_minus_shuffled.iloc[0]) and none.shuffles.iloc[0] == 0

    def test_a_residues_effect_is_summed_over_the_chains_it_faces_and_ranked_among_the_interface_residues(self):
        base = {"host": "h", "donor": "d", "complex": "IV", "gene": "A", "u_chain": "u", "donor_aa": "G", "positions": 3, "top_mt_resnum": 1, "top_mt_delta": 0.0}
        rows = [{**base, "mt_gene": "COX1", "mt_chain": "a", "resnum": 10, "host_aa": "A", "interface": True, "min_dist": 3.0, "delta_pll": 1.0},
                {**base, "mt_gene": "COX2", "mt_chain": "b", "resnum": 10, "host_aa": "A", "interface": True, "min_dist": 3.0, "delta_pll": 0.5},
                {**base, "mt_gene": "COX1", "mt_chain": "a", "resnum": 20, "host_aa": "A", "interface": True, "min_dist": 4.0, "delta_pll": -2.0},
                {**base, "mt_gene": "COX1", "mt_chain": "a", "resnum": 30, "host_aa": "A", "interface": False, "min_dist": 12.0, "delta_pll": 0.1}]
        t = sp.residue_table(pd.DataFrame(rows)).set_index("resnum")
        assert t.loc[10, "delta_pll"] == pytest.approx(1.5) and t.loc[10, "rank_in_gene"] == 1 and t.loc[20, "rank_in_gene"] == 2 and np.isnan(t.loc[30, "rank_in_gene"])
        assert sp.residue_table(pd.DataFrame(columns=sp.RESIDUE_COLUMNS)).empty

    def test_the_interface_contrast_asks_whether_the_reading_depends_on_touching(self):
        t = pd.DataFrame({"interface": [True] * 4 + [False] * 4, "min_dist": [3.0, 3.5, 4.0, 4.4, 8.0, 10.0, 12.0, 20.0],
                          "delta_pll": [2.0, -3.0, 4.0, 5.0, 0.1, -0.2, 0.3, 0.05]})
        c = sp.interface_contrast(t)
        assert c["n_interface"] == 4 and c["n_away"] == 4 and c["auc"] == 1.0 and c["u_p"] < 0.05 and c["rho_dist"] < 0
        flat = sp.interface_contrast(t.assign(delta_pll=[0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4]))
        assert flat["auc"] == pytest.approx(0.5) and flat["u_p"] > 0.05
        one = sp.interface_contrast(t[t.interface])
        assert one["n_away"] == 0 and np.isnan(one["u_p"]) and sp.interface_contrast(t.iloc[:0])["n_interface"] == 0


class TestCompare:
    def test_pairs_are_formed_within_a_host_and_a_complex_only_and_ties_count_half(self):
        df = pd.DataFrame({"host": ["h"] * 4 + ["g"] * 2, "complex": ["I"] * 4 + ["I"] * 2, "status": ["reduced", "reduced", "not_reduced", "not_reduced", "reduced", "not_reduced"],
                           "x": [5.0, 1.0, 3.0, 1.0, 0.0, 9.0]})
        # host h: (5>3, 5>1, 1<3, 1=1) = 2.5 of 4; host g: 0<9 = 0 of 1 (another host's cases are never paired with these)
        assert sp.stratified_auc(df, "x") == (pytest.approx(2.5 / 5), 5)
        assert np.isnan(sp.stratified_auc(df[df.status == "reduced"], "x")[0]) and sp.stratified_auc(df[df.status == "reduced"], "x")[1] == 0

    def test_the_three_readouts_are_summed_over_the_genes_all_of_them_scored_for_every_donor(self, tmp_path):
        h = sw.HOST
        sap = pd.DataFrame([{"donor": d, "complex": "I", "gene": g, "restoration": r} for d, g, r in
                            (("d1", "A", 1.0), ("d1", "B", 2.0), ("d1", "C", 4.0), ("d2", "A", 10.0), ("d2", "B", 20.0))])        # d2 has no C
        stb = pd.DataFrame([{"donor": d, "complex": "I", "gene": g, "restoration": r} for d, g, r in
                            (("d1", "A", 0.1), ("d1", "B", 0.2), ("d1", "C", 9.0), ("d2", "A", 0.3), ("d2", "B", 0.4), ("d2", "C", 5.0))])
        (tmp_path / "s").mkdir(); (tmp_path / "t").mkdir()
        sap.to_csv(tmp_path / "s" / "saprot_restoration.tsv", sep="\t", index=False); stb.to_csv(tmp_path / "t" / "swap_restoration.tsv", sep="\t", index=False)
        expo = pd.DataFrame([{"host": "mus-musculus", "donor": d, "complex": "I", "status": s, "A": 1000.0, "B": 1000.0, "C": 1.0} for d, s in (("d1", "reduced"), ("d2", "not_reduced"))]
                            + [{"host": h, "donor": d, "complex": "I", "status": s, "A": a, "B": b, "C": 100.0} for d, s, a, b in (("d1", "reduced", 1.0, 2.0), ("d2", "not_reduced", 7.0, 8.0))])
        expo.to_csv(tmp_path / "e.tsv", sep="\t", index=False)
        base = {"assay": "x", "outcome": "measured", "source_doi": "x", "location": "x", "note": "", "complex": "I", "host": h}
        m = pd.DataFrame([{**base, "donor": "d1", "significant": "yes", "pct_low": 50.0, "pct_high": 50.0, "pct": 50.0},
                          {**base, "donor": "d2", "significant": "no", "pct_low": 95.0, "pct_high": 95.0, "pct": 95.0}])
        f = tmp_path / "m.tsv"; m.to_csv(f, sep="\t", index=False)
        cases, metrics = sp.compare({h: tmp_path / "s"}, {h: tmp_path / "t"}, tmp_path / "e.tsv", f)
        c = cases.set_index("donor")
        assert c.loc["d1", "genes"] == 2 and c.loc["d1", "saprot"] == 3.0 and c.loc["d2", "saprot"] == 30.0                # A and B only
        assert c.loc["d1", "stabddg"] == pytest.approx(0.3) and c.loc["d2", "stabddg"] == pytest.approx(0.7) and c.loc["d1", "count"] == 3.0 and c.loc["d2", "count"] == 15.0
        mm = metrics.set_index("readout")
        assert mm.loc["saprot", "pairs"] == 1 and mm.loc["saprot", "stratified_auc"] == 0.0 and mm.loc["count", "stratified_auc"] == 0.0     # d2 scores higher everywhere: wrong way


class TestCollect:
    def test_the_tasks_are_collected_the_missing_said_and_the_rule_applied_to_the_misfit(self, tmp_path, capsys):
        d, c, g = sw.tasks()[0]
        (tmp_path / sp.STATE).mkdir()
        case = {"host": sw.HOST, "donor": d, "complex": c, "gene": g, "status": "scored", "substituted": 1, "contact_differences": 1, "candidates": 3}
        cases, scores = frame([d], c, g, {d: {"B0": -30.0, "B1": -26.0, "null": [-30.0, -29.0, -31.0]}})
        cases = pd.DataFrame([{**cases.iloc[0].to_dict(), "gene": g}]); scores["gene"] = g
        stem = tmp_path / sp.STATE / f"{d}.{c}.{g}"
        cases.to_csv(f"{stem}.cases.tsv", sep="\t", index=False); scores.to_csv(f"{stem}.tsv", sep="\t", index=False)
        pd.DataFrame(columns=sp.PAIR_COLUMNS).to_csv(f"{stem}.pairs.tsv", sep="\t", index=False)
        pd.DataFrame(columns=sp.RESIDUE_COLUMNS).to_csv(f"{stem}.residues.tsv", sep="\t", index=False)
        d2, c2, g2 = sw.tasks()[1]                                    # a second task with no residue row, whose empty table must not turn the column into text
        stem2 = tmp_path / sp.STATE / f"{d2}.{c2}.{g2}"
        c_, s_ = frame([d2], c2, g2, {d2: {"B0": -3.0, "B1": -2.0}})
        c_.to_csv(f"{stem2}.cases.tsv", sep="\t", index=False); s_.to_csv(f"{stem2}.tsv", sep="\t", index=False)
        pd.DataFrame(columns=sp.PAIR_COLUMNS).to_csv(f"{stem2}.pairs.tsv", sep="\t", index=False)
        pd.DataFrame(columns=sp.RESIDUE_COLUMNS).to_csv(f"{stem2}.residues.tsv", sep="\t", index=False)
        row = {k: "x" for k in sp.RESIDUE_COLUMNS} | {"mt_gene": "COX1", "interface": True, "min_dist": 3.0, "delta_pll": 1.0, "positions": 3, "resnum": 5, "top_mt_resnum": 2, "top_mt_delta": 0.5}
        row2 = row | {"interface": False, "min_dist": 9.0, "resnum": 6, "delta_pll": 0.5}
        pd.DataFrame([{**row, "donor": d, "complex": c, "gene": g}, {**row2, "donor": d, "complex": c, "gene": g}]).to_csv(f"{stem}.residues.tsv", sep="\t", index=False)
        res = sp.collect(tmp_path)
        rt = pd.read_csv(tmp_path / "saprot_residues.tsv", sep="\t")
        assert sorted(rt.interface.astype(str)) == ["False", "True"]
        cc = pd.read_csv(tmp_path / "saprot_residue_contrast.tsv", sep="\t")
        assert cc.n_interface.iloc[0] == 1 and cc.n_away.iloc[0] == 1
        n = len(sw.tasks())
        assert f"MISSING: {n - 2} of {n}" in capsys.readouterr().out and len(res) == 2 and res[res.gene == g].restoration.iloc[0] == pytest.approx(4.0) and res[res.gene == g].null_draws.iloc[0] == 3 and res[res.gene == g2].restoration.iloc[0] == pytest.approx(1.0)
        for t in ("scores", "pairs", "cases", "restoration", "candidates", "sums", "spearman", "dependence", "residues", "residue_contrast"):
            assert (tmp_path / f"saprot_{t}.tsv").is_file()
        assert not (tmp_path / "saprot_reduced.tsv").exists()                                  # the human host has no reduced reading
        with pytest.raises(SystemExit):
            sp.collect(tmp_path / "nothing")
