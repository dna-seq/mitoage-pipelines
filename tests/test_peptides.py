"""H15 (RM16): the published peptides' frames across species, the nulls and the tree.

Real data: the build's alignments (`download build`). Expected values are computed at run time from the data or
from a fixture built here; no conservation result is asserted.
"""
import io
import itertools
import math
from collections import Counter

import numpy as np
import polars as pl
import pytest
from Bio import Phylo

from mitoage_pipelines import paths
from mitoage_pipelines.peptides import catalogue, codes, frames, nulls, shuffle, tree

needs_build = pytest.mark.skipif(not (paths.BUILD / "species.parquet").is_file(), reason="needs `download build`")
needs_tree = pytest.mark.skipif(not paths.DATED_TREE.is_file(), reason="needs `longevity dated-tree`")
BY_NAME = {p.name: p for p in catalogue.PEPTIDES}


@pytest.fixture(scope="module")
def meta():
    sp = frames.species_table(paths.BUILD)
    return {r["slug"]: r for r in sp.iter_rows(named=True)}


@pytest.fixture(scope="module")
def code_of(meta):
    return {s: int(r["transl_table"]) for s, r in meta.items()}


@pytest.fixture(scope="module")
def hosts():
    return {}


def placed(name, meta, hosts):
    return frames.placed(paths.BUILD, BY_NAME[name], set(meta), hosts)


def slow_state(seq: str, n: int, table: codes.Table) -> tuple[str, int]:
    """An independent reading of one species' frame from its strand-oriented ungapped sequence `seq` (from the
    window's first base): plain string translation, no arrays."""
    m = math.ceil(0.9 * n - 1e-9)
    aa = codes.translate(seq, table)
    for k, a in enumerate(aa[:m]):
        if a == "*":
            return "lost", k
        if a == "X":
            return "not assessed", k
    if len(aa) < m:
        return "not assessed", len(aa)
    return ("kept" if table.is_start(seq[:3]) else "lost"), m


# --- the catalogue -----------------------------------------------------------------------------------------------
@needs_build
@pytest.mark.parametrize("pep", catalogue.PEPTIDES, ids=lambda p: p.name)
def test_the_human_row_reads_the_published_peptide_and_then_a_stop(pep):
    from Bio.Seq import Seq
    from mitoage_pipelines.seqio import read_fasta
    r = catalogue.place(pep, paths.BUILD)
    h = read_fasta(paths.BUILD / "alignments" / f"{pep.alignment}.fasta.gz")["homo-sapiens"].upper().replace("-", "")
    s = h if pep.strand == "+" else str(Seq(h).reverse_complement())
    i = r["human_index"] if pep.strand == "+" else len(h) - 1 - r["human_index"]
    t = codes.tables(2, pep.reading)
    n = len(pep.sequence)
    assert t.is_start(s[i:i + 3])
    assert codes.translate(s[i + 3:i + 3 * n], t) == pep.sequence[1:]
    assert t.is_stop(s[i + 3 * n:i + 3 * n + 3]) and r["stop_after_is_stop"]
    if pep.nucleotides:
        assert s[i:i + len(pep.nucleotides)] == pep.nucleotides


@needs_build
def test_placements_agree_with_the_papers_coordinates():
    """MTALTND4, MTALTCO1 at the published start; MENTSH's start codon holds m.4977; the SHLPs and humanin are 2 and
    1 above rCRS (the paper's numbering), each SHLP by the same offset."""
    got = {p.name: catalogue.place(p, paths.BUILD) for p in catalogue.PEPTIDES}
    assert got["MTALTND4"]["rcrs_start"] == 11557 and got["MTALTCO1"]["rcrs_start"] == 6089
    assert got["MENTSH"]["rcrs_start"] <= 4977 <= got["MENTSH"]["rcrs_start"] + 2
    offsets = {k: int(v["published_coordinates"].split("-")[0]) - v["rcrs_start"]
               for k, v in got.items() if k.startswith("SHLP") or k == "humanin"}
    assert {offsets[k] for k in offsets if k.startswith("SHLP")} == {2}
    assert offsets["humanin"] == 1


def test_a_sequence_not_in_the_human_row_is_refused():
    if not (paths.BUILD / "species.parquet").is_file():
        pytest.skip("needs `download build`")
    bad = catalogue.Peptide(**{**BY_NAME["SHLP6"].__dict__, "sequence": "MLDQDIPMVQPLLKVRLFNW"})
    with pytest.raises(ValueError, match="read 0 times"):
        catalogue.place(bad, paths.BUILD)


def test_agr_is_arginine_only_in_the_agr_reading():
    assert codes.tables(2, "mito").is_stop("AGA") and not codes.tables(2, "mito_agr_r").is_stop("AGA")
    assert codes.tables(2, "own_agr_r").aa[codes.CODONS.index("AGG")] == "R"
    assert codes.tables(5, "own_agr_r").aa[codes.CODONS.index("AGG")] == codes.tables(5, "own").aa[codes.CODONS.index("AGG")]
    std = codes.tables(2, "standard")
    assert std.is_stop("TGA") and not std.is_stop("AGA") and {c for c in codes.CODONS if std.is_start(c)} == {"ATG"}
    assert codes.species_reading("own", "mito_agr_r") == "own_agr_r" and codes.species_reading("own", "mito") == "own"


# --- the frames --------------------------------------------------------------------------------------------------
@needs_build
@pytest.mark.parametrize("name", ["humanin", "SHLP3", "MTALTCO1", "MENTSH"])
@pytest.mark.parametrize("reading", frames.READINGS)
def test_the_array_reading_equals_a_plain_translation_in_every_species(name, reading, meta, code_of, hosts):
    pc = placed(name, meta, hosts)
    st, tabs = frames.observe(pc, reading, code_of)
    label = {frames.KEPT: "kept", frames.LOST: "lost", frames.NOT_ASSESSED: "not assessed"}
    checked = 0
    for i in range(len(pc.host.slugs)):
        if st["coverage"][i] < frames.COVERAGE_MIN:
            assert st["state"][i] == frames.NOT_ASSESSED and st["reason"][i] == "gap"
            continue
        p = int(st["p"][i])
        want, _ = slow_state(pc.host.sequence(i)[p:], pc.n, tabs[i])
        assert label[int(st["state"][i])] == want, pc.host.slugs[i]
        checked += 1
    assert checked > 1000


@needs_build
@pytest.mark.parametrize("pep", catalogue.PEPTIDES, ids=lambda p: p.name)
def test_human_keeps_every_peptide_under_its_published_reading(pep, meta, code_of, hosts):
    pc = placed(pep.name, meta, hosts)
    reading = "standard" if pep.reading == "standard" else "own"
    st, tabs = frames.observe(pc, reading, code_of)
    h = pc.host.row("homo-sapiens")
    assert st["state"][h] == frames.KEPT and st["orf"][h] == pc.n
    ident, aa = frames.identity(pc, h, int(st["p"][h]), tabs[h])
    assert ident == 1.0 and aa[1:] == pep.sequence[1:]


@needs_build
def test_a_stop_an_unknown_base_and_a_gap_in_the_human_window(meta, code_of, hosts):
    """Edit the human row in memory: a stop inside the frame makes it lost, an N before the decision makes it not
    assessed, gapping the window makes it not assessed; the other species do not move."""
    pc = placed("MTALTND4", meta, hosts)
    h = pc.host.row("homo-sapiens")
    base, _ = frames.observe(pc, "own", code_of)
    host = pc.host
    keep = (host.nt.copy(), host.nongap.copy(), host.cum.copy(), host.length.copy())
    try:
        p = int(base["p"][h])
        host.nt[h, p + 30:p + 33] = codes.to_index("TAA")              # codon 11: an early stop
        st, _ = frames.observe(pc, "own", code_of)
        assert st["state"][h] == frames.LOST and st["reason"][h] == "early stop" and st["orf"][h] == 10
        host.nt[h] = keep[0][h]
        host.nt[h, p + 31] = 4                                        # an unknown base in codon 11
        st, _ = frames.observe(pc, "own", code_of)
        assert st["state"][h] == frames.NOT_ASSESSED and st["reason"][h] == "ambiguous base"
        host.nt[h] = keep[0][h]
        host.nongap[h, pc.cols()[: 3 * pc.n // 2 + 3]] = False          # over half the window gapped
        st, _ = frames.observe(pc, "own", code_of)
        assert st["state"][h] == frames.NOT_ASSESSED and st["reason"][h] == "gap"
        others = np.arange(len(host.slugs)) != h
        assert (st["state"][others] == base["state"][others]).all()
    finally:
        host.nt[:], host.nongap[:], host.cum[:], host.length[:] = keep


def test_not_assessed_is_never_counted_as_kept_or_lost():
    st = pl.DataFrame({
        "peptide": ["x"] * 6, "reading": ["own"] * 6, "slug": list("abcdef"), "class": ["A"] * 6,
        "state": ["kept", "lost", "not assessed", "not assessed", "kept", "lost"],
        "reason": ["kept", "no start", "gap", "no host row", "kept", "early stop"],
        "open": [True, True, False, False, True, False], "start": [True, False, False, False, True, True],
        "identity": [1.0, None, None, None, 0.5, None]})
    k = frames.summarise(st).filter(pl.col("class") == "A").row(0, named=True)
    assert (k["kept"], k["lost"], k["not_assessed"], k["assessed"]) == (2, 2, 2, 4)
    assert k["kept"] + k["lost"] + k["not_assessed"] == k["species"]
    assert k["kept_share"] == 2 / 4 and k["open_share"] == 3 / 4 and k["start"] == 3 and k["with_region"] == 5
    assert k["status"].startswith("not assessed")                    # 4 assessed species, under 20


# --- the shuffles ------------------------------------------------------------------------------------------------
@needs_build
@pytest.mark.parametrize("slug", ["homo-sapiens", "gallus-gallus", "danio-rerio"])
def test_the_synonymous_shuffle_keeps_the_protein_and_the_codon_counts(slug, meta, code_of):
    from mitoage_pipelines.seqio import read_fasta
    seq = read_fasta(paths.BUILD / "alignments" / "COX1.codon.fasta.gz")[slug].upper().replace("-", "")
    cod = codes.codon_index(codes.to_index(seq))[0::3][:len(seq) // 3]
    t = codes.tables(code_of[slug], "own")
    aa = shuffle.residues(cod, t)
    draws = shuffle.synonymous(cod, aa, 50, shuffle.rng_for(1, "test", slug))
    protein = "".join(t.aa[c] if c >= 0 else "X" for c in cod)
    for row in draws:
        assert "".join(t.aa[c] if c >= 0 else "X" for c in row) == protein
        assert Counter(row.tolist()) == Counter(cod.tolist())
    assert (draws[:, aa < 0] == cod[aa < 0]).all()                    # fixed codons stay
    assert (draws != cod).any(1).all()                                # and every shuffle moved something


def test_the_shuffles_repeat_with_the_seed_and_differ_without():
    cod = np.array([0, 1, 2, 3, 16, 17, 18, 19] * 10, dtype=np.int16)
    t = codes.tables(2, "own")
    aa = shuffle.residues(cod, t)
    a = shuffle.synonymous(cod, aa, 20, shuffle.rng_for(7, "p", "s"))
    b = shuffle.synonymous(cod, aa, 20, shuffle.rng_for(7, "p", "s"))
    c = shuffle.synonymous(cod, aa, 20, shuffle.rng_for(7, "p", "t"))
    assert a.tobytes() == b.tobytes() and a.tobytes() != c.tobytes()
    s = "ATGGCATTAGCCTTAGGATAAACGT" * 3
    x = shuffle.dinucleotide_many(s, 30, shuffle.rng_for(7, "p", "s"))
    y = shuffle.dinucleotide_many(s, 30, shuffle.rng_for(7, "p", "s"))
    assert x.tobytes() == y.tobytes()


def dinucs(s):
    return Counter(zip(s, s[1:]))


def test_the_dinucleotide_shuffle_keeps_counts_ends_and_is_uniform():
    s = "AACGTTACGA"
    rng = shuffle.rng_for(3, "uniform")
    draws = [shuffle.dinucleotide(s, rng) for _ in range(6000)]
    assert all(dinucs(d) == dinucs(s) and d[0] == s[0] and d[-1] == s[-1] for d in draws)
    # every arrangement with the same dinucleotides, first and last base, by brute force
    pool = {"".join(p) for p in itertools.permutations(s)}
    target = {p for p in pool if dinucs(p) == dinucs(s) and p[0] == s[0] and p[-1] == s[-1]}
    seen = Counter(draws)
    assert set(seen) == target and len(target) > 2
    expected = len(draws) / len(target)
    chi2 = sum((seen[t] - expected) ** 2 / expected for t in target)
    from scipy import stats
    assert stats.chi2.sf(chi2, len(target) - 1) > 1e-4


def test_short_and_single_letter_sequences_are_returned_as_they_are():
    rng = shuffle.rng_for(1, "edge")
    assert shuffle.dinucleotide("", rng) == "" and shuffle.dinucleotide("AC", rng) == "AC"
    assert shuffle.dinucleotide("AAAAAA", rng) == "AAAAAA"


@needs_build
@pytest.mark.parametrize("name", ["MTALTND4", "MENTSH", "SHLP6", "SHLP2"])
def test_the_null_reads_the_unshuffled_frame_as_the_observation_does(name, meta, code_of, hosts):
    """The shuffle machinery, applied to the real codons, gives each species the state `frames` gives it."""
    pc = placed(name, meta, hosts)
    d = nulls.shuffle_draws(pc, code_of, 3, 1)
    for reading in frames.READINGS:
        st, _ = frames.observe(pc, reading, code_of)
        idx = [pc.host.slugs.index(s) for s in d["slugs"]]
        sub = st["state"][idx]
        ok = sub != frames.NOT_ASSESSED
        assert ok.sum() > 0.9 * len(idx)
        assert ((sub[ok] == frames.KEPT) == d[reading]["kept_obs"][ok]).all()
        assert (st["open"][idx][ok] == d[reading]["open_obs"][ok]).all()
    assert d["own"]["kept"].shape == (3, len(d["slugs"]))


@needs_build
def test_window_null_windows_are_orfs_in_human_and_never_the_peptide(meta, code_of, hosts):
    for name in ("humanin", "MTALTND4", "SHLP1"):
        pc = placed(name, meta, hosts)
        reading = "standard" if pc.pep.reading == "standard" else "own"
        t = codes.tables(2, codes.species_reading(reading, pc.pep.reading))
        seq = pc.host.sequence(pc.host.row("homo-sapiens"))
        kept = nulls.window_starts(pc, reading, 2, "kept")
        opened = nulls.window_starts(pc, reading, 2, "open")
        assert set(kept) <= set(opened) and pc.q0 not in opened
        for q in opened:
            assert "*" not in codes.translate(seq[q:q + 3 * pc.n], t)
            assert not ((q - pc.q0) % 3 == 0 and abs(q - pc.q0) < 3 * pc.n)
            if pc.pep.alignment.endswith(".codon"):
                assert q % 3 == pc.q0 % 3
        assert all(t.is_start(seq[q:q + 3]) for q in kept)


def test_representatives_are_one_per_family_and_fixed():
    meta = {s: {"family": f, "genus": "g"} for s, f in [("a", "F1"), ("b", "F1"), ("c", "F2"), ("d", None)]}
    r = nulls.representatives(["a", "b", "c", "d"], meta)
    assert len(r) == 3 and "c" in r and len(r & {"a", "b"}) == 1
    assert r == nulls.representatives(["d", "c", "b", "a"], meta)


def test_p_counts_ties_and_adds_one():
    pep = catalogue.PEPTIDES[0]
    meta = {f"s{i}": {"class": "A", "family": f"F{i}", "genus": "g"} for i in range(20)}
    obs = np.array([True] * 10 + [False] * 10)
    null = np.zeros((4, 20), dtype=bool); null[0, :10] = True; null[1, :12] = True     # 0.5, 0.6, 0, 0
    draws = {"slugs": list(meta), "own": {"kept_obs": obs, "open_obs": obs, "kept": null, "open": null}}
    draws["standard"] = draws["own"]
    rows = nulls.shuffle_rows(pep, draws, meta, ["A"], 4)
    r = next(x for x in rows if x["reading"] == "own" and x["statistic"] == "kept" and x["level"] == "species")
    assert r["observed"] == 0.5 and r["p"] == (1 + 2) / (1 + 4)


def test_q_is_benjamini_hochberg_within_null_statistic_level_and_reading():
    from statsmodels.stats.multitest import multipletests
    df = pl.DataFrame({"null": ["w", "w", "w", "s"], "statistic": ["kept"] * 4, "level": ["species"] * 4,
                       "reading": ["own"] * 4, "p": [0.01, 0.04, None, 0.03],
                       "status": ["tested", "tested", "not assessed", "tested"]})
    q = nulls.bh(df)["q"].to_list()
    assert q[:2] == pytest.approx(list(multipletests([0.01, 0.04], method="fdr_bh")[1]))
    assert q[2] is None and q[3] == pytest.approx(0.03)


# --- the tree ----------------------------------------------------------------------------------------------------
NEWICK = "((a:1,b:1,c:1):1,(d:1,(e:1,f:1):1):1,g:2);"


def brute_force(newick, states):
    """Minimum changes over every assignment of states to the internal nodes."""
    t = Phylo.read(io.StringIO(newick), "newick")
    internal = [c for c in t.find_clades() if c.clades]
    best = None
    for assign in itertools.product([0, 1], repeat=len(internal)):
        s = {id(c): v for c, v in zip(internal, assign)}
        for c in t.get_terminals():
            s[id(c)] = states.get(c.name)
        cost = 0
        for c in internal:
            for k in c.clades:
                if s[id(k)] is not None:
                    cost += s[id(k)] != s[id(c)]
        best = cost if best is None else min(best, cost)
    return best


@pytest.mark.parametrize("states", [
    {"a": 1, "b": 1, "c": 0, "d": 0, "e": 0, "f": 1, "g": 1},
    {"a": 0, "b": 0, "c": 0, "d": 1, "e": 1, "f": 1, "g": 0},
    {"a": 1, "b": 0, "c": 1, "e": 0, "g": 1},                                   # d and f missing
])
def test_sankoff_is_the_minimum_and_losses_plus_gains_equal_it(states):
    t = tree.Tree(NEWICK)
    best = brute_force(NEWICK, states)
    assert tree.sankoff(t, {k: np.array([v]) for k, v in states.items()})[0] == best
    changes, losses, gains, _ = tree.reconstruct(t, states)
    assert changes == best == losses + gains


def test_missing_tips_cost_nothing_and_clustering_is_seen():
    t = tree.Tree(NEWICK)
    assert tree.reconstruct(t, {"a": 1, "b": 1})[0] == 0
    clustered = {"a": 0, "b": 0, "c": 0, "d": 1, "e": 1, "f": 1, "g": 1}
    r = tree.test_group(t, clustered, 200, tree.seed_for(1, "x"))
    assert r["min_changes"] == 1 and r["perm_mean_changes"] > 1
    assert r["kept"] + r["lost"] == r["tips"] == 7
    r2 = tree.test_group(t, clustered, 200, tree.seed_for(1, "x"))
    assert r == r2


# --- reruns ------------------------------------------------------------------------------------------------------
@needs_build
def test_nulls_rerun_gives_identical_bytes_in_a_stable_order(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    nulls.main(paths.BUILD, a, replicates=5, peptides=["MENTSH", "SHLP6"])
    nulls.main(paths.BUILD, b, replicates=5, peptides=["SHLP6", "MENTSH"], jobs=2)
    assert (a / "nulls.tsv").read_bytes() == (b / "nulls.tsv").read_bytes()
    df = pl.read_csv(a / "nulls.tsv", separator="\t")
    keys = ["null", "statistic", "level", "reading", "peptide", "class"]
    assert df.select(keys).rows() == df.sort(keys).select(keys).rows()
    assert df.filter(pl.col("status") == "tested")["p"].is_between(0, 1).all()
    assert df.filter(pl.col("status") != "tested")["p"].null_count() == df.filter(pl.col("status") != "tested").height


def test_each_replicate_is_its_own_shuffle():
    cod = np.array(list(range(0, 64, 2)) * 6, dtype=np.int16)
    t = codes.tables(2, "own")
    d = shuffle.synonymous(cod, shuffle.residues(cod, t), 30, shuffle.rng_for(5, "r"))
    assert len({row.tobytes() for row in d}) == 30
    x = shuffle.dinucleotide_many("ATGCGTACGTTAGCATGCAAGT" * 3, 30, shuffle.rng_for(5, "d"))
    assert len({row.tobytes() for row in x}) > 25


# --- lifespan (H12's model) --------------------------------------------------------------------------------------
COMPOSITION = paths.LONGEVITY_OUT / "composition.tsv"
FOURFOLD = paths.LONGEVITY_OUT / "fourfold.tsv"
needs_composition = pytest.mark.skipif(not COMPOSITION.is_file() or not FOURFOLD.is_file() or not paths.DATED_TREE.is_file(),
                                       reason="needs `longevity controls`, `fourfold` and `dated-tree`")


@needs_composition
def test_lifespan_fits_kept_against_lost_and_leaves_out_not_assessed(tmp_path):
    """A states table made here (not a result): Aves species kept or lost by a hash of the slug, one in five not
    assessed; and one peptide with only three lost species."""
    import hashlib
    import pandas as pd
    from mitoage_pipelines.longevity import aa_pgls
    from mitoage_pipelines.peptides import lifespan
    sp, _ = aa_pgls.select(paths.BUILD, COMPOSITION, FOURFOLD, paths.DATED_TREE, "Aves")
    rows = []
    for i, s in enumerate(sorted(sp.slug)):
        h = int(hashlib.sha256(s.encode()).hexdigest(), 16)
        state = "not assessed" if h % 5 == 0 else ("kept" if h % 2 else "lost")
        rows.append({"peptide": "humanin", "reading": "standard", "slug": s, "state": state})
        rows.append({"peptide": "SHLP1", "reading": "standard", "slug": s, "state": "lost" if i < 3 else "kept"})
    f = tmp_path / "states.tsv"
    pd.DataFrame(rows).to_csv(f, sep="\t", index=False)
    res = lifespan.main(paths.BUILD, f, COMPOSITION, FOURFOLD, paths.DATED_TREE, tmp_path / "pgls.tsv", classes=("Aves",))
    st = pd.DataFrame(rows)
    hn = res[res.peptide == "humanin"].iloc[0]
    want = st[(st.peptide == "humanin") & (st.state != "not assessed")]
    assert hn.n == len(want) == hn.n_kept + hn.n_lost and hn.n_kept == (want.state == "kept").sum()
    assert hn.status == "tested" and hn.lo95 < hn.slope < hn.hi95 and 0 <= hn["lambda"] <= 1
    assert res[res.peptide == "SHLP1"].iloc[0].status.startswith("not assessed")
    assert set(res.peptide) == {p.name for p in catalogue.PEPTIDES}


def test_on_a_tie_the_child_keeps_its_parents_state():
    """Root and the polytomy (a, b, e) are both ties; the rule (root kept, a child keeps its parent's state on a tie)
    puts the changes at the tips: a and b lose the frame, nothing is gained."""
    t = tree.Tree("((a:1,b:1,e:1):1,c:2);")
    assert tree.reconstruct(t, {"a": 0, "b": 0, "e": 1, "c": 1}) == (2, 2, 0, 1)


@needs_build
def test_published_nucleotides_that_differ_from_the_human_row_are_refused():
    pep = BY_NAME["MENTSH"]
    wrong = pep.nucleotides[:48] + "CCC" + pep.nucleotides[51:]       # the residue stays P, the bases do not
    with pytest.raises(ValueError, match="published nucleotides"):
        catalogue.place(catalogue.Peptide(**{**pep.__dict__, "nucleotides": wrong}), paths.BUILD)
