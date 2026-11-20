"""H21-H25's code (adaptations/): trait clades, pairs per origin, control draws, the RELAX and aBSREL inputs and
readers, the coevolution statistics and the readout, on the downloaded build, its dated tree and the trait sources.
No test prints an estimate: these run before the pre-registered analysis does."""
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio import Phylo

from mitoage_pipelines import paths
from mitoage_pipelines.adaptations import clades, coupling, fits, pairs, run, selection, traits
from mitoage_pipelines.build2 import pins

FB = Path(os.environ.get("MITOAGE_FISHBASE", paths.B2_FISHBASE / "v26.06"))
NAMES = Path(os.environ.get("MITOAGE_TAXDUMP", paths.B2_TAXDUMP / pins.TAXDUMP_BUILD2)) / "names.dmp"
DATA = Path(__file__).parent / "data"
need = pytest.mark.skipif(not all(p.is_file() for p in (paths.BUILD / "species.parquet", paths.DATED_TREE, FB / "species.parquet")),
                          reason="needs `download build`, `longevity dated-tree` and `download fishbase-temperature`")


@pytest.fixture(scope="module")
def species():
    sp = pd.read_parquet(paths.BUILD / "species.parquet")
    return sp[sp.record_status.ne("retired")].reset_index(drop=True)


@pytest.fixture(scope="module")
def tr(species):
    return traits.classify(species, FB, paths.HIBERNATION, NAMES if NAMES.is_file() else None)


@pytest.fixture(scope="module")
def tree():
    return clades.read(paths.DATED_TREE)


def known_tree(tree, tr, trait):
    t = tr[tr.trait.eq(trait) & tr.state.notna()]
    tips = set(clades.tips(tree))
    state = {s: int(v) for s, v in zip(t.slug, t.state) if s in tips}
    root = clades.prune(tree, set(state))
    return root, {s: state[s] for s in clades.tips(root)}


@need
def test_traits_follow_their_sources(tr, species):
    fb = pd.read_parquet(FB / "species.parquet", columns=["SpecCode", "AnaCat", "DemersPelag"])
    assert tr.groupby("trait").slug.apply(lambda s: s.is_unique).all()
    d = tr[tr.trait.eq("diadromy")]
    assert set(d["class"]) == {"Actinopterygii"}
    assert (d.state.eq(1) == d.how.str.removeprefix("AnaCat ").isin(traits.DIADROMOUS)).all()
    assert d.state[d.how.str.startswith("AnaCat empty") | d.how.str.startswith("no FishBase")].isna().all()
    dd = tr[tr.trait.eq("deep_depth")]
    assert dd.state[dd.how.str.contains("intermediate")].isna().all()
    b = tr[tr.trait.eq("bats")]
    assert b.state.sum() == species[species["class"].eq("Mammalia") & species.order.eq("Chiroptera")].shape[0]


@need
def test_prune_keeps_distances(tree):
    names = sorted(clades.tips(tree))[::40]
    p = clades.prune(tree, set(names))
    assert sorted(clades.tips(p)) == names
    ph = Phylo.read(paths.DATED_TREE, "newick")
    a, b = names[0], names[-1]
    node = {n.name: n for n in clades.walk(p) if n.is_tip}
    up = []
    n = node[a]
    while n is not None: up.append(n); n = n.parent
    n, d = node[b], 0.0
    while n not in up: d += n.length; n = n.parent
    d += clades.depth_to(node[a], n)
    assert d == pytest.approx(ph.distance(a, b), rel=1e-6)


@need
@pytest.mark.parametrize("trait", ["diadromy", "hibernation", "deep_habitat"])
def test_trait_clades_are_maximal_and_cover(tree, tr, trait):
    root, state = known_tree(tree, tr, trait)
    cl = clades.trait_clades(root, state)
    covered = [s for c in cl for s in clades.tips(c)]
    assert sorted(covered) == sorted(s for s, v in state.items() if v == 1)
    for c in cl:
        assert all(state[s] == 1 for s in clades.tips(c))
        assert c.parent is None or any(state[s] == 0 for s in clades.tips(c.parent))


@need
@pytest.mark.parametrize("trait", ["diadromy", "deep_habitat"])
def test_origin_pairs(tree, tr, trait):
    root, state = known_tree(tree, tr, trait)
    rows = clades.origin_pairs(root, state)
    cl = {clades.label(c): set(clades.tips(c)) for c in clades.trait_clades(root, state)}
    assert len(rows) == len(cl)
    node = {n.name: n for n in clades.walk(root) if n.is_tip}
    def edges(a, b):
        A, B = node[a], node[b]; up = set()
        n = A
        while n is not None: up.add(id(n)); n = n.parent
        n = B
        while id(n) not in up: n = n.parent
        return clades.path_edges(A, B, n), n
    taken = [r for r in rows if r["taken"]]
    used = []
    for r in taken:
        assert r["sp_trait"] in cl[r["clade"]] and r["sp_other"] not in cl[r["clade"]]
        assert state[r["sp_other"]] == 0
        e, anc = edges(r["sp_trait"], r["sp_other"])
        assert r["outgroup"] not in clades.tips(anc)
        assert all(not (e & u) for u in used)
        used.append(e)
    for r in rows:
        if not r["taken"] and r["reason"].startswith("path shares"):
            e, _ = edges(r["sp_trait"], r["sp_other"])
            assert any(e & u for u in used)


@need
def test_control_draws(tree, tr):
    root, state = known_tree(tree, tr, "diadromy")
    cl = clades.trait_clades(root, state)
    draws, short = clades.control_matches(root, state, cl, 5)
    again, _ = clades.control_matches(root, state, cl, 5)
    lab = lambda d: [(clades.label(c), clades.label(z)) for c, z in d]
    assert [lab(d) for d in draws] == [lab(d) for d in again]
    assert len({tuple(sorted(lab(d))) for d in draws}) > 1
    zeros = clades.zero_clades(root, state)
    for d, s in zip(draws, short):
        assert len(d) + s == len(cl)
        seen = set()
        for c, z in d:
            t, n = set(clades.tips(z)), len(clades.tips(c))
            assert all(state[x] == 0 for x in t) and not (t & seen)
            seen |= t
            assert n / 2 <= len(t) <= 2 * n
            # among the size-matched candidates, the drawn one is within the POOL nearest in stem age
            near = sorted((abs(clades.stem_age(y) - clades.stem_age(c)), clades.label(y)) for y in zeros
                          if n / 2 <= len(clades.tips(y)) <= 2 * n)
            assert (abs(clades.stem_age(z) - clades.stem_age(c)), clades.label(z)) in near[:clades.POOL]
    assert [[z for _, z in d] for d in draws] == clades.control_draws(root, state, cl, 5)[0]


@need
def test_selection_inputs(tree, tr, species, tmp_path):
    root, state = known_tree(tree, tr, "bat_hibernation")
    cl = clades.trait_clades(root, state)
    keep = selection.draw_tips(root, state, cl)
    inside = {s for c in cl for s in clades.tips(c)}
    assert inside <= keep and all(state[s] == 0 for s in keep - inside)
    assert len(keep - inside) <= selection.N_REFERENCE * len(cl)
    # in a control draw the clades are all-0, so the trait species are outside them: none may be a reference tip
    ctl = clades.control_draws(root, state, cl, 3)[0]
    for d in ctl:
        k2 = selection.draw_tips(root, state, d)
        in2 = {s for c in d for s in clades.tips(c)}
        assert all(state[s] == 0 for s in k2 - in2)
    t, nwk, n_test = selection.labelled(root, cl, keep, stems_only=False)
    expect = sum(1 for n in clades.walk(t) if set(clades.tips(n)) <= inside and any(set(clades.tips(n)) <= set(clades.tips(c)) for c in cl))
    assert nwk.count("{Test}") == n_test == expect
    assert nwk.count("{Reference}") == sum(1 for _ in clades.walk(t)) - n_test
    _, anwk, n_stem = selection.labelled(root, cl, keep, stems_only=True)
    assert anwk.count("{Test}") == n_stem == len(cl) and "{Reference}" not in anwk
    info = selection.write_codons(paths.BUILD, sorted(keep), tmp_path / "c.fna")
    seqs = [l for l in (tmp_path / "c.fna").read_text().split("\n") if l and not l.startswith(">")]
    assert len(seqs) == len(keep) and {len(s) for s in seqs} == {3 * info["codons"]}
    stops = {"TAA", "TAG", "AGA", "AGG"}
    assert not any(s[i:i + 3] in stops for s in seqs for i in range(0, len(s), 3))


@need
def test_planted_stop_codon_is_masked(tree, tr, monkeypatch, tmp_path):
    """A stop codon planted in a copy of a real alignment is written as a gap, and counted."""
    from mitoage_pipelines.seqio import read_fasta as real_read
    tips = sorted(clades.tips(tree))[:6]
    def planted(path):
        a = real_read(path)
        if "COX1" in str(path):
            s = a[tips[0]]
            i = next(j for j in range(0, len(s), 3) if s[j:j + 3] not in ("---",))
            a[tips[0]] = s[:i] + "TAA" + s[i + 3:]
        return a
    monkeypatch.setattr(selection, "read_fasta", planted)
    info = selection.write_codons(paths.BUILD, tips, tmp_path / "c.fna")
    seq = (tmp_path / "c.fna").read_text().split("\n")[1]
    assert info["masked"] >= 1 and not any(seq[i:i + 3] == "TAA" for i in range(0, len(seq), 3))


def test_read_relax_real_output(tmp_path):
    r = selection.read_relax(DATA / "relax_ATP8_h3.json")
    assert r == pytest.approx({"LRT": 1.674722064417438, "p": 0.1956272448343103, "K": 1.896716246758435})
    bad = json.loads((DATA / "relax_ATP8_h3.json").read_text()); del bad["test results"]
    (tmp_path / "b.json").write_text(json.dumps(bad))
    with pytest.raises(selection.Refused):
        selection.read_relax(tmp_path / "b.json")


def test_summarise_reads_the_draw_folders_beside_the_draws_table(tmp_path):
    """The folder `selection-prep` writes holds draws.tsv beside draw00..draw20: only the folders are draws."""
    (tmp_path / "draws.tsv").write_text("trait\tdraw\n")
    for k in (0, 3):
        (tmp_path / f"draw{k:02d}").mkdir()
    (tmp_path / "draw00" / "relax.json").write_bytes((DATA / "relax_ATP8_h3.json").read_bytes())
    res = selection.summarise(tmp_path)
    assert res.draw.tolist() == [0, 3]
    assert res.status.tolist()[1] == "missing" and res.K[0] == selection.read_relax(DATA / "relax_ATP8_h3.json")["K"]


@pytest.mark.skipif(not (DATA / "absrel_nd3_fixture.json").is_file(), reason="the aBSREL fixture is not in tests/data")
def test_read_absrel_real_output(tmp_path):
    t = selection.read_absrel(DATA / "absrel_nd3_fixture.json")
    raw = json.loads((DATA / "absrel_nd3_fixture.json").read_text())
    assert set(t.branch) == {b for b, v in raw["tested"]["0"].items() if v == "test"}
    assert ((t.p_corrected >= t.p_uncorrected - 1e-12) & t.p_corrected.between(0, 1)).all()
    raw["tested"]["0"] = {b: "background" for b in raw["tested"]["0"]}
    (tmp_path / "b.json").write_text(json.dumps(raw))
    with pytest.raises(selection.Refused):
        selection.read_absrel(tmp_path / "b.json")


def test_r1_needs_every_control():
    res = pd.DataFrame({"draw": range(21), "LRT": [9.0] + list(range(20)), "p": 0.01, "K": 1.3})
    r = selection.r1(res)
    assert r["status"] == "assessed" and r["controls_below"] == 9      # LRTs 0..8 of 0..19 are below 9.0
    assert selection.r1(res.iloc[:15])["status"].startswith("not assessed")


@need
def test_pair_rates(species, tr, tree):
    dist = pairs.Distances(paths.BUILD, species)
    p = pairs.run(paths.DATED_TREE, tr, "diadromy", set(clades.tips(tree)) & dist.aligned, dist)
    ok = p[p.usable.fillna(False).astype(bool)]
    assert len(ok) > 0
    assert np.allclose(ok.d_omega, ok.d_dN - ok.d_dS)
    assert np.allclose(ok.d_dS, np.log10(ok.dS_trait) - np.log10(ok.dS_other))
    assert (ok.saturated == (ok[["dS_trait", "dS_other"]].max(axis=1) > 1.0)).all()
    r = ok.iloc[0]
    dt, do = dist.to_outgroup(r.sp_trait, r.sp_other, r.outgroup)
    assert dt["dS"] == pytest.approx(r.dS_trait)
    assert do["dS"] == pytest.approx(r.dS_other)
    assert dt["codons"] == do["codons"]
    assert dist(r.outgroup, r.sp_trait) == dist(r.sp_trait, r.outgroup)


def test_pair_test_minimum_and_saturation():
    p = pd.DataFrame({"usable": True, "saturated": [False] * 7 + [True] * 5, "d_omega": np.linspace(-1, 1, 12)})
    assert pairs.test(p, "d_omega")["status"].startswith("not assessed")
    r = pairs.test(p, "d_omega", with_saturated=True)
    assert r["n_pairs"] == 12 and r["status"] == "assessed" and r["ci_low"] <= r["median"] <= r["ci_high"]


def _bl_on(topology: clades.Node, genes, fg_keys, shift, rng):
    allt = frozenset(clades.tips(topology))
    keys = []
    for n in clades.walk(topology):
        if n.parent is None: continue
        s = frozenset(clades.tips(n)); keys.append(min(s, allt - s, key=lambda z: (len(z), sorted(z))))
    out = {}
    for g in genes:
        base = rng.uniform(0.5, 1.5, len(keys))
        out[g] = {k: b * (shift.get(g, 1.0) if k in fg_keys else 1.0) for k, b in zip(keys, base)}
    return out


@need
def test_coupling_statistics_find_a_planted_shift(tree, tr):
    mam = set(tr.slug[tr.trait.eq("bats")])
    sub = sorted(mam & set(clades.tips(tree)))[::3]
    full = clades.prune(tree, set(sub))
    state = {s: int(v) for s, v in zip(tr.slug[tr.trait.eq("bats")], tr.state[tr.trait.eq("bats")]) if s in sub}
    cl = clades.trait_clades(full, state)
    genes = [coupling.mt_name(g) for g in ("COX1", "CYTB", "ND1")] + coupling.SUBUNITS[:12] + coupling.CYTOSOLIC_CONTROLS[:12]
    rng = np.random.default_rng(1)
    allt = frozenset(clades.tips(full))
    fgk = set()
    for c in cl:
        for n in clades.walk(c):
            s = frozenset(clades.tips(n)); fgk.add(min(s, allt - s, key=lambda z: (len(z), sorted(z))))
    planted = _bl_on(full, genes, fgk, {g: 3.0 for g in coupling.SUBUNITS[:12]}, rng)
    logr, logm, zeros = coupling.relative_rates(planted, genes)
    fg = coupling.foreground(logr.attrs["keys"], full, cl, set(state))
    assert fg.sum() == len(fgk & set(logr.attrs["keys"]))
    r = coupling.c1(logr, fg, "C1")
    assert r["status"] == "assessed" and r["estimate"] > 0.5 and r["p"] < 0.01
    flat = _bl_on(full, genes, fgk, {}, np.random.default_rng(2))
    logr2, _, _ = coupling.relative_rates(flat, genes)
    assert abs(coupling.c1(logr2, fg, "C1")["estimate"]) < 0.5
    stem = coupling.foreground(logr.attrs["keys"], full, cl, set(state), stems_only=True)
    assert stem.sum() == len(cl)
    # a control protein's background leaves its own length out
    L = pd.DataFrame({g: [planted[g][k] for k in logr.attrs["keys"]] for g in genes})
    sc = L / L.mean(); ctl = [g for g in genes if g in coupling.CYTOSOLIC_CONTROLS]
    g0 = ctl[0]; bg = sc[[c for c in ctl if c != g0]].mean(axis=1)
    assert np.allclose(logr[g0], np.log(sc[g0] / bg))


def test_l2_reading():
    feats = [f"f{i}" for i in range(8)]
    def tab(slopes, q):
        return pd.DataFrame({"status": "assessed", "term": "ll", "feature": feats, "n": 100, "slope_sd": slopes, "q": q})
    a = tab([1.0] * 8, [0.01] * 8)
    assert fits.l2(a, tab([0.3] * 8, [0.01] * 2 + [0.5] * 6))["reading"] == "largely hibernation's"
    assert fits.l2(a, tab([0.9] * 8, [0.01] * 8))["reading"] == "not hibernation's"
    assert fits.l2(a, tab([0.6] * 8, [0.01] * 8))["reading"] == "partly"
    assert fits.l2(a, tab([0.3] * 8, [0.01] * 8))["reading"] == "partly"      # small slopes, count not halved
    assert fits.l2(tab([1.0] * 8, [0.5] * 8), a)["status"].startswith("not assessed")
    with pytest.raises(ValueError):
        fits.l2(a, tab([0.3] * 8, [0.01] * 8).assign(n=99))


def test_readout_readings():
    R = lambda **k: type("R", (), {"status": "assessed", "kind": "test", "primary": True, "p_holm": np.nan, "detail": {}, **k})
    assert run.reading(R(readout="R1", p=0.01, p_holm=0.01, estimate=1.4, predicted=1, detail={"controls_below": 19}), True) == "intensified: as predicted"
    assert run.reading(R(readout="R1", p=0.01, p_holm=0.01, estimate=0.7, predicted=1, detail={"controls_below": 20}), True) == "relaxed: opposite to the prediction"
    assert run.reading(R(readout="R1", p=0.01, p_holm=0.01, estimate=1.4, predicted=1, detail={"controls_below": 18}), True) == "not separated from the control clades"
    assert run.reading(R(readout="R1", p=0.01, p_holm=0.08, estimate=1.4, predicted=1, detail={"controls_below": 20}), True) == "not detected"
    assert run.reading(R(readout="R2", p=0.001, p_holm=0.004, estimate=-0.1, predicted=-1), True) == "negative: as predicted"
    assert run.reading(R(readout="C1", p=0.001, p_holm=0.004, estimate=0.2, predicted=0), False).startswith("not read")
    assert run.reading(R(readout="P", kind="positive control", p=0.2, estimate=0.3, predicted=1), True) == "fails"
