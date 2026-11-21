"""H18: the bookkeeping of the interface energy (which chain is which protein, what is substituted, what is counted as
not modelled, the tasks, the collection) by hand, and, where OpenMM is installed (`uv sync --group energy`), the cutting
of chains at breaks, the free and restrained atoms and the energy of two residues pulled apart."""
import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.mitoport import cybrid, cybrid_energy as ce

P = "MKVLAWYRTEQNDSGHIFPCMKVLAWYRTEQNDSGHIFPC"          # 40 residues


class TestBookkeeping:
    def test_a_chain_is_matched_to_its_protein(self):
        frag = P[5:35]
        g, m = ce.match_protein(frag, {"X": P, "Y": P[::-1]})
        assert g == "X" and m == {i: i + 5 for i in range(30)}
        mutated = frag[:10] + "AAAA" + frag[14:]                 # 26 of 30 identical: 0.867, under 0.9
        assert ce.match_protein(mutated, {"X": P}) is None
        assert ce.match_protein(frag[:10], {"X": P}) is None      # too short to tell

    def test_substitutions_and_what_is_not_modelled(self):
        host = "MKVLAWYRTE"
        struct = "KVLAWYRT"; to_protein = {i: i + 1 for i in range(8)}     # residues 2 to 9 resolved
        donor = list("MRVL-WYKTD")     # K2R; A5 deleted; R8K; E10D not resolved; M1 same
        subs, unmodelled = ce.substitutions(struct, to_protein, host, donor)
        assert subs == {0: "R", 6: "K"}
        assert unmodelled == 2          # the deletion, and the unresolved E10D
        subs, unmodelled = ce.substitutions(struct, to_protein, host, list("MKVLAWY?TE"))
        assert subs == {} and unmodelled == 1

    def test_the_structure_already_carrying_the_donors_residue_is_left(self):
        subs, _ = ce.substitutions("KR", {0: 0, 1: 1}, "KK", list("KR"))
        assert subs == {}

    def test_insertions(self):
        assert ce.insertions("MKVLAWYRTE", "MKVLAGGWYRTE") == 2
        assert ce.insertions("MKVLAWYRTE", "MKVLAWYRTE") == 0

    def test_tasks_are_the_pairs_with_a_structure(self):
        meas = pd.DataFrame({"host": ["mus-musculus", "homo-sapiens", "homo-sapiens"],
                             "donor": ["rattus-norvegicus", "pan-paniscus", "pan-paniscus"]})
        t = ce.triples(meas)
        assert t == [("homo-sapiens", "pan-paniscus", c) for c in ("I", "III", "IV", "V")] + \
            [("mus-musculus", "rattus-norvegicus", c) for c in ("I", "III", "IV")]
        real = ce.triples(cybrid.read_measurements())
        assert len(real) == len(set(real)) == 5 * 4 + 8 * 3

    def test_collect_says_what_is_missing(self, tmp_path, capsys):
        with pytest.raises(SystemExit):
            ce.collect(tmp_path)
        ce.write({"host": "homo-sapiens", "donor": "pan-paniscus", "complex": "IV", "ddE": 1.5, "status": "computed"}, tmp_path)
        res = ce.collect(tmp_path)
        assert len(res) == 1 and "NOT ASSESSED: 43 of 44" in capsys.readouterr().out

    def test_repeat_build_compares_only_tasks_computed_in_both_folders(self, tmp_path):
        a, b = tmp_path / "a", tmp_path / "b"
        base = {"host": "homo-sapiens", "donor": "pan-paniscus", "substituted": 17, "e_int_reference": -100.0, "status": "computed"}
        with pytest.raises(SystemExit):
            ce.repeat_build(a, b)
        ce.write({**base, "complex": "IV", "ddE": -7.0, "ddE_sd": 3.0}, a)
        ce.write({**base, "complex": "IV", "ddE": 5.0, "ddE_sd": 4.0, "e_int_reference": -90.0}, b)
        ce.write({**base, "complex": "V", "ddE": 88.0, "ddE_sd": 1.0}, a)
        ce.write({**base, "complex": "V", "ddE": 80.5, "ddE_sd": 2.0}, b)
        ce.write({**base, "complex": "III", "ddE": 1.0, "ddE_sd": 1.0}, a)                       # not run a second time
        ce.write({**base, "complex": "I", "ddE": 1.0, "ddE_sd": 1.0}, a)
        ce.write({**base, "complex": "I", "ddE": 0.0, "ddE_sd": 0.0, "status": "no substitution to model"}, b)
        res = ce.repeat_build(a, b).set_index("complex")
        assert sorted(res.index) == ["IV", "V"]
        assert res.loc["IV", "difference"] == 12.0 and res.loc["V", "difference"] == -7.5
        assert res.loc["IV", "sd_of_difference_by_the_starts"] == 5.0
        assert not res.loc["IV", "same_sign"] and res.loc["V", "same_sign"]
        assert (res.loc["IV", "e_int_reference_first"], res.loc["IV", "e_int_reference_second"]) == (-100.0, -90.0)
        assert (b / "h18_repeat_build.tsv").is_file() and not (a / "h18_repeat_build.tsv").exists()
        # one name, two tasks: refused
        ce.write({**base, "complex": "V", "ddE": 80.5, "ddE_sd": 2.0, "substituted": 18}, b)
        with pytest.raises(SystemExit):
            ce.repeat_build(a, b)


try:
    import openmm
except ImportError:
    openmm = None
needs_openmm = pytest.mark.skipif(openmm is None, reason="needs `uv sync --group energy`")


def peptide(seqs, gap_after=None):
    """Chains of extended residues (N, CA, C, O, CB) 0.38 nm apart along x; `gap_after` leaves a hole in chain 0."""
    from openmm import app, unit, Vec3
    top = app.Topology(); pos = []
    el = app.element
    for ci, seq in enumerate(seqs):
        ch = top.addChain(chr(65 + ci)); x = 0.0
        for k, name in enumerate(seq):
            if ci == 0 and gap_after is not None and k == gap_after + 1:
                x += 2.0
            r = top.addResidue(name, ch, str(k + 1))
            for an, e, d in (("N", el.nitrogen, (0, 0, 0)), ("CA", el.carbon, (0.12, 0.08, 0)), ("C", el.carbon, (0.25, 0, 0)),
                             ("O", el.oxygen, (0.25, -0.12, 0)), ("CB", el.carbon, (0.12, 0.2, 0.05))):
                if name == "GLY" and an == "CB":
                    continue
                top.addAtom(an, e, r); pos.append(Vec3(x + d[0], d[1] + ci * 0.6, d[2]))
            x += 0.38
    return top, unit.Quantity(pos, unit.nanometer)


@needs_openmm
class TestStructure:
    def test_a_chain_is_cut_at_a_break_and_nowhere_else(self):
        top, pos = peptide([["ALA"] * 6, ["ALA"] * 3], gap_after=2)
        new, origin = ce.split_at_breaks(top, pos)
        assert [len(list(c.residues())) for c in new.chains()] == [3, 3, 3]
        assert origin == {0: 0, 1: 0, 2: 1}
        assert [r.id for r in new.residues()] == [r.id for r in top.residues()]
        peptide_bonds = [(a.residue.index, b.residue.index) for a, b in new.bonds() if a.residue != b.residue]
        assert sorted(peptide_bonds) == [(0, 1), (1, 2), (3, 4), (4, 5), (6, 7), (7, 8)]

    def test_free_atoms_are_side_chains_near_a_substitution(self):
        top, pos = peptide([["ALA"] * 8, ["ALA"] * 8])
        restrained, free = ce.masks(top, pos, {(0, 0)})
        atoms = list(top.atoms())
        free_names = {(a.residue.chain.index, a.residue.index % 8, a.name) for a, f in zip(atoms, free) if f}
        # within 0.8 nm of residue 1 of chain A: residues 1 to 3 of both chains (0.6 nm apart), side chains only
        assert free_names == {(c, k, "CB") for c in (0, 1) for k in (0, 1, 2)}
        assert (restrained | free).all() and not (restrained & free).any()

    def test_near(self):
        pos = np.array([[0, 0, 0], [0.3, 0, 0], [2, 0, 0]], float)
        assert ce.near(pos, np.array([1, 2]), np.array([0]), 0.45).tolist() == [True, False]
        assert ce.near(pos, np.array([1]), np.array([], int), 0.45).tolist() == [False]


@needs_openmm
def test_on_cox2_of_5z62_the_system_the_substitutions_and_the_energy_of_parted_chains():
    cif, gb = paths.PDB / "5Z62.cif", paths.CYBRID_GB / "records.gb"
    if not cif.is_file() or not gb.is_file():
        pytest.skip("needs `mitoage download pdb 5Z62` and `mitoage download cybrid`")
    from openmm import unit
    rec = cybrid.read_proteins(gb)
    hp, dp = cybrid.species_proteins(rec, "homo-sapiens"), cybrid.species_proteins(rec, "gorilla-gorilla")
    mt = {"COX2": hp["COX2"][0]}
    top, pos, mtc = ce.prepare(cif, mt, {g: p for g, (p, _, _) in hp.items() if g != "COX2"})
    assert list(mtc) == [0] and mtc[0][0] == "COX2" and len(mtc[0][1]) == 227
    seqs = [ce.chain_sequence(c)[1] for c in top.chains()]
    assert seqs[0] == mt["COX2"]
    # COX1 and COX3 touch COX2 but are mtDNA-encoded proteins of the complex not asked for: dropped
    assert hp["COX1"][0] not in seqs and all(ce.match_protein(s, {"COX1": hp["COX1"][0], "COX3": hp["COX3"][0]}) is None for s in seqs)
    subs, unmodelled = ce.substitutions(seqs[0], mtc[0][1], mt["COX2"], cybrid.on_human(mt["COX2"], dp["COX2"][0], False))
    want = {i: d for i, (h, d) in enumerate(zip(mt["COX2"], dp["COX2"][0])) if h != d}      # equal lengths, no indel
    assert subs == want and unmodelled == 0 and len(subs) > 0
    t, p = ce.build(top, pos, {(0, k): aa for k, aa in subs.items()})
    got = ce.chain_sequence(list(t.chains())[0])[1]
    assert got == dp["COX2"][0] and [ce.chain_sequence(c)[1] for c in t.chains()][1:] == seqs[1:]
    # the two groups parted by 30 nm do not interact. In double precision: the structure as built is not relaxed and
    # its hydrogens are placed anew at each build, so its energy can reach 1e8 kJ/mol, where single precision is off
    # by thousands (the check of 2026-11-21: 102 kJ/mol in one build, -3,608 in another, 0.000 in double precision)
    xyz = np.array(p.value_in_unit(unit.nanometer))
    first = np.array([a.residue.chain.index == 0 for a in t.atoms()])
    xyz[first] += np.array([30.0, 0, 0])
    assert abs(ce.interaction(t, unit.Quantity(xyz, unit.nanometer), {0}, platform="Reference")) < 0.01
