"""H14: the contact maps of complexes I, III, IV and V, and the readouts fixed in advance (`mitoport readouts`).

The structures, the build's mitochondrial alignments and the dated tree are real; expected values are computed in the
test from those files. Where nuclear ortholog alignments would carry a result, the tests build alignments in which every
species carries the human residue: then a pair is discordant exactly when its mitochondrial residue is off its
consensus, which the test counts on its own."""
import dataclasses
import hashlib
import shutil
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy import stats
from statsmodels.stats.multitest import multipletests

from mitoage_pipelines import paths
from mitoage_pipelines.genes import CYTOSOLIC_CONTROLS, MT_GENES
from mitoage_pipelines.mitoport import burial, complexes, contacts, iface, partners, readouts
from mitoage_pipelines.seqio import read_fasta

REPO = Path(__file__).resolve().parents[1]
CIFS = {n: paths.PDB / f"{cx.pdb}.cif" for n, cx in complexes.COMPLEXES.items()}
needs_pdb = pytest.mark.skipif(not all(p.is_file() for p in CIFS.values()),
                               reason="needs `download pdb 5XTD 5Z62 8H9S 5XTE`")
needs_build = pytest.mark.skipif(not (paths.BUILD / "alignments" / "ATP8.aa.fasta.gz").is_file() or not paths.DATED_TREE.is_file()
                                 or not paths.QUALITY_OK.is_file(), reason="needs the build and `longevity dated-tree`")

# complex IV's contact map as the complex-IV-only code wrote it (commit 7692a55) from build c525adb and the site's tree at
# 6fe8e72; sites.tsv's digest is the one H12's jobs printed on haumea
IV_DIGESTS = {"chains.tsv": "853d29226eccebf977fc6d9dd6973107b7e1e515b5917a18f5c1753a34a44df5",
              "contacts.tsv": "dc712a025912e5d7fb418d94a88bfa23016533c74729849b54c4c8dcbc17395c",
              "control_pairs.tsv": "e6c1923d2d051f24ac302535bcd003c90cc02f7f1e8696569aa6811b59f80298",
              "sites.tsv": "6f219d3a7d70e10fce1a5fdf8220b46eb813291eba1159376c17f7b1d0f5e4c7",
              "species.txt": "54706d4569b7e4bac1029d507e6cd182e7ef5af565aa61f9b57c8abfe76ec37a"}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(work: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(work / name, sep="\t", dtype={"chain": str, "partner_chain": str})


# --- the table of entries --------------------------------------------------------------------------
@needs_pdb
class TestEntries:
    def test_complex_iv_names_agree_with_the_descriptions_h5_read(self):
        """5Z62's chains by UniProt reference are the genes H5's description rules give."""
        d = contacts.entity_descriptions(CIFS["IV"])
        got = complexes.chain_genes(CIFS["IV"], complexes.get("IV"))
        assert got == {c: contacts.mt_gene(v) or burial.gene_of(v) for c, v in d.items()}
        assert len(got) == 14

    @pytest.mark.parametrize("name", list(complexes.COMPLEXES))
    def test_every_chain_is_named_and_each_mt_gene_has_its_reference(self, name):
        cx = complexes.get(name)
        got = complexes.chain_genes(CIFS[name], cx)
        model = MMCIFParser(QUIET=True).get_structure("x", str(CIFS[name]))[0]
        assert set(got) == {c.id for c in model}                                 # every chain of the model
        assert {g for g in got.values() if g in MT_GENES} == set(cx.mt)
        assert set(got.values()) - set(MT_GENES) == set(cx.nuclear)
        for g, c in cx.mt.items():
            assert got[c] == g

    def test_the_second_cytochrome_b_of_the_dimer_is_not_the_reference(self):
        got = complexes.chain_genes(CIFS["III"], complexes.get("III"))
        assert sorted(c for c, g in got.items() if g == "CYTB") == ["J", "V"]
        assert complexes.get("III").mt == {"CYTB": "J"}
        assert sum(g == "UQCRFS1" for g in got.values()) == 4                  # presequence and mature, twice

    def test_an_entity_the_table_does_not_name_is_refused(self):
        cx = complexes.get("V")
        short = dataclasses.replace(cx, uniprot={k: v for k, v in cx.uniprot.items() if v != "ATP5IF1"})
        with pytest.raises(SystemExit, match="refusing to guess"):
            complexes.chain_genes(CIFS["V"], short)

    def test_a_reference_chain_of_another_gene_is_refused(self):
        cx = complexes.get("I")
        wrong = dataclasses.replace(cx, mt={**cx.mt, "ND1": "i"})              # chain i is ND2
        with pytest.raises(SystemExit, match="not ND1"):
            complexes.chain_genes(CIFS["I"], wrong)


# --- the contact maps, made here from the structures ----------------------------------------------
@pytest.fixture(scope="module")
def maps(tmp_path_factory):
    out = {}
    for name in ("III", "IV", "V"):
        cx = complexes.get(name)
        w = tmp_path_factory.mktemp(f"map{name}")
        contacts.main(CIFS[name], paths.BUILD, w, cx)
        burial.main(CIFS[name], w / "contacts.tsv", paths.BUILD, paths.DATED_TREE, paths.QUALITY_OK, w, cx)
        out[name] = w
    return out


@needs_pdb
@needs_build
class TestContactMaps:
    def test_complex_iv_is_byte_for_byte_the_published_map(self, maps):
        assert {f: sha(maps["IV"] / f) for f in IV_DIGESTS} == IV_DIGESTS

    def test_the_other_copy_of_cytochrome_b_is_a_mitochondrial_partner(self, maps):
        ct, ch = read(maps["III"], "contacts.tsv"), read(maps["III"], "chains.tsv")
        assert set(ct.chain) == {"J"}
        other = ct[ct.partner_chain.eq("V")]
        assert len(other) > 0 and other.partner_is_mt.eq(1).all()
        assert dict(zip(ch.chain, ch.is_mt))["V"] == 1
        mt_chains = set(ch.chain[ch.gene.isin(MT_GENES)])
        assert set(ct.partner_chain[ct.partner_is_mt.eq(1)]) <= mt_chains
        assert not set(ct.partner_chain[ct.partner_is_mt.eq(0)]) & mt_chains

    def test_contacts_count_nuclear_partners_only(self, maps):
        s, ct = read(maps["III"], "sites.tsv"), read(maps["III"], "contacts.tsv")
        nuc = ct[ct.partner_is_mt.eq(0)].groupby("pdb_resnum").size()
        got = s.set_index("pdb_resnum").n_contacts
        assert (got[got > 0].sort_index() == nuc.sort_index()).all() and set(got[got > 0].index) == set(nuc.index)
        assert set(s.pdb_resnum[s.group.eq("contact")]) == set(nuc.index)

    def test_species_hold_every_mitochondrial_alignment_of_the_complex(self, maps):
        cx = complexes.get("V")
        alns = {g: set(read_fasta(paths.BUILD / "alignments" / f"{g}.aa.fasta.gz")) for g in cx.mt}
        want = sorted(set.intersection(*alns.values()) & burial.tree_tips(paths.DATED_TREE.read_text())
                      & set(paths.QUALITY_OK.read_text().split()))
        assert (maps["V"] / "species.txt").read_text().split() == want
        assert alns["ATP6"] != alns["ATP8"]                       # the two alignments do differ: the intersection matters

    def test_a_gene_that_runs_out_of_controls_keeps_its_contacts(self, maps, capsys):
        s, m = read(maps["V"], "sites.tsv"), read(maps["V"], "control_pairs.tsv")
        a8 = s[s.gene.eq("ATP8")]
        con, ctl = a8.group.eq("contact").sum(), a8.group.eq("control").sum()
        assert con > ctl > 0 and a8.group.eq("other").sum() == 0                # every non-contact became a control
        assert len(m[m.gene.eq("ATP8")]) == ctl
        got = partners.control_of(s, m)
        assert len(got) == len(m) == s.group.eq("control").sum()
        assert f"{con - ctl} of {con} contact residues have no control" in capsys.readouterr().out
        assert {(g, c) for (g, _), c in got.items()} == set(zip(s.gene[s.group.eq("control")], s.aln_col[s.group.eq("control")]))


def entity_rows(cif: Path, genes: dict[str, str]) -> dict[str, str]:
    """{gene: its entities' canonical sequences joined, in entity order}: a stand-in for the human ortholog."""
    d = MMCIF2Dict(str(cif))
    seq = {}
    for ent, chains, s in zip(d["_entity_poly.entity_id"], d["_entity_poly.pdbx_strand_id"], d["_entity_poly.pdbx_seq_one_letter_code_can"]):
        g = genes[chains.split(",")[0].strip()]
        seq[g] = seq.get(g, "") + s.replace("\n", "")
    return seq


@pytest.fixture(scope="module")
def synthetic(maps, tmp_path_factory):
    """Nuclear alignments of complexes III and V in which 40 real species all carry the human residue."""
    out = {}
    for name in ("III", "V"):
        cx = complexes.get(name)
        species = (maps[name] / "species.txt").read_text().split()[:: 50][:40]
        rows = entity_rows(CIFS[name], complexes.chain_genes(CIFS[name], cx))
        folder = tmp_path_factory.mktemp(f"aln{name}")
        for g in cx.nuclear:
            (folder / f"{g}.afa").write_text("".join(f">{s}\n{rows[g]}\n" for s in ["homo-sapiens", *species]))
        out[name] = folder, species
    return out


@needs_pdb
@needs_build
class TestNuclearSide:
    @pytest.mark.parametrize("name", ["III", "V"])
    def test_a_gene_in_several_chains_takes_the_union_of_their_interfaces(self, maps, synthetic, name):
        folder, _ = synthetic[name]
        ch, ct = read(maps[name], "chains.tsv"), read(maps[name], "contacts.tsv")
        ct = ct[ct.partner_is_mt.eq(0)]
        model = MMCIFParser(QUIET=True).get_structure("x", str(CIFS[name]))[0]
        alns, full, other = iface.interface_columns(model, ch, ct, folder)
        multi = [g for g, n in ch[ch.is_mt.eq(0)].gene.value_counts().items() if n > 1]
        assert multi
        for g in multi:
            parts = [iface.interface_columns(model, ch[ch.chain.eq(c)], ct, folder)[1][g] for c in ch.chain[ch.gene.eq(g)]]
            assert set(full[g]) == set().union(*map(set, parts))
            assert set(full[g]).isdisjoint(other[g])
            assert set(full[g]) | set(other[g]) == set(range(len(alns[g]["homo-sapiens"])))
        for g in ch.gene[ch.is_mt.eq(0)].unique():                              # one chain: one column per touching residue
            cids = list(ch.chain[ch.gene.eq(g)])
            if len(cids) == 1:
                assert len(full[g]) == ct[ct.partner_chain.eq(cids[0])].partner_resnum.nunique()

    def test_the_usual_size_leaves_out_subunits_that_touch_nothing(self):
        cols = {"a": np.arange(10), "b": np.array([], dtype=int), "c": np.arange(20)}
        assert iface.usual_size(cols) == 15
        assert iface.usual_size({"b": np.array([], dtype=int)}) == 30

    @pytest.mark.parametrize("name", ["III", "V"])
    def test_with_identical_partners_a_pair_is_discordant_when_its_mito_residue_is_off(self, maps, synthetic, name):
        cx = complexes.get(name)
        folder, species = synthetic[name]
        nuc = {g: read_fasta(folder / f"{g}.afa") for g in cx.nuclear}
        mt = {g: read_fasta(paths.BUILD / "alignments" / f"{g}.aa.fasta.gz") for g in cx.mt}
        s = read(maps[name], "sites.tsv")
        pr = partners.pairs(read(maps[name], "contacts.tsv"), s,
                            partners.nuclear_columns(CIFS[name], read(maps[name], "chains.tsv"), nuc),
                            partners.control_of(s, read(maps[name], "control_pairs.tsv")))
        assert set(pr.nuc_gene) <= set(cx.nuclear) and set(pr.mt_gene) <= set(cx.mt)
        got = partners.discordance(pr, mt, nuc, species)

        def mode(g, col):
            return Counter(sorted(mt[g][x][col] for x in species if mt[g][x][col] != "-")).most_common(1)[0][0]
        for sp in species:
            for which, col in (("contact", "mt_col"), ("control", "ctl_col")):
                filled = [(r.mt_gene, getattr(r, col)) for r in pr.itertuples() if mt[r.mt_gene][sp][getattr(r, col)] != "-"]
                off = sum(mt[g][sp][c] != mode(g, c) for g, c in filled)
                assert got.loc[sp, f"{which}_pairs"] == len(filled)
                if got.loc[sp, "contact_pairs"] >= partners.MIN_SHARE * len(pr):
                    assert got.loc[sp, which] == pytest.approx(off / len(filled), abs=1e-12)


@needs_pdb
@needs_build
class TestCoupling:
    def test_a_contact_without_a_column_is_refused(self, maps, tmp_path):
        work = tmp_path / "work"; shutil.copytree(maps["IV"], work)
        s = read(work, "sites.tsv")
        s.loc[s.index[s.group.eq("contact")][0], "aln_col"] = 0
        s.to_csv(work / "sites.tsv", sep="\t", index=False)
        with pytest.raises(SystemExit, match="no alignment column"):
            iface.main(work, tmp_path, work / "contacts.tsv", CIFS["IV"], paths.BUILD, None, paths.DATED_TREE,
                       tmp_path / "out", skip_esm=True)

    @pytest.mark.skipif(not (paths.MITOPORT_IFACE_ALN / "COX4I1.afa").is_file(), reason="needs `mitoport iface-align`")
    def test_alignments_of_other_complexes_are_no_cytosolic_controls(self, maps, tmp_path):
        """Complex IV on H5's 62 alignments and on the same folder with one more (a complex I subunit) give the same
        bytes: only the listed cytosolic proteins are controls."""
        from mitoage_pipelines.genes import NUCLEAR
        h5, more = tmp_path / "h5", tmp_path / "more"
        h5.mkdir(); more.mkdir()
        for g in NUCLEAR:
            if (paths.MITOPORT_IFACE_ALN / f"{g}.afa").is_file():
                for d in (h5, more):
                    (d / f"{g}.afa").symlink_to(paths.MITOPORT_IFACE_ALN / f"{g}.afa")
        (more / "NDUFS2.afa").write_text((paths.MITOPORT_IFACE_ALN / "PGK1.afa").read_text())
        for d in (h5, more):
            iface.main(maps["IV"], d, maps["IV"] / "contacts.tsv", CIFS["IV"], paths.BUILD, None, paths.DATED_TREE,
                       tmp_path / f"out-{d.name}", skip_esm=True)
        for f in ("coupling.tsv", "lifespan.tsv"):
            assert (tmp_path / "out-h5" / f).read_bytes() == (tmp_path / "out-more" / f).read_bytes()
        cp = pd.read_csv(tmp_path / "out-more" / "coupling.tsv", sep="\t")
        assert set(cp.gene[cp.group.eq("cytosolic control")]) <= set(CYTOSOLIC_CONTROLS)
        assert "NDUFS2" not in set(cp.gene)


# --- the readouts --------------------------------------------------------------------------------
class TestReadouts:
    def test_the_intercept_on_independent_tips_is_the_mean(self):
        rng = np.random.default_rng(3)
        y = rng.normal(0.02, 0.05, 12)
        r = readouts.intercept(np.eye(12), y)
        se = y.std(ddof=1) / np.sqrt(12)
        assert r["estimate"] == pytest.approx(y.mean())
        h = stats.t.ppf(0.975, 11) * se
        assert (r["ci95_low"], r["ci95_high"]) == pytest.approx((y.mean() - h, y.mean() + h))
        assert r["p"] == pytest.approx(stats.ttest_1samp(y, 0).pvalue)

    def test_shared_history_widens_the_interval(self):
        """Three clades of ten tips that share most of their history, with values set mostly by the clade: lambda comes
        out high and the interval of the mean is wider than if the 30 tips were independent."""
        C = np.kron(np.eye(3), np.full((10, 10), 0.9)); np.fill_diagonal(C, 1.0)
        y = np.repeat([2.0, -1.0, 0.5], 10) + np.random.default_rng(4).normal(0, 0.1, 30)
        star, tree = readouts.intercept(np.eye(30), y), readouts.intercept(C, y)
        assert 0.5 < tree["lambda"] <= 1
        assert tree["ci95_high"] - tree["ci95_low"] > 1.2 * (star["ci95_high"] - star["ci95_low"]) and tree["p"] > star["p"]
        assert tree["estimate"] == pytest.approx(y.mean())        # equal clades: the GLS mean is the plain mean

    def test_group_test_is_two_sided_mann_whitney_and_drops_unknowns(self):
        a = pd.Series([0.3, 0.1, np.nan, 0.5, 0.2]); b = pd.Series([0.0, -0.1, 0.05, 0.02, np.nan, 0.01])
        r = readouts.group_test({}, a, b, "positive")
        assert (r["n_a"], r["n_b"]) == (4, 5)
        assert r["p"] == pytest.approx(stats.mannwhitneyu(a.dropna(), b.dropna(), alternative="two-sided").pvalue)
        assert r["estimate"] == pytest.approx(a.median() - b.median())

    def test_too_few_proteins_is_not_assessed_not_negative(self):
        r = readouts.group_test({}, pd.Series([0.1, np.nan]), pd.Series([0.0, 0.1, 0.2]), "positive")
        assert r["status"].startswith("not assessed") and "p" not in r

    def test_r1b_and_r2_take_subunits_with_an_interface_r0_takes_all(self):
        cx = complexes.get("III")
        rng = np.random.default_rng(5)
        part = pd.DataFrame({"group": cx.partner_group, "gene": cx.nuclear, "iface_cols": [0, 2, 9, 12, 30, 7, 4, 8, 20]})
        ctl = pd.DataFrame({"group": "cytosolic control", "gene": CYTOSOLIC_CONTROLS, "iface_cols": 25})
        cp = pd.concat([part, ctl], ignore_index=True)
        for c in ("seq_control_other", "seq_contact_iface", "seq_excess_mt", "seq_excess_nuc"):
            cp[c] = rng.normal(size=len(cp))
        rows = {r["measure"].split(":")[0]: r for r in readouts.coupling(cx, cp)}
        touch = part.iface_cols.ge(iface.MIN_IFACE).sum()
        assert rows["seq_control_other"]["n_a"] == len(cx.nuclear) and rows["seq_excess_mt"]["n_a"] == touch == 6
        assert rows["seq_excess_nuc"]["n_a"] == touch and rows["seq_excess_mt"]["n_b"] == len(CYTOSOLIC_CONTROLS)

    def test_holm_covers_r1b_r2_and_mammalian_r3_only(self):
        res = pd.DataFrame([
            {"readout": "R1b", "set": "mammals", "status": "tested", "p": 0.01},
            {"readout": "R2", "set": "mammals", "status": "tested", "p": 0.04},
            {"readout": "R3", "set": "Mammalia", "status": "tested", "p": 0.002},
            {"readout": "R3", "set": "Aves", "status": "tested", "p": 0.001},
            {"readout": "R0", "set": "mammals", "status": "tested", "p": 1e-6},
            {"readout": "R2", "set": "mammals", "status": "not assessed: 2 and 30 proteins with a value", "p": np.nan}])
        got = readouts.holm(res)
        want = multipletests([0.01, 0.04, 0.002], method="holm")[1]
        assert got[:3].tolist() == pytest.approx(want.tolist())
        assert got[3:].isna().all()

    def test_a_class_without_its_table_is_not_assessed(self):
        rows = readouts.concordance(complexes.get("I"), None, paths.DATED_TREE)
        assert len(rows) == 4 and all(r["status"].startswith("not assessed") for r in rows)

    @needs_pdb
    @needs_build
    def test_conservation_pairs_each_contact_with_its_own_control(self, maps):
        s, m = read(maps["IV"], "sites.tsv"), read(maps["IV"], "control_pairs.tsv")
        group, paired = readouts.conservation(complexes.get("IV"), s, m)
        con, ctl = s.identity_to_human[s.group.eq("contact")], s.identity_to_human[s.group.eq("control")]
        assert group["p"] == pytest.approx(stats.mannwhitneyu(con, ctl, alternative="two-sided").pvalue)
        by = s.set_index(["gene", "aln_col"]).identity_to_human
        d = [by[(r.gene, r.contact_aln_col)] - by[(r.gene, r.control_aln_col)] for r in m.itertuples()]
        assert paired["n_a"] == len(m) and paired["estimate"] == pytest.approx(np.median(d))
        assert paired["p"] == pytest.approx(stats.wilcoxon(d).pvalue)

    def test_missing_tables_are_refused(self, tmp_path, monkeypatch):
        monkeypatch.setattr(complexes.Complex, "work", property(lambda self: tmp_path / "none"))
        with pytest.raises(SystemExit, match="missing"):
            readouts.main(("I",), paths.DATED_TREE, tmp_path / "r.tsv")
        assert not (tmp_path / "r.tsv").exists()
