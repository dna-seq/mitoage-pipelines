"""H19: `mitoport classes`, the coupling class by class, the distance from the core and the mitochondrial side.

The structures, the build's alignments (proteins and codons) and the dated tree are real. Nuclear ortholog alignments
are not in the repository's inputs, so the tests make them from the build's mitochondrial alignments: a made "nuclear"
protein carries, at chosen columns, the residues real species have at chosen columns of a mitochondrial protein. Its
distances are then real distances, and what the code should find is known: a protein that carries complex III's
control residues is coupled with complex III at exactly 1. Expected values are computed in the test from the files."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from scipy import stats
from scipy.spatial.distance import cdist
from statsmodels.stats.multitest import multipletests

from mitoage_pipelines import paths
from mitoage_pipelines.genes import COMPLEX_II, CYTOSOLIC_CONTROLS, INNER_MEMBRANE, MATRIX, MT_GENES, MT_HEAVY_STRAND
from mitoage_pipelines.mitoport import burial, classes, complexes, contacts, iface
from mitoage_pipelines.seqio import read_fasta

NAMES = ("III", "IV")
CIFS = {n: paths.PDB / f"{complexes.get(n).pdb}.cif" for n in NAMES}
pytestmark = pytest.mark.skipif(
    not all(p.is_file() for p in CIFS.values()) or not (paths.BUILD / "alignments" / "ATP8.codon.fasta.gz").is_file()
    or not paths.DATED_TREE.is_file() or not paths.QUALITY_OK.is_file(),
    reason="needs `download build`, `download pdb 5Z62 5XTE` and `longevity dated-tree`")


@pytest.fixture(scope="module")
def mt():
    return {g: read_fasta(paths.BUILD / "alignments" / f"{g}.aa.fasta.gz") for g in MT_GENES}


@pytest.fixture(scope="module")
def codons():
    return {g: read_fasta(paths.BUILD / "alignments" / f"{g}.codon.fasta.gz") for g in MT_HEAVY_STRAND}


@pytest.fixture(scope="module")
def tree():
    return classes.shared_paths(paths.DATED_TREE)


@pytest.fixture(scope="module")
def maps(tmp_path_factory):
    out = {}
    for name in NAMES:
        cx = complexes.get(name)
        w = tmp_path_factory.mktemp(f"map{name}")
        contacts.main(CIFS[name], paths.BUILD, w, cx)
        burial.main(CIFS[name], w / "contacts.tsv", paths.BUILD, paths.DATED_TREE, paths.QUALITY_OK, w, cx)
        out[name] = w
    return out


def read(work: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(work / name, sep="\t", dtype={"chain": str, "partner_chain": str})


def entity_rows(cif: Path, genes: dict[str, str]) -> dict[str, str]:
    """{gene: its entities' canonical sequences joined, in entity order}: a stand-in for the human ortholog."""
    d = MMCIF2Dict(str(cif))
    seq = {}
    for chains, s in zip(d["_entity_poly.pdbx_strand_id"], d["_entity_poly.pdbx_seq_one_letter_code_can"]):
        g = genes[chains.split(",")[0].strip()]
        seq[g] = seq.get(g, "") + s.replace("\n", "")
    return seq


def of_class(cls: str, mt, codons, tips) -> list[str]:
    sp = pd.read_parquet(paths.BUILD / "species.parquet")
    have = set.intersection(*(set(a) for a in mt.values()), *(set(a) for a in codons.values())) & set(tips)
    return sorted(set(sp.slug[sp["class"].eq(cls)]) & have)


COPY, NOISY, FEW = "CYC1", "COX4I1", "COX5A"      # carries complex III's controls; noise at its interface; 39 mammals
HEADLESS = "UQCRH"                                 # its alignment holds no human protein
CONTROLS = CYTOSOLIC_CONTROLS[:6]
SOURCES = ["ND1", "ND2", "ND4", "ND5", "ATP6", "ND3", "COX1", "COX2", "CYTB", "COX3", "ND4L", "ATP8", "ND6"]


@pytest.fixture(scope="module")
def made(maps, mt, codons, tree, tmp_path_factory):
    """The made alignments: folder, the species of each class, the interface columns, and what each protein carries."""
    tips, _ = tree
    sites = read(maps["III"], "sites.tsv")
    cytb_controls = sites[sites.gene.eq("CYTB") & sites.group.eq("control")].aln_col.to_numpy() - 1

    def whole(cls):          # the class's species that fill every control column of cytochrome b
        return [s for s in of_class(cls, mt, codons, tips) if all(mt["CYTB"][s][c] != "-" for c in cytb_controls)]
    species = {"Mammalia": whole("Mammalia")[::12][:48], "Aves": whole("Aves")[::8][:44],
               "Reptilia": whole("Reptilia")[::7][:12]}
    assert [len(v) for v in species.values()] == [48, 44, 12]
    if "homo-sapiens" not in species["Mammalia"]:
        species["Mammalia"][0] = "homo-sapiens"
        species["Mammalia"].sort()
    everyone = [s for v in species.values() for s in v]
    human, cols = {}, {}
    folder = tmp_path_factory.mktemp("nucaln")
    for name in NAMES:
        cx = complexes.get(name)
        rows = entity_rows(CIFS[name], complexes.chain_genes(CIFS[name], cx))
        for g in cx.nuclear:
            human[g] = rows[g]
            (folder / f"{g}.afa").write_text(f">homo-sapiens\n{rows[g]}\n")
        ct = read(maps[name], "contacts.tsv")
        model = MMCIFParser(QUIET=True).get_structure("x", str(CIFS[name]))[0]
        cols.update(iface.interface_columns(model, read(maps[name], "chains.tsv"), ct[ct.partner_is_mt.eq(0)], folder)[1])
    carried = {}
    few_out = [s for s in species["Mammalia"] if s != "homo-sapiens"][-9:]         # FEW keeps 39 of the 48 mammals
    rng = np.random.default_rng(1)
    for k, g in enumerate(sorted(human)):
        free = [c for c in range(len(human[g])) if c not in set(cols[g])]
        src = "CYTB" if g == COPY else SOURCES[k % len(SOURCES)]
        hrow = mt[src]["homo-sapiens"]
        from_cols = cytb_controls if g == COPY else np.array([c for c in range(len(hrow)) if hrow[c] != "-"])
        assert g != COPY or len(free) >= len(from_cols)          # the copy must carry every control residue
        m = min(len(free), len(from_cols)); free, from_cols = free[:m], from_cols[:m]
        carried[g] = (src, from_cols, np.array(free))
        with open(folder / f"{g}.afa", "w") as f:
            for s in everyone:
                if (g == FEW and s in few_out) or (g == HEADLESS and s == "homo-sapiens"):
                    continue
                row = list(human[g])
                for c, p in zip(from_cols, free):
                    row[p] = mt[src][s][c]
                if g == NOISY and s != "homo-sapiens":          # anything at the interface: it must not be read
                    for c in cols[g]:
                        row[c] = "ACDEFGHIKLMNPQRSTVWY"[rng.integers(20)]
                f.write(f">{s}\n{''.join(row)}\n")
    others = [*CONTROLS, *INNER_MEMBRANE[:5], *COMPLEX_II, *MATRIX[:5]]
    for k, g in enumerate(others):                       # the build's alignment of a mitochondrial gene, as it is
        src = SOURCES[(k + 3) % len(SOURCES)]
        carried[g] = (src, None, None)
        (folder / f"{g}.afa").write_text("".join(f">{s}\n{mt[src][s]}\n" for s in everyone))
    return folder, species, cols, carried


@pytest.fixture(scope="module")
def ran(made, maps, tmp_path_factory):
    folder, *_ = made
    out = tmp_path_factory.mktemp("out")
    classes.main(paths.BUILD, folder, paths.DATED_TREE, out, names=NAMES, work=maps)
    return out


def table(out: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(out / name, sep="\t", keep_default_na=True)


# --- the pieces --------------------------------------------------------------------------------------
@pytest.fixture(scope="module")
def real(mt, codons, tree):
    """Real distances over pairs of 70 mammals: COX1, CYTB, ND5 and the tree's."""
    tips, H = tree
    sp = of_class("Mammalia", mt, codons, tips)[::8][:70]
    iu = np.triu_indices(len(sp), 1)
    x = iface.hamd(mt["COX1"], sp, np.arange(60, 400))[iu]
    y = iface.hamd(mt["CYTB"], sp, np.arange(30, 330))[iu]
    z = iface.hamd(mt["ND5"], sp, np.arange(50, 450))[iu]
    return x, y, z, classes.patristic(H, np.array([tips[s] for s in sp]))[iu]


class TestPartialCorrelation:

    def test_without_a_covariate_it_is_the_coupling_of_h14(self, real):
        x, y, _, T = real
        r, n = classes.pcorr(x, y, T)
        assert r == pytest.approx(iface.pcorr(x, y, T), abs=1e-12) and n == len(x) and 0 < r < 1

    def test_a_covariate_enters_as_its_rank_beside_time(self, real):
        x, y, z, T = real
        Z = sm.add_constant(np.column_stack([T, np.log1p(T), stats.rankdata(z)]))
        want = stats.pearsonr(sm.OLS(stats.rankdata(x), Z).fit().resid, sm.OLS(stats.rankdata(y), Z).fit().resid)[0]
        got = classes.pcorr(x, y, T, (z,))[0]
        assert got == pytest.approx(want, abs=1e-9)
        assert got < classes.pcorr(x, y, T)[0]                 # another mitochondrial gene takes part of it

    def test_a_side_the_covariate_already_holds_is_not_a_correlation(self, real):
        x, y, _, T = real
        assert np.isnan(classes.pcorr(x, y, T, (x * 3 + 1,))[0])

    def test_too_few_pairs_are_not_assessed_and_counted(self, real):
        x, y, _, T = real
        r, n = classes.pcorr(x[:199], y[:199], T[:199])
        assert np.isnan(r) and n == 199
        x2 = x.copy(); x2[: len(x) - 150] = np.nan
        assert classes.pcorr(x2, y, T) [1] == 150


class TestFourfoldTransversions:
    def test_only_shared_fourfold_families_count_and_only_transversions_differ(self):
        aln = {"a": "CTA" * 120 + "ATA" + "GGC" * 120 + "TTA",
               "b": "CTG" * 120 + "ATC" + "GGA" * 120 + "CTA",   # transitions; not fourfold; transversions; another family
               "c": "CTC" * 120 + "---" + "GGN" * 120 + "CTA"}
        D = classes.fourfold_tv({"g": aln}, ["a", "b", "c"])
        assert D[0, 1] == pytest.approx(120 / 240)                # 120 transitions, 120 transversions, 240 sites
        assert D[0, 2] == pytest.approx(120 / 120)                # GGN is no site: 120 sites, all transversions
        assert D[1, 2] == pytest.approx(120 / 121) and (D == D.T).all()     # the last codon, CTA in both, is a site
        short = {k: v[:150] for k, v in aln.items()}              # 50 shared sites: under the minimum
        assert np.isnan(classes.fourfold_tv({"g": short}, ["a", "b"])[0, 1])

    def test_a_codon_alignment_out_of_frame_is_refused(self):
        with pytest.raises(SystemExit, match="multiple of three"):
            classes.fourfold_tv({"g": {"a": "CTAC", "b": "CTAC"}}, ["a", "b"])

    def test_the_build_gives_what_counting_codon_by_codon_gives(self, mt, codons, tree):
        sp = of_class("Aves", mt, codons, tree[0])[::40][:6]
        got = classes.fourfold_tv(codons, sp)
        purine = set("AG")
        for i in range(len(sp)):
            for j in range(i + 1, len(sp)):
                sites = tv = 0
                for aln in codons.values():
                    a, b = aln[sp[i]].upper(), aln[sp[j]].upper()
                    for k in range(0, len(a), 3):
                        if a[k:k + 2] == b[k:k + 2] and a[k:k + 2] in classes.FOURFOLD_PREFIXES and a[k + 2] in "ACGT" and b[k + 2] in "ACGT":
                            sites += 1; tv += (a[k + 2] in purine) != (b[k + 2] in purine)
                assert sites > 1000 and got[i, j] == pytest.approx(tv / sites) and got[j, i] == got[i, j]
        assert 0 < np.nanmax(got) < 0.5 and (np.diag(got) == 0).all()


class TestCoreDistance:
    def test_it_is_the_median_distance_to_the_nearest_mitochondrial_atom(self):
        cx = complexes.get("IV")
        d = classes.core_distance(cx).set_index("gene")
        gene_of = complexes.chain_genes(cx.cif, cx)
        model = MMCIFParser(QUIET=True).get_structure("x", str(cx.cif))[0]
        core = np.array([a.coord for ch in model if gene_of[ch.id] in cx.mt for a in ch.get_atoms() if a.element != "H"])
        chain = next(c for c, g in gene_of.items() if g == "COX5A")
        mine = np.array([a.coord for a in model[chain].get_atoms() if a.element != "H"])
        nearest = cdist(mine, core).min(axis=1)
        assert d.loc["COX5A", "median_dist"] == pytest.approx(np.median(nearest), abs=1e-3)
        assert d.loc["COX5A", "min_dist"] == pytest.approx(nearest.min(), abs=1e-3)
        assert set(d.index) == set(cx.nuclear) and (d.min_dist <= d.median_dist).all()

    def test_a_gene_in_two_chains_takes_the_copy_nearer_the_core(self):
        cx = complexes.get("III")
        d = classes.core_distance(cx).set_index("gene")
        gene_of = complexes.chain_genes(cx.cif, cx)
        model = MMCIFParser(QUIET=True).get_structure("x", str(cx.cif))[0]
        core = np.array([a.coord for ch in model if gene_of[ch.id] == "CYTB" for a in ch.get_atoms() if a.element != "H"])
        medians = {c: np.median(cdist(np.array([a.coord for a in model[c].get_atoms() if a.element != "H"]), core).min(axis=1))
                   for c, g in gene_of.items() if g == "UQCRC2"}
        assert len(medians) == 2
        assert d.loc["UQCRC2", "chain"] == min(medians, key=medians.get)
        assert d.loc["UQCRC2", "median_dist"] == pytest.approx(min(medians.values()), abs=1e-3)

    @pytest.mark.parametrize("name", NAMES)
    def test_the_subunits_the_contact_map_names_are_within_the_cutoff(self, maps, name):
        """The contact map reads the reference chain of each mitochondrial gene; the distance counts every copy. So
        every subunit of the map is within the cutoff, and one beyond the map touches the other copy only."""
        cx = complexes.get(name)
        d = classes.core_distance(cx)
        ct, ch = read(maps[name], "contacts.tsv"), read(maps[name], "chains.tsv")
        gene = dict(zip(ch.chain, ch.gene))
        touching = {gene[c] for c in ct.partner_chain[ct.partner_is_mt.eq(0)]}
        within = set(d.gene[d.min_dist <= contacts.CUTOFF])
        assert 0 < len(touching) and touching <= within
        model = MMCIFParser(QUIET=True).get_structure("x", str(cx.cif))[0]
        reference = np.array([a.coord for c in cx.mt.values() for a in model[c].get_atoms() if a.element != "H"])
        for g in within - touching:
            nearest = min(cdist(np.array([a.coord for a in model[c].get_atoms() if a.element != "H"]), reference).min()
                          for c in ch.chain[ch.gene.eq(g)])
            assert nearest > contacts.CUTOFF
        if name == "IV":                                        # one copy of each gene: the two sets are one
            assert within == touching


class TestStratified:
    X = np.array([1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 7, 8], float)
    S = np.array(["a"] * 5 + ["b"] * 6 + ["c"] * 2)

    def test_ranks_are_taken_inside_each_complex(self):
        y = np.where(self.S == "a", -self.X, 1000 - self.X)       # falls with x inside each; across them it rises
        got = classes.stratified_corr(self.X, y, self.S)
        assert got["r"] == pytest.approx(-1) and got["n"] == 11 and got["strata"] == 2      # c has two items: left out
        assert got["p"] < 0.01
        assert stats.spearmanr(self.X[:11], y[:11])[0] > -0.9

    def test_the_p_value_is_the_share_of_permutations_at_least_as_far(self):
        y = np.array([3, 5, 1, 4, 2, 4, 1, 6, 2, 5, 3, 0, 0], float)
        got = classes.stratified_corr(self.X, y, self.S, seed=3)
        rx, ry = classes.within_ranks(self.X[:11], self.S[:11]), classes.within_ranks(y[:11], self.S[:11])
        assert got["r"] == pytest.approx(np.corrcoef(rx, ry)[0, 1])
        rng = np.random.default_rng(3); hits = 0
        for _ in range(classes.PERMUTATIONS):
            p = ry.copy()
            for s in ("a", "b"):
                i = np.flatnonzero(self.S[:11] == s); p[i] = ry[rng.permutation(i)]
            hits += abs(np.corrcoef(rx, p)[0, 1]) >= abs(got["r"]) - 1e-12
        assert got["p"] == (1 + hits) / (1 + classes.PERMUTATIONS) and 0.05 < got["p"] <= 1

    def test_a_covariate_that_holds_one_side_leaves_nothing(self):
        y = np.arange(13, dtype=float)[::-1]
        assert classes.stratified_corr(self.X, y, self.S, covariate=self.X * 2)["n"] == 0
        z = np.array([3, 1, 2, 5, 4, 1, 3, 2, 4, 6, 5, 1, 2], float)
        with_z = classes.stratified_corr(self.X, y, self.S, covariate=z)
        assert with_z["n"] == 11 and -1 <= with_z["r"] < 0

    def test_missing_values_and_small_complexes_are_left_out(self):
        y = -self.X; y[[0, 1, 2]] = np.nan                         # a keeps two items: out
        got = classes.stratified_corr(self.X, y, self.S)
        assert got["n"] == 6 and got["strata"] == 1
        assert classes.stratified_corr(self.X[:4], y[:4], self.S[:4]) == {"n": 0}

    def test_the_difference_is_between_ranks_inside_each_complex(self):
        v = np.array([5, 6, 1, 2, 3, 50, 60, 70, 10, 20, 30, 1, 2], float)
        lab = np.array([1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0], bool)
        got = classes.stratified_diff(v, lab, self.S)
        r = classes.within_ranks(v[:11], self.S[:11])
        assert got["diff"] == pytest.approx(r[lab[:11]].mean() - r[~lab[:11]].mean()) and got["diff"] > 0
        assert (got["n_a"], got["n_b"], got["strata"]) == (5, 6, 2) and got["p"] < 0.05
        assert classes.stratified_diff(v, np.ones(13, bool), self.S) == {"n": 0}


class TestSpecies:
    def test_a_species_needs_every_alignment_a_tip_and_a_nuclear_sequence(self, mt, codons, tree):
        tips, _ = tree
        sp = pd.read_parquet(paths.BUILD / "species.parquet")
        birds = of_class("Aves", mt, codons, tips)
        nuc = {"x": {s: "A" for s in birds[:30]}, "y": {s: "A" for s in birds[20:50]}}
        assert classes.class_species(sp, "Aves", mt, codons, tips, nuc) == birds[:50]
        drop = birds[5]
        less = {**mt, "ND6": {k: v for k, v in mt["ND6"].items() if k != drop}}
        assert drop not in classes.class_species(sp, "Aves", less, codons, tips, nuc)
        less = {**codons, "ND1": {k: v for k, v in codons["ND1"].items() if k != drop}}
        assert drop not in classes.class_species(sp, "Aves", mt, less, tips, nuc)
        assert drop not in classes.class_species(sp, "Aves", mt, codons, {k: v for k, v in tips.items() if k != drop}, nuc)
        assert not set(classes.class_species(sp, "Mammalia", mt, codons, tips, nuc)) & set(birds)


# --- the whole run on the made alignments --------------------------------------------------------------
class TestRun:
    def test_a_protein_carrying_complex_iiis_controls_is_coupled_with_it_at_one(self, ran):
        t = table(ran, "class_coupling.tsv")
        for cls in ("Mammalia", "Aves"):
            row = t[t["class"].eq(cls) & t.gene.eq(COPY) & t.mt.eq("III")].iloc[0]
            assert row.status == "computed" and row.coupling == pytest.approx(1, abs=1e-4)
            assert row.given_other == pytest.approx(1, abs=1e-4) and np.isnan(row.other_given)
            other = t[t["class"].eq(cls) & t.gene.eq(COPY) & t.mt.eq("IV")].iloc[0]
            assert other.coupling < 0.99

    def test_the_interface_columns_are_not_read(self, made, maps, ran, tmp_path):
        folder, species, cols, _ = made
        assert len(cols[NOISY]) >= iface.MIN_IFACE
        clean = tmp_path / "aln"; clean.mkdir()
        for p in folder.glob("*.afa"):
            (clean / p.name).write_bytes(p.read_bytes())
        human = read_fasta(folder / f"{NOISY}.afa")["homo-sapiens"]
        rows = read_fasta(folder / f"{NOISY}.afa")
        for s, row in rows.items():                              # the human residue back at every interface column
            row = list(row)
            for c in cols[NOISY]:
                row[c] = human[c]
            rows[s] = "".join(row)
        (clean / f"{NOISY}.afa").write_text("".join(f">{s}\n{r}\n" for s, r in rows.items()))
        classes.main(paths.BUILD, clean, paths.DATED_TREE, tmp_path / "out", classes=("Mammalia",), names=NAMES, work=maps)
        a, b = table(ran, "class_coupling.tsv"), table(tmp_path / "out", "class_coupling.tsv")
        a = a[a["class"].eq("Mammalia") & a.gene.eq(NOISY)].reset_index(drop=True)
        b = b[b.gene.eq(NOISY)].reset_index(drop=True)
        pd.testing.assert_frame_equal(a, b)
        assert (a.iface_cols == len(cols[NOISY])).all() and a.coupling.notna().all()

    def test_a_control_proteins_coupling_is_what_the_files_give(self, made, maps, mt, codons, tree, ran):
        _, species, _, carried = made
        tips, H = tree
        g = CONTROLS[2]; src = carried[g][0]
        sp = species["Aves"]
        hrow = mt[src]["homo-sapiens"]; filled = [c for c in range(len(hrow)) if hrow[c] != "-"]
        sites = read(maps["IV"], "sites.tsv")

        def share(aln, a, b, cols):
            both = [c for c in cols if aln[a][c] != "-" and aln[b][c] != "-"]
            return sum(aln[a][c] != aln[b][c] for c in both) / len(both) if len(both) >= 5 else np.nan
        per_gene = []
        for cox in ("COX1", "COX2", "COX3"):
            cols = sites[sites.gene.eq(cox) & sites.group.eq("control")].aln_col.to_numpy() - 1
            D = np.array([[share(mt[cox], a, b, cols) for b in sp] for a in sp])
            per_gene.append(D / np.nanmean(D))
        iu = np.triu_indices(len(sp), 1)
        dm = np.nanmean(np.stack(per_gene), axis=0)[iu]
        dn = np.array([[share(mt[src], a, b, filled) for b in sp] for a in sp])[iu]
        T = classes.patristic(H, np.array([tips[s] for s in sp]))[iu]
        Z = sm.add_constant(np.column_stack([T, np.log1p(T)]))
        want = stats.pearsonr(sm.OLS(stats.rankdata(dn), Z).fit().resid, sm.OLS(stats.rankdata(dm), Z).fit().resid)[0]
        t = table(ran, "class_coupling.tsv")
        row = t[t["class"].eq("Aves") & t.gene.eq(g) & t.mt.eq("IV")].iloc[0]
        assert row.coupling == pytest.approx(want, abs=1e-4)

        def partial(x, y, covariate):                           # x and y given time and the rank of the covariate
            Zc = np.column_stack([Z, stats.rankdata(covariate)])
            return stats.pearsonr(sm.OLS(stats.rankdata(x), Zc).fit().resid, sm.OLS(stats.rankdata(y), Zc).fit().resid)[0]
        s3 = read(maps["III"], "sites.tsv")                    # the other complex read here is III: cytochrome b
        cols = s3[s3.gene.eq("CYTB") & s3.group.eq("control")].aln_col.to_numpy() - 1
        do = np.array([[share(mt["CYTB"], a, b, cols) for b in sp] for a in sp])[iu]
        assert row.given_other == pytest.approx(partial(dn, dm, do), abs=1e-4)
        assert row.other_given == pytest.approx(partial(dn, do, dm), abs=1e-4)
        assert abs(row.given_other - row.coupling) > 0.01
        both = np.nanmean(np.stack([*per_gene, (lambda D: D / np.nanmean(D))(
            np.array([[share(mt["CYTB"], a, b, cols) for b in sp] for a in sp]))]), axis=0)[iu]
        s4 = classes.fourfold_tv({k: {x: v[x] for x in sp} for k, v in codons.items()}, sp)[iu]
        everything = t[t["class"].eq("Aves") & t.gene.eq(g) & t.mt.eq("all")].iloc[0]
        assert everything.tv4 == pytest.approx(stats.pearsonr(sm.OLS(stats.rankdata(dn), Z).fit().resid,
                                                              sm.OLS(stats.rankdata(s4), Z).fit().resid)[0], abs=1e-4)
        assert everything.aa_given_tv4 == pytest.approx(partial(dn, both, s4), abs=1e-4)
        assert everything.tv4_given_aa == pytest.approx(partial(dn, s4, both), abs=1e-4)
        assert abs(everything.aa_given_tv4 - everything.coupling) > 0.01
        assert (row["species"], row["pairs"], row["columns"], row["iface_cols"]) == (len(sp), len(iu[0]), len(filled), 0)
        assert row.divergence == pytest.approx(np.nanmean(dn), abs=1e-4)

    def test_what_cannot_be_computed_says_why_and_holds_no_number(self, ran):
        t = table(ran, "class_coupling.tsv")
        few = t[t["class"].eq("Mammalia") & t.gene.eq(FEW)]
        assert len(few) == 3 and (few.status == "not assessed: 39 species").all() and few[classes.VALUES].isna().all().all()
        rep = t[t["class"].eq("Reptilia")]
        assert len(rep) and rep.status.str.startswith("not assessed: 12 species of the class").all()
        assert rep[classes.VALUES].isna().all().all()
        absent = t[t.gene.eq(CYTOSOLIC_CONTROLS[10]) & t["class"].isin(["Mammalia", "Aves"])]
        assert len(absent) == 6 and (absent.status == "not assessed: no alignment").all()
        headless = t[t.gene.eq(HEADLESS) & t["class"].isin(["Mammalia", "Aves"])]
        assert len(headless) == 6 and (headless.status == "not assessed: the alignment holds no human protein").all()
        assert headless[classes.VALUES].isna().all().all()
        fish = t[t["class"].eq("Actinopterygii")]
        assert fish.status.str.startswith("not assessed: 0 species of the class").all()
        r = table(ran, "class_readouts.tsv")
        rr = r[r["class"].eq("Reptilia")]
        assert len(rr) == 1 and rr.status.iloc[0].startswith("not assessed") and rr.p.isna().all()
        o1 = r[r.readout.eq("O1") & r["class"].eq("Aves")].iloc[0]          # complex II has four proteins
        assert o1.status == "not assessed: 4 and 6 proteins with a value" and np.isnan(o1.p)

    def test_every_protein_has_a_row_per_class_and_set_in_its_group(self, ran):
        t = table(ran, "class_coupling.tsv")
        want = {g: classes.subunit_group(n) for n in NAMES for g in complexes.get(n).nuclear}
        for grp, genes in classes.OTHER_GROUPS.items():
            want.update({g: grp for g in genes})
        assert dict(zip(t.gene, t.group)) == want
        assert t.groupby(["class", "gene"]).mt.apply(lambda v: sorted(v) == ["III", "IV", "all"]).all()
        assert set(t["class"]) == {"Mammalia", "Aves", "Actinopterygii", "Reptilia"}
        sub = t[t.complex.notna()]
        assert set(sub.complex) == set(NAMES) and sub.group.str.endswith("subunit").all()
        rate = t[t.status.eq("computed")]
        assert rate[rate.mt.eq("all")][["tv4", "aa_given_tv4", "tv4_given_aa"]].notna().all().all()
        assert rate[~rate.mt.eq("all")][["tv4", "aa_given_tv4", "tv4_given_aa"]].isna().all().all()
        assert rate[rate.mt.eq("all")][["given_other", "other_given"]].isna().all().all()

    def test_the_readouts_are_the_tests_of_the_table(self, ran):
        t, r = table(ran, "class_coupling.tsv"), table(ran, "class_readouts.tsv")
        d = table(ran, "core_distance.tsv")
        for cls in ("Mammalia", "Aves"):
            x = t[t["class"].eq(cls)]
            a = x[x.group.eq("complex IV subunit") & x.mt.eq("IV")].given_other.dropna()
            b = x[x.group.eq(classes.CONTROL) & x.mt.eq("IV")].given_other.dropna()
            got = r[r["class"].eq(cls) & r.readout.eq("S1") & r.set.eq("complex IV")].iloc[0]
            assert got.p == pytest.approx(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue, rel=1e-4)
            assert (got.n_a, got.n_b) == (len(a), len(b)) and got.estimate == pytest.approx(a.median() - b.median(), abs=1e-4)
            own = x[x.complex.eq(x.mt)].merge(d, on=["complex", "gene"])
            want = classes.stratified_corr(own.median_dist, own.coupling, own.complex.to_numpy(),
                                           seed=classes.SEED + 10 * classes.CLASSES.index(cls))
            g1 = r[r["class"].eq(cls) & r.readout.eq("G1")].iloc[0]
            assert g1.estimate == pytest.approx(want["r"], abs=1e-4) and g1.p == pytest.approx(want["p"], rel=1e-4)
            assert g1.n_a == want["n"] == own.coupling.notna().sum()
            want = classes.stratified_corr(own.median_dist, own.coupling, own.complex.to_numpy(), own.divergence.to_numpy(),
                                           seed=classes.SEED + 10 * classes.CLASSES.index(cls) + 1)
            g2 = r[r["class"].eq(cls) & r.readout.eq("G2")].iloc[0]
            assert g2.estimate == pytest.approx(want["r"], abs=1e-4) and g2.p == pytest.approx(want["p"], rel=1e-4)
            assert abs(g2.estimate - g1.estimate) > 1e-3
            two = own[(own.iface_cols >= iface.MIN_IFACE) | (own.iface_cols == 0)]
            assert 0 < len(two) < len(own)                                      # some subunits have one to four columns
            want = classes.stratified_diff(two.coupling, two.iface_cols >= iface.MIN_IFACE, two.complex.to_numpy(),
                                           seed=classes.SEED + 10 * classes.CLASSES.index(cls) + 2)
            g3 = r[r["class"].eq(cls) & r.readout.eq("G3")].iloc[0]
            assert g3.estimate == pytest.approx(want["diff"], abs=1e-4) and (g3.n_a, g3.n_b) == (want["n_a"], want["n_b"])
            every = x[x.mt.eq("all")]
            a = every[every.complex.notna()].aa_given_tv4.dropna()
            b = every[every.group.eq(classes.CONTROL)].aa_given_tv4.dropna()
            t2 = r[r["class"].eq(cls) & r.readout.eq("T2")].iloc[0]
            assert t2.p == pytest.approx(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue, rel=1e-4)
            a = every[every.complex.notna() & (every.iface_cols >= iface.MIN_IFACE)].coupling.dropna()
            b = every[every.group.eq("inner membrane, not respiratory")].coupling.dropna()
            m2 = r[r["class"].eq(cls) & r.readout.eq("M2")].iloc[0]
            assert (m2.n_a, m2.n_b) == (len(a), len(b)) and m2.p == pytest.approx(
                stats.mannwhitneyu(a, b, alternative="two-sided").pvalue, rel=1e-4)
            c1 = r[r["class"].eq(cls) & r.readout.eq("C1") & r.set.eq("complex III")].iloc[0]
            a = x[x.group.eq("complex III subunit") & x.mt.eq("III")].coupling.dropna()
            assert c1.median_a == pytest.approx(a.median(), abs=1e-4) and c1.n_a == len(a)
            tested = r[r["class"].eq(cls) & r.status.eq("tested")]
            prim = tested[tested.readout.isin(classes.PRIMARY)]
            assert len(prim) == 2 + 1 + 2 + 2                                   # C1 and S1 per complex, G1, T2, T3
            assert np.allclose(prim.p_holm, multipletests(prim.p, method="holm")[1], rtol=1e-4)
            assert tested[~tested.readout.isin(classes.PRIMARY)].p_holm.isna().all()

    def test_the_matrix_holds_the_medians(self, ran):
        t, m = table(ran, "class_coupling.tsv"), table(ran, "class_matrix.tsv")
        ok = t[t.status.eq("computed")]
        for row in m.itertuples():
            x = ok[ok["class"].eq(row._1) & ok.group.eq(row.group) & ok.mt.eq(row.mt)]
            assert row.proteins == len(x) and row.median_coupling == pytest.approx(x.coupling.median(), abs=1e-4)
        assert len(m) == ok.groupby(["class", "group", "mt"]).ngroups and "Reptilia" not in set(m["class"])

    def test_a_second_run_writes_the_same_bytes(self, made, maps, ran, tmp_path):
        classes.main(paths.BUILD, made[0], paths.DATED_TREE, tmp_path, names=NAMES, work=maps)
        for f in ("core_distance.tsv", "class_coupling.tsv", "class_matrix.tsv", "class_readouts.tsv"):
            assert (tmp_path / f).read_bytes() == (ran / f).read_bytes()

    def test_a_missing_contact_map_is_refused(self, made, tmp_path):
        with pytest.raises(SystemExit, match="needs complex III's contact map"):
            classes.main(paths.BUILD, made[0], paths.DATED_TREE, tmp_path, names=NAMES, work={n: tmp_path for n in NAMES})
