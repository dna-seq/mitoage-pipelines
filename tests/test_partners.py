"""H12: the mito-nuclear contact pairs of complex IV and their discordance (`mitoport partners`)."""
from collections import Counter

import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.genes import COMPLEX_IV_NUCLEAR
from mitoage_pipelines.longevity.aa_pgls import select
from mitoage_pipelines.mitoport import partners as pt
from mitoage_pipelines.seqio import read_fasta

WORK = paths.MITOPORT_IFACE
ALN = paths.PARTNERS / "aln"
CIF = paths.PDB / "5Z62.cif"
needs_sites = pytest.mark.skipif(not (WORK / "control_pairs.tsv").is_file(), reason="needs `mitoport contacts` and `mitoport burial`")
needs_partners = pytest.mark.skipif(not (ALN / "COX4I1.afa").is_file(),
                                    reason="needs the vertebrate partner alignments (`mitoport iface-align` on partners/sel)")


@pytest.fixture(scope="module")
def sites():
    return pd.read_csv(WORK / "sites.tsv", sep="\t")


@pytest.fixture(scope="module")
def matched():
    return pd.read_csv(WORK / "control_pairs.tsv", sep="\t")


@needs_sites
class TestControlPairs:
    def test_every_contact_gets_its_own_control(self, sites, matched):
        ctl = pt.control_of(sites, matched)
        con = sites[sites.group.eq("contact")]
        assert set(ctl) == set(zip(con.gene, con.aln_col))
        assert len(set((g, c) for (g, _), c in ctl.items())) == len(ctl)              # no control used twice
        kept = sites[sites.group.eq("control")]
        assert set((g, c) for (g, _), c in ctl.items()) == set(zip(kept.gene, kept.aln_col))

    def test_a_control_of_another_gene_is_never_taken(self, sites, matched):
        by = sites.set_index(["gene", "aln_col"])
        for r in matched.itertuples():
            assert by.loc[(r.gene, r.contact_aln_col), "group"] == "contact"
            assert by.loc[(r.gene, r.control_aln_col), "group"] == "control"

    def test_a_sites_table_that_disagrees_with_the_pairing_is_refused(self, sites, matched):
        bad = sites.copy()
        i = bad.index[bad.group.eq("control")][0]
        j = bad.index[bad.group.eq("other") & bad.gene.eq(bad.gene[i]) & bad.aln_col.gt(0)][0]
        bad.loc[i, "group"], bad.loc[j, "group"] = "other", "control"
        with pytest.raises(SystemExit, match="one to one"):
            pt.control_of(bad, matched)


@pytest.fixture(scope="module")
def built():
    nuc = {g: read_fasta(ALN / f"{g}.afa") for g in COMPLEX_IV_NUCLEAR if (ALN / f"{g}.afa").is_file()}
    mt = {g: read_fasta(paths.BUILD / "alignments" / f"{g}.aa.fasta.gz") for g in pt.COX}
    s = pd.read_csv(WORK / "sites.tsv", sep="\t")
    col = pt.nuclear_columns(CIF, pd.read_csv(WORK / "chains.tsv", sep="\t"), nuc)
    pr = pt.pairs(pd.read_csv(WORK / "contacts.tsv", sep="\t"), s, col,
                  pt.control_of(s, pd.read_csv(WORK / "control_pairs.tsv", sep="\t")))
    return nuc, mt, col, pr


@needs_sites
@needs_partners
class TestPairs:
    def test_the_human_residue_at_each_nuclear_column_is_the_structures(self, built):
        from Bio.PDB import MMCIFParser
        from Bio.PDB.Polypeptide import protein_letters_3to1
        nuc, _, col, _ = built
        model = MMCIFParser(QUIET=True).get_structure("x", str(CIF))[0]
        same = [protein_letters_3to1.get(model[ch][(" ", rn, " ")].get_resname().upper(), "X") == nuc[g]["homo-sapiens"][c]
                for (ch, rn), (g, c) in col.items()]
        assert len(same) > 500 and np.mean(same) > 0.95

    def test_pairs_are_contacts_on_the_mito_side_and_controls_are_not(self, built, sites):
        _, _, _, pr = built
        con = set(zip(sites.gene[sites.group.eq("contact")], sites.aln_col[sites.group.eq("contact")] - 1))
        ctl = set(zip(sites.gene[sites.group.eq("control")], sites.aln_col[sites.group.eq("control")] - 1))
        assert set(zip(pr.mt_gene, pr.mt_col)) <= con and set(zip(pr.mt_gene, pr.ctl_col)) <= ctl
        assert not pr.duplicated().any() and len(pr) > 300


@pytest.fixture(scope="module")
def birds(built):
    nuc, mt, _, pr = built
    sp, _ = select(paths.BUILD, paths.LONGEVITY_OUT / "composition.tsv", paths.LONGEVITY_OUT / "fourfold.tsv",
                   paths.DATED_TREE, "Aves")
    slugs = [s for s in sp.slug if any(s in a for a in nuc.values())]
    return slugs, pt.discordance(pr, mt, nuc, slugs)


@needs_sites
@needs_partners
class TestDiscordance:

    def test_discordance_by_hand(self, built, birds):
        """Exactly one side off its consensus, counted with Counter over the birds, for a sample of species."""
        nuc, mt, _, pr = built
        slugs, got = birds
        def mode(aln, col):
            v = sorted(aln[s][col] for s in slugs if s in aln and aln[s][col] != "-")
            return Counter(v).most_common(1)[0][0] if v else "-"
        for s in slugs[:: max(1, len(slugs) // 8)]:
            d = n = 0
            for r in pr.itertuples():
                if s not in nuc[r.nuc_gene]:
                    continue
                a, b = mt[r.mt_gene][s][r.mt_col], nuc[r.nuc_gene][s][r.nuc_col]
                if a == "-" or b == "-":
                    continue
                n += 1; d += (a != mode(mt[r.mt_gene], r.mt_col)) ^ (b != mode(nuc[r.nuc_gene], r.nuc_col))
            assert got.loc[s, "contact_pairs"] == n
            if n >= pt.MIN_SHARE * len(pr):
                assert got.loc[s, "contact"] == pytest.approx(d / n, abs=1e-12)
            else:
                assert np.isnan(got.loc[s, "contact"])

    def test_a_species_filling_under_half_the_pairs_is_not_scored(self, built, birds):
        """Real species and alignments, with the partners of one subunit only: every species then fills fewer than
        half the contact pairs, and none may be scored."""
        nuc, mt, _, pr = built
        slugs, _ = birds
        one = pr.nuc_gene.value_counts().index[-1]                      # the subunit with the fewest pairs
        assert (pr.nuc_gene == one).mean() < pt.MIN_SHARE
        few = {g: (a if g == one else {}) for g, a in nuc.items()}
        d = pt.discordance(pr, mt, few, slugs[:40])
        assert d.contact_pairs.gt(0).any() and d.contact_pairs.lt(pt.MIN_SHARE * len(pr)).all()
        assert d.contact.isna().all() and d["contact minus control"].isna().all()

    def test_the_test_variable_is_contact_minus_control(self, birds):
        _, got = birds
        assert np.allclose(got["contact minus control"], got.contact - got.control, equal_nan=True)
        assert got.contact.dropna().between(0, 1).all()
