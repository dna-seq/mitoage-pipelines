"""H36: the RefSeq orthologs of the cybrid hosts and donors (`mitoport/donor_orthologs.py`): what NCBI's answers are read as, the
choice of the ortholog gene and of its protein, the pinned table, the download that refuses a changed sequence and the reading
the swap makes of it.

The fixtures in tests/data/donor_orthologs/ are NCBI's own answers of 2026-11-24 for the human genes NDUFS2 (4720), ATP5MC1 (516) and
NDUFB1 (4707) (NCBI Datasets `gene/id/<id>/orthologs` for the mouse, rat, M. caroli, M. pahari and the Chinese hamster, and `product_report`
of the genes found), `gene/symbol/NDUFS2/taxon/human`, `genome/taxon/<id>/dataset_report?filters.assembly_source=refseq` for the Chinese
hamster (two RefSeq assemblies, one of them the reference), M. spretus and M. terricolor (none), and the 38 RefSeq proteins of those
genes from efetch, with the fields these tests do not read (exons, ontology, names) taken out. NDUFB1 has no ortholog in the hamster's
group; mouse ATP5MC1 has three isoforms of one length. NCBI's ortholog group holds one gene of a species at most (no human gene of the
72 had two), so the ambiguous case is made from two real gene records put in one group, as the test says."""
import copy
import json
from pathlib import Path

import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.mitoport import donor_orthologs as do
from mitoage_pipelines.mitoport.complexes import COMPLEXES

FIX = Path(__file__).parent / "data" / "donor_orthologs"
ASSEMBLY = {"mus-musculus": "GCF_000001635.27", "rattus-norvegicus": "GCF_036323735.1", "mus-caroli": "GCF_978045245.1",
            "mus-pahari": "GCF_900095145.1", "cricetulus-griseus": "GCF_003668045.3"}
TAXA = {s: do.SPECIES[s] for s in ASSEMBLY}


def load(name):
    return json.loads((FIX / name).read_text())


def found(sym):
    return do.ortholog_genes([load(f"{sym}.orthologs.json")])


def products(sym):
    return {int(r["product"]["gene_id"]): r["product"] for r in load(f"{sym}.product_report.json")["reports"]}


PROTEINS = do.parse_fasta((FIX / "proteins.faa").read_text())


def rows_of(sym, hid):
    f = found(sym)
    return do.gene_rows(sym, hid, f, do.chosen_genes(f, TAXA, ASSEMBLY), products(sym), TAXA, ASSEMBLY)


def fixture_table():
    t = pd.DataFrame(rows_of("NDUFS2", 4720) + rows_of("ATP5MC1", 516) + rows_of("NDUFB1", 4707), columns=do.COLUMNS)
    return do.fill_sequences(t, fetch=lambda accs: {a: PROTEINS[a] for a in accs})


class TestAssembly:
    def test_the_reference_genome_among_the_refseq_assemblies_is_the_source(self):
        rep = load("genome_10029.json")
        assert len(rep["reports"]) == 2 and do.reference_assembly(rep) == "GCF_003668045.3"          # CriGri_1.0 is no reference

    def test_a_species_without_a_refseq_annotation_has_no_source_and_nothing_stands_in(self):
        assert do.reference_assembly(load("genome_10096.json")) is None and do.reference_assembly(load("genome_254704.json")) is None
        assert do.reference_assembly({}) is None                                      # NCBI answers {} for a taxon with none

    def test_a_refseq_assembly_that_is_not_the_reference_and_a_genbank_one_are_not_taken(self):
        rep = copy.deepcopy(load("genome_10029.json"))
        for r in rep["reports"]:
            r["assembly_info"]["refseq_category"] = None
        assert do.reference_assembly(rep) is None
        rep = copy.deepcopy(load("genome_10029.json"))
        for r in rep["reports"]:
            r["accession"] = r["accession"].replace("GCF_", "GCA_")
        assert do.reference_assembly(rep) is None

    def test_two_reference_genomes_are_refused(self):
        rep = copy.deepcopy(load("genome_10029.json"))
        for r in rep["reports"]:
            r["assembly_info"]["refseq_category"] = "reference genome"
        with pytest.raises(ValueError, match="more than one"):
            do.reference_assembly(rep)


class TestOrthologGenes:
    def test_the_human_gene_is_the_one_with_exactly_the_symbol(self):
        assert do.human_gene_id(load("NDUFS2.human.json"), "NDUFS2") == 4720
        with pytest.raises(ValueError, match="found 0"):
            do.human_gene_id(load("NDUFS2.human.json"), "NDUFS")

    def test_the_group_is_split_by_taxon(self):
        f = found("NDUFS2")
        assert set(f) == set(TAXA.values()) and all(len(g) == 1 for g in f.values())
        assert {int(g[0]["gene_id"]) for g in f.values()} == {226646, 289218, 110305317, 110321406, 100768067}

    def test_a_species_with_no_gene_in_the_group_has_no_ortholog(self):
        f = found("NDUFB1")
        assert TAXA["cricetulus-griseus"] not in f
        assert do.choose_gene(f.get(TAXA["cricetulus-griseus"], []), ASSEMBLY["cricetulus-griseus"]) == (None, "no ortholog")
        assert do.choose_gene(f[TAXA["mus-musculus"]], ASSEMBLY["mus-musculus"])[1] == do.OK

    def test_two_genes_of_a_species_in_the_group_are_ambiguous_and_none_is_taken(self):
        mouse_a, mouse_b = found("NDUFS2")[10090][0], found("ATP5MC1")[10090][0]          # two real mouse genes, put in one group
        group = do.ortholog_genes([{"reports": [{"gene": mouse_a}, {"gene": mouse_b}]}])
        assert do.choose_gene(group[10090], ASSEMBLY["mus-musculus"]) == (None, "ambiguous (paralogs)")
        rows = do.gene_rows("NDUFS2", 4720, group, do.chosen_genes(group, {"mus-musculus": 10090}, {"mus-musculus": ASSEMBLY["mus-musculus"]}),
                            {}, {"mus-musculus": 10090}, {"mus-musculus": ASSEMBLY["mus-musculus"]})
        assert rows[0]["status"] == "ambiguous (paralogs)" and rows[0]["protein"] == ""
        assert rows[0]["candidates"] == "11951:Atp5mc1:PROTEIN_CODING;226646:Ndufs2:PROTEIN_CODING"

    def test_the_candidates_are_in_gene_id_order(self):
        a, b = found("NDUFS2")[10090][0], found("ATP5MC1")[10090][0]
        got = do.ortholog_genes([{"reports": [{"gene": a}, {"gene": b}]}])[10090]
        assert [int(g["gene_id"]) for g in got] == sorted(int(g["gene_id"]) for g in (a, b)) and int(got[0]["gene_id"]) == 11951

    def test_a_gene_that_is_not_protein_coding_or_not_on_the_assembly_is_not_taken(self):
        g = copy.deepcopy(found("NDUFS2")[10090][0])
        assert do.choose_gene([g], ASSEMBLY["mus-musculus"])[1] == do.OK
        assert do.choose_gene([g], ASSEMBLY["rattus-norvegicus"]) == (None, "ortholog not annotated on the assembly")
        g["type"] = "PSEUDO"
        assert do.choose_gene([g], ASSEMBLY["mus-musculus"]) == (None, "ortholog is not protein coding")


class TestProteins:
    def test_the_proteins_of_a_gene_on_each_assembly_it_is_annotated_on_are_told_apart(self):
        h = found("NDUFS2")[10029][0]
        prod = products("NDUFS2")[100768067]
        on = {a["assembly_accession"]: do.assembly_accessions(h, a["assembly_accession"]) for a in h["annotations"]}
        assert set(on) == {"GCF_000223135.1", "GCF_003668045.3"} and on["GCF_000223135.1"].isdisjoint(on["GCF_003668045.3"])
        assert [p["accession"] for p in do.refseq_proteins(prod, on["GCF_003668045.3"])] == ["XP_027274810.1"]
        assert [p["accession"] for p in do.refseq_proteins(prod, on["GCF_000223135.1"])] == ["XP_016833292.1"]
        assert do.refseq_proteins(prod, {"NC_000000.1"}) == []
        assert do.assembly_accessions(h, "GCF_000001635.27") == set()

    def test_the_longest_isoform_is_taken_and_each_protein_is_counted_once(self):
        m = found("NDUFS2")[10090][0]
        prots = do.refseq_proteins(products("NDUFS2")[226646], do.assembly_accessions(m, ASSEMBLY["mus-musculus"]))
        assert len(prots) == 3 and len({p["accession"] for p in prots}) == 3
        top, tied = do.choose_protein(prots)
        assert top["length"] == max(p["length"] for p in prots) == 463 and top["accession"] == "NP_694704.1" and tied == ["XP_030109284.1"]
        assert all(p["length"] == len(PROTEINS[p["accession"]]) for p in prots)          # NCBI's lengths are the sequences'

    def test_isoforms_of_one_length_are_told_apart_by_accession_whatever_the_order(self):
        m = found("ATP5MC1")[10090][0]
        prots = do.refseq_proteins(products("ATP5MC1")[11951], do.assembly_accessions(m, ASSEMBLY["mus-musculus"]))
        assert len(prots) == 3 and {p["length"] for p in prots} == {136}                 # three transcripts, one length
        for order in (prots, prots[::-1], prots[1:] + prots[:1]):
            top, tied = do.choose_protein(order)
            assert (top["accession"], tied) == ("NP_031532.2", ["NP_001154891.1", "NP_001404256.1"])
        assert len({PROTEINS[p["accession"]] for p in prots}) == 1                        # and one sequence

    def test_a_shorter_protein_of_the_gene_is_never_taken_and_a_longer_one_always_is(self):
        m = found("NDUFS2")[10090][0]
        prots = do.refseq_proteins(products("NDUFS2")[226646], do.assembly_accessions(m, ASSEMBLY["mus-musculus"]))
        longer = {"accession": "XP_999999999.1", "length": 464, "genomic": "x"}
        assert do.choose_protein([*prots, longer])[0] == longer and do.choose_protein(prots)[0]["length"] == 463
        assert min(p["length"] for p in prots) < 463

    def test_curated_before_predicted_and_numbers_in_order_not_as_text(self):
        mk = lambda *accs: [{"accession": a, "length": 50} for a in accs]
        assert do.choose_protein(mk("XP_000100.1", "NP_999999.1"))[0]["accession"] == "NP_999999.1"
        assert do.choose_protein(mk("NP_10.1", "NP_9.1"))[0]["accession"] == "NP_9.1"            # as text, "NP_10.1" comes first
        assert do.choose_protein(mk("NP_9.2", "NP_9.1"))[0]["accession"] == "NP_9.1"
        assert do.choose_protein([{"accession": "NP_1.1", "length": 10}, {"accession": "XP_1.1", "length": 11}])[0]["accession"] == "XP_1.1"

    def test_no_protein_to_choose_from_is_refused(self):
        with pytest.raises(ValueError, match="no protein to choose"):
            do.choose_protein([])

    def test_a_gene_whose_transcripts_have_no_protein_has_no_protein(self):
        prod = copy.deepcopy(products("NDUFS2")[100768067])
        for t in prod["transcripts"]:
            t.pop("protein")
        assert do.refseq_proteins(prod, {"NC_048598.1", "NW_003614301.1"}) == []
        f = found("NDUFS2")
        rows = do.gene_rows("NDUFS2", 4720, f, do.chosen_genes(f, TAXA, ASSEMBLY), {**products("NDUFS2"), 100768067: prod}, TAXA, ASSEMBLY)
        assert [r["status"] for r in rows if r["species"] == "cricetulus-griseus"] == ["ortholog has no RefSeq protein on the assembly"]


class TestRows:
    def test_one_row_per_species_and_a_gene_of_each_species_with_its_longest_protein(self):
        rows = rows_of("NDUFS2", 4720)
        assert [r["species"] for r in rows] == list(ASSEMBLY) and {r["status"] for r in rows} == {do.OK}
        by = {r["species"]: r for r in rows}
        assert by["mus-musculus"]["protein"] == "NP_694704.1" and by["mus-musculus"]["gene_id"] == "226646" and by["mus-musculus"]["proteins"] == 3
        assert by["cricetulus-griseus"]["protein"] == "XP_027274810.1" and by["cricetulus-griseus"]["proteins"] == 1     # the other assembly's protein is not counted
        assert all(r["human_gene_id"] == 4720 and r["length"] == 463 for r in rows)

    def test_a_missing_ortholog_is_a_row_with_the_reason_and_no_protein(self):
        by = {r["species"]: r for r in rows_of("NDUFB1", 4707)}
        assert by["cricetulus-griseus"]["status"] == "no ortholog" and by["cricetulus-griseus"]["protein"] == "" and by["cricetulus-griseus"]["gene_id"] == ""
        assert by["mus-musculus"]["status"] == do.OK and by["mus-musculus"]["protein"].startswith(("NP_", "XP_"))

    def test_sequences_are_checked_against_the_length_and_their_digest_is_kept(self):
        t = fixture_table()
        ok = t[t.status == do.OK]
        assert (ok.sha256.str.len() == 64).all() and (t[t.status != do.OK].sha256 == "").all()
        assert all(do.sequence_digest(PROTEINS[r.protein]) == r.sha256 and len(PROTEINS[r.protein]) == int(r.length) for r in ok.itertuples())
        with pytest.raises(ValueError, match="residues"):
            do.fill_sequences(t, fetch=lambda accs: {a: PROTEINS[a][:-1] for a in accs})

    def test_tied_isoforms_are_counted_as_the_same_or_another_sequence(self):
        t = fixture_table()
        r = t[(t.species == "mus-musculus") & (t.gene == "ATP5MC1")].iloc[0]
        assert r.tied_with == "NP_001154891.1;NP_001404256.1" and r.tied_differ == "0"
        alt = dict(PROTEINS)
        alt["NP_001404256.1"] = alt["NP_001404256.1"][:-1] + ("A" if alt["NP_001404256.1"][-1] != "A" else "G")
        r = do.fill_sequences(t, fetch=lambda accs: {a: alt[a] for a in accs})
        r = r[(r.species == "mus-musculus") & (r.gene == "ATP5MC1")].iloc[0]
        assert r.tied_differ == "1"


class TestDownload:
    def table_file(self, tmp_path):
        p = tmp_path / "table.tsv"
        do.write_table(fixture_table(), p)
        return p

    def test_the_table_is_written_in_the_order_of_the_species_and_of_the_genes(self, tmp_path):
        p = tmp_path / "t.tsv"
        do.write_table(fixture_table().sample(frac=1, random_state=3), p)
        t = do.read_table(p)
        order = sorted(zip(t.species, t.gene), key=lambda k: (list(do.SPECIES).index(k[0]), do.genes().index(k[1])))
        assert list(zip(t.species, t.gene)) == order and len(t) == 15

    def test_the_fasta_of_a_species_holds_its_proteins_by_accession_and_the_manifest_the_digests(self, tmp_path):
        p = self.table_file(tmp_path)
        files = do.download(tmp_path / "in", tmp_path / "MANIFEST.json", p, fetch=lambda accs: {a: PROTEINS[a] for a in accs})
        assert sorted(f.name for f in files) == sorted(f"{s}.faa" for s in ASSEMBLY)
        t = do.read_table(p)
        for slug in ASSEMBLY:
            got = do.parse_fasta((tmp_path / "in" / f"{slug}.faa").read_text())
            assert got == {r.protein: PROTEINS[r.protein] for r in t[(t.species == slug) & (t.status == do.OK)].itertuples()}
        rec = json.loads((tmp_path / "MANIFEST.json").read_text())[0]
        import hashlib
        assert rec["sha256"] == {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
        assert rec["table_sha256"] == hashlib.sha256(p.read_bytes()).hexdigest()

    def test_a_sequence_that_is_not_the_pinned_one_refuses_and_writes_nothing(self, tmp_path):
        p = self.table_file(tmp_path)
        acc = "NP_694704.1"
        for bad in (PROTEINS[acc][:-1], PROTEINS[acc][:-1] + ("A" if PROTEINS[acc][-1] != "A" else "G")):      # shorter; the same length, one residue changed
            with pytest.raises(SystemExit, match="not the pinned"):
                do.download(tmp_path / "in", tmp_path / "MANIFEST.json", p, fetch=lambda accs: {a: (bad if a == acc else PROTEINS[a]) for a in accs})
            assert not (tmp_path / "in").exists() and not (tmp_path / "MANIFEST.json").exists()

    def test_the_gaps_are_not_fetched(self, tmp_path):
        p = self.table_file(tmp_path)
        asked = []
        do.download(tmp_path / "in", tmp_path / "MANIFEST.json", p, fetch=lambda accs: (asked.extend(accs), {a: PROTEINS[a] for a in accs})[1])
        t = do.read_table(p)
        assert sorted(asked) == sorted(t[t.status == do.OK].protein) and "" not in asked and not do.gaps(p).empty
        assert list(do.gaps(p).columns) == ["species", "gene", "status", "candidates"] and set(do.gaps(p).status) == {"no ortholog"}

    def test_the_proteins_read_back_are_the_pinned_ones_and_a_changed_file_is_refused(self, tmp_path):
        p = self.table_file(tmp_path)
        do.download(tmp_path / "in", tmp_path / "MANIFEST.json", p, fetch=lambda accs: {a: PROTEINS[a] for a in accs})
        got = do.read(tmp_path / "in", p)
        t = do.read_table(p)
        assert set(got) == {(r.species, r.gene) for r in t[t.status == do.OK].itertuples()}
        o = got[("mus-musculus", "NDUFS2")]
        assert (o.accession, o.gene_id, o.protein) == ("NP_694704.1", "226646", PROTEINS["NP_694704.1"])
        assert ("cricetulus-griseus", "NDUFB1") not in got
        f = tmp_path / "in" / "mus-musculus.faa"
        f.write_text(f.read_text().replace(PROTEINS["NP_694704.1"][:10], "AAAAAAAAAA", 1))
        with pytest.raises(SystemExit, match="not the pinned"):
            do.read(tmp_path / "in", p)
        f.write_text(f.read_text() + ">NP_000001.1 gene=X gene_id=1\nMKV\n")
        with pytest.raises(SystemExit, match="does not hold the table's proteins"):
            do.read(tmp_path / "in", p)
        f.unlink()
        with pytest.raises(SystemExit, match="is missing"):
            do.read(tmp_path / "in", p)

    def test_a_header_twice_in_a_fasta_is_refused_and_wrapped_lines_are_joined(self):
        assert do.parse_fasta(">A x\nmkv\nLA\n>B\nW\n") == {"A": "MKVLA", "B": "W"}
        with pytest.raises(ValueError, match="twice"):
            do.parse_fasta(">A\nM\n>A\nK\n")


class TestSource:
    def source(self, tmp_path):
        p = tmp_path / "table.tsv"
        do.write_table(fixture_table(), p)
        do.download(tmp_path / "in", tmp_path / "MANIFEST.json", p, fetch=lambda accs: {a: PROTEINS[a] for a in accs})
        sp = tmp_path / "species.tsv"
        pd.DataFrame({"species": ["mus-musculus", "cricetulus-griseus", "mus-spretus"], "taxon": [10090, 10029, 10096],
                      "assembly": ["GCF_000001635.27", "GCF_003668045.3", ""], "source": ["RefSeq", "RefSeq", "none"]}).to_csv(sp, sep="\t", index=False)
        return do.read_source(tmp_path / "in", p, sp)

    def test_a_covered_species_is_read_from_the_table_a_gap_and_a_species_without_a_source_are_not_assessed_and_an_uncovered_one_is_left_to_the_alignments(self, tmp_path):
        s = self.source(tmp_path)
        o, reason = s.get("mus-musculus", "NDUFS2")
        assert reason is None and o.accession == "NP_694704.1"
        o, reason = s.get("cricetulus-griseus", "NDUFB1")
        assert o is None and reason == "not assessed: no ortholog of the donor in RefSeq (no ortholog)"
        o, reason = s.get("mus-spretus", "NDUFS2")
        assert o is None and reason.startswith("not assessed: no ortholog source")
        assert s.get("gorilla-gorilla", "NDUFS2") == (None, None)
        with pytest.raises(SystemExit, match="not in the pinned table"):
            s.get("mus-musculus", "COX5B")


class TestPinned:
    def test_the_pinned_table_has_every_gene_of_the_complexes_for_each_species_with_a_source(self):
        t = do.read_table()
        sp = pd.read_csv(do.SPECIES_TABLE, sep="\t", dtype=str, keep_default_na=False)
        have = set(sp.species[sp.assembly != ""])
        assert have == {"mus-musculus", "rattus-norvegicus", "mus-caroli", "mus-pahari", "cricetulus-griseus"}
        assert set(sp.species) == set(do.SPECIES) and set(t.species) == have
        genes = do.genes()
        assert len(genes) == len(set(genes)) == 72 and set(genes) == {g for c in COMPLEXES.values() for g in c.nuclear}
        assert len(t) == len(have) * len(genes) and set(zip(t.species, t.gene)) == {(s, g) for s in have for g in genes}
        assert list(t.columns) == do.COLUMNS

    def test_each_row_is_either_a_protein_with_its_digest_or_a_stated_gap(self):
        t = do.read_table()
        ok = t[t.status == do.OK]
        assert (ok.protein.str.match(r"^(NP|XP)_\d+\.\d+$")).all() and (ok.sha256.str.len() == 64).all() and (ok.length.astype(int) > 0).all()
        assert (ok.gene_id != "").all() and (ok.proteins.astype(int) >= 1).all()
        bad = t[t.status != do.OK]
        assert (bad.protein == "").all() and (bad.sha256 == "").all()
        assert set(bad.status) <= {"no ortholog", "ambiguous (paralogs)", "ortholog is not protein coding", "ortholog not annotated on the assembly",
                                   "ortholog has no RefSeq protein on the assembly"}
        assert not ok.duplicated(["species", "gene_id"]).any()

    def test_a_species_without_a_refseq_genome_is_stated_as_having_no_source(self):
        sp = pd.read_csv(do.SPECIES_TABLE, sep="\t", dtype=str, keep_default_na=False).set_index("species")
        for s in ("mus-spretus", "mus-terricolor", "otomys-irroratus"):
            assert sp.loc[s, "assembly"] == "" and sp.loc[s, "source"].startswith("no ortholog source")

    def test_the_copy_in_the_results_is_the_pin(self):
        r = Path(__file__).parents[1] / "results" / "h36" / "orthologs" / "donor_orthologs.tsv"
        if not r.is_file():
            pytest.skip("results/h36/orthologs/ not in this checkout")
        assert r.read_bytes() == do.TABLE.read_bytes()

    def test_the_downloaded_proteins_agree_with_the_vertebrate_alignments(self):
        d = paths.DONOR_ORTHOLOGS
        aln = paths.MITOPORT / "partners" / "aln"
        if not (d / "mus-musculus.faa").is_file() or not aln.is_dir():
            pytest.skip("needs `mitoage download donor-orthologs` and the vertebrate alignments")
        from mitoage_pipelines.mitoport import cybrid_swap as sw
        src = do.read(d)
        same = total = 0
        for (sp, g), o in src.items():
            f = aln / f"{g}.afa"
            if f.is_file() and sp in (rows := sw.read_alignment(f)):
                total += 1
                same += rows[sp].replace("-", "") == o.protein
        assert total >= 90 and same / total > 0.8              # 90 of 100 identical on 2026-11-24; the rest are another isoform of the same gene
