"""H31's scores from a likelihood table: the drop of a change, a stretch on either strand, the cybrid rule, the
variants read from the two files; and, where `download cybrid` has put the records, the stretches of the human record.
Expected values are worked out by hand in each test."""
import json

import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.mitoport import cybrid, cybrid_dna, kron

GB = paths.CYBRID_GB / "records.gb"
needs_records = pytest.mark.skipif(not GB.is_file(), reason="needs `mitoage download cybrid`")
AL = cybrid_dna.aligner()
GENOME = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"


def table(genome=GENOME):
    """log-probabilities: the genome's own base 0, the others -1, -2, -3 in the order A, C, G, T after leaving it out."""
    lp = np.zeros((len(genome), 4))
    for i, b in enumerate(genome):
        others = [x for x in "ACGT" if x != b]
        for k, x in enumerate(others, 1):
            lp[i, "ACGT".index(x)] = -k
    return lp


class TestDrop:
    def test_the_reference_minus_the_other_and_unknowns(self):
        lp = table()
        assert kron.drop(lp, 0, "A", "A") == 0
        assert kron.drop(lp, 0, "A", "C") == 1 and kron.drop(lp, 0, "A", "T") == 3
        assert np.isnan(kron.drop(lp, 0, "A", "N")) and np.isnan(kron.drop(lp, 0, "N", "A"))
        lp[3] = np.nan
        assert np.isnan(kron.drop(lp, 3, "T", "A"))


class TestTableOf:
    def write(self, tmp_path, genome):
        lp = table(genome)
        t = pd.DataFrame({"pos": range(1, len(genome) + 1), "ref": list(genome), **{f"lp_{b}": lp[:, i] for i, b in enumerate("ACGT")}})
        t.to_csv(tmp_path / "t.tsv", sep="\t", index=False)
        return tmp_path / "t.tsv", lp

    def test_a_table_of_the_genome_and_of_another(self, tmp_path):
        p, lp = self.write(tmp_path, GENOME)
        assert np.array_equal(kron.table_of(p, GENOME, "mito"), lp)
        with pytest.raises(ValueError, match="not of the record"):
            kron.table_of(p, GENOME[1:] + GENOME[0], "mito")

    def test_the_shuffled_arm_takes_another_order_of_the_same_bases_only(self, tmp_path):
        turned = GENOME[1:] + GENOME[0]
        p, lp = self.write(tmp_path, turned)
        got = kron.table_of(p, GENOME, "shuffled")
        assert np.array_equal(got, lp)
        # read at the record's place with the record's bases: place 0 holds A in the record, C in the shuffle
        assert kron.drop(got, 0, "A", "C") == lp[0, 0] - lp[0, 1] == -1.0
        with pytest.raises(ValueError, match="not a shuffle"):
            kron.table_of(p, turned, "shuffled")                      # the genome itself is no shuffle
        with pytest.raises(ValueError, match="not a shuffle"):
            kron.table_of(p, "A" + GENOME[1:-1] + "A", "shuffled")    # other bases


class TestStretchScore:
    def test_plus_strand_two_changes_and_a_gap(self):
        seq = GENOME[4:24]
        host = (seq, np.arange(4, 24), 1)
        donor = seq[:2] + "T" + seq[3:10] + seq[11:15] + "A" + seq[16:]      # G>T at place 6, one base lost, T>A at 19
        assert seq[2] == "G" and seq[15] == "T"
        s, n, gaps = kron.stretch_score(host, donor, GENOME, table(), AL)
        # at a G the others are A -1, C -2, T -3; at a T they are A -1, C -2, G -3
        assert (s, n, gaps) == (3 + 1, 2, 1)

    def test_minus_strand_bases_are_complemented_and_places_run_backwards(self):
        places = np.arange(23, 3, -1)
        seq = "".join(GENOME[i] for i in places).translate(kron.COMPLEMENT)
        host = (seq, places, -1)
        donor = seq[:5] + ("A" if seq[5] != "A" else "C") + seq[6:]
        s, n, gaps = kron.stretch_score(host, donor, GENOME, table(), AL)
        ref, alt = GENOME[places[5]], donor[5].translate(kron.COMPLEMENT)
        assert (s, n, gaps) == (kron.drop(table(), int(places[5]), ref, alt), 1, 0) and s > 0

    def test_a_host_base_that_is_not_the_tables_is_refused(self):
        seq = GENOME[4:24]
        with pytest.raises(ValueError):
            kron.stretch_score((seq, np.arange(5, 25), 1), seq[:2] + "T" + seq[3:], GENOME, table(), AL)

    def test_the_same_stretch_scores_nothing(self):
        seq = GENOME[4:24]
        assert kron.stretch_score((seq, np.arange(4, 24), 1), seq, GENOME, table(), AL) == (0.0, 0, 0)


class TestVerdict:
    def summary(self, t1, t2, t3):
        rows = []
        for s, (a, b, c) in {"kron_cds": (t1, t2, t3), "kron_cds_mean": (-0.1, 0.5, 2), "d_all": (-0.36, 0.74, 4),
                             "f_all": (-0.59, 0.61, 1), "nt_cds_d": (-0.04, 0.83, 4), "nt_cds_f": (-0.63, 0.76, 1)}.items():
            rows += [{"test": "T1", "score": s, "value": a}, {"test": "T2", "score": s, "value": b}]
            rows += [{"test": "T3", "score": s, "value": float(i < c)} for i in range(6)]
        return pd.DataFrame(rows)

    def test_all_three_two_and_fewer(self):
        v = kron.verdict(self.summary(-0.7, 0.9, 5)).set_index("score")
        assert v.loc["kron_cds", "verdict"] == "counts" and v.loc["kron_cds_mean", "verdict"] == "does not count"
        assert (v.loc["kron_cds", "T1_bar"], v.loc["kron_cds", "T2_bar"], v.loc["kron_cds", "T3_bar"]) == (-0.63, 0.83, 4)
        assert kron.verdict(self.summary(-0.7, 0.9, 4)).set_index("score").loc["kron_cds", "verdict"].startswith("partial")
        assert kron.verdict(self.summary(-0.63, 0.83, 5)).set_index("score").loc["kron_cds", "verdict"] == "does not count"


class TestReadVariants:
    def test_labels_and_what_is_left_out(self, tmp_path):
        clinvar = [
            {"type": "single nucleotide variant", "classification": "Pathogenic", "review": kron.variant_stab.PANEL, "spdi": "NC_012920.1:9:A:G"},
            {"type": "single nucleotide variant", "classification": "Likely pathogenic", "review": "criteria provided, single submitter", "spdi": "NC_012920.1:19:C:T"},
            {"type": "single nucleotide variant", "classification": "Benign", "review": kron.variant_stab.PANEL, "spdi": "NC_012920.1:29:C:T"},
            {"type": "Deletion", "classification": "Pathogenic", "review": kron.variant_stab.PANEL, "spdi": "NC_012920.1:39:CC:C"},
            {"type": "single nucleotide variant", "classification": "Pathogenic", "review": kron.variant_stab.PANEL, "spdi": None},
        ]
        (tmp_path / "c.json").write_text(json.dumps(clinvar))
        pd.DataFrame({"position": [10, 20, 30, 40, 50, 60], "ref": list("ACCGTA"), "alt": ["G", "T", "T", "A", "C", "AT"],
                      "filters": ["PASS"] * 4 + ["npg", "PASS"], "AF_hom": [0.5, 0.5, 0.002, 0.0001, 0.5, 0.5]}).to_csv(tmp_path / "g.tsv", sep="\t", index=False)
        v, counts = kron.read_variants(tmp_path / "c.json", tmp_path / "g.tsv")
        assert [(r.pos, r.ref, r.alt, r.label) for r in v.itertuples()] == [(10, "A", "G", "harmful"), (30, "C", "T", "harmless")]
        assert counts == {"harmless left out: pathogenic in ClinVar": 2, "harmful without a sequence change given": 1}


@needs_records
class TestOnTheRecords:
    def test_the_stretches_of_the_human_record(self):
        st = kron.stretches(GB, cybrid.GENOMES[cybrid.HUMAN])
        from Bio import SeqIO
        rec = next(r for r in SeqIO.parse(str(GB), "genbank") if r.id == "NC_012920.1")
        genome = str(rec.seq).upper()
        assert len([k for k in st if k[0] == "cds"]) == 13 and len([k for k in st if k[0] == "trna"]) == 22
        seq, places, strand = st[("cds", "ND6")]
        assert strand == -1 and "".join(genome[i] for i in places).translate(kron.COMPLEMENT) == seq and seq.startswith("ATG")
        seq, places, strand = st[("cds", "COX1")]
        assert strand == 1 and "".join(genome[i] for i in places) == seq
        seq, places, _ = st[("dloop",)]
        assert places[0] == 16023 and places[-1] == 575 and "".join(genome[i] for i in places) == seq

    def test_the_chinese_hamsters_listed_corrections(self):
        # RECORD_FIXES: two tRNA genes on the minus strand, read 5' to 3' with their places running backwards, and
        # the control region as the stretch after tRNA-Pro
        s = kron.stretches(GB, cybrid.GENOMES["cricetulus-griseus"])
        for product in ("tRNA-Asn", "tRNA-Cys"):
            seq, places, strand = s[("trna", (product, 0))]
            assert strand == -1 and list(places) == sorted(places, reverse=True) and len(seq) == len(places)
        seq, places, strand = s[("dloop",)]
        assert (int(places[0]), int(places[-1]), len(seq), strand) == (15403, 16283, 881, 1)

    def test_a_trna_gene_on_different_strands_is_refused(self, monkeypatch):
        real = kron.stretches

        def crossed(gb, accession):
            out = dict(real(gb, accession))
            if accession == cybrid.GENOMES["mus-spretus"]:
                k = ("trna", ("tRNA-Asn", 0)); seq, places, strand = out[k]
                out[k] = (seq, places, -strand)
            return out
        monkeypatch.setattr(kron, "stretches", crossed)
        genome = next(str(r.seq).upper() for r in __import__("Bio.SeqIO", fromlist=["x"]).parse(str(GB), "genbank") if r.id == cybrid.GENOMES["mus-musculus"])
        with pytest.raises(ValueError, match="different strands"):
            kron.pair_rows(GB, "mus-musculus", "mus-spretus", genome, np.zeros((len(genome), 4)))

    def test_the_kinds_of_place(self):
        seq, kind, _ = kron.kinds(GB)
        assert kind[3306] == "protein" and kind[100] == "control region" and kind[700] == "rRNA" and kind[3229] == "tRNA"
        assert set(kind) <= set(kron.KINDS)

    def test_a_pair_against_a_flat_table_counts_its_scored_columns(self):
        genome = kron.stretches(GB, "NC_005089.1")
        from Bio import SeqIO
        rec = next(r for r in SeqIO.parse(str(GB), "genbank") if r.id == "NC_005089.1")
        g = str(rec.seq).upper()
        lp = np.full((len(g), 4), -1.0)
        lp[np.arange(len(g)), [("ACGT".index(b) if b in "ACGT" else 0) for b in g]] = 0.0      # every change drops by 1
        rows = kron.pair_rows(GB, "mus-musculus", "mus-spretus", g, lp)
        dna = {r["complex"]: r for r in cybrid_dna.pair_rows(cybrid_dna.read_features(GB), "mus-musculus", "mus-spretus")}
        for r in rows:
            assert r["kron_cds"] == r["kron_cds_n"] and r["kron_cds_mean"] == 1
            assert r["kron_cds_n"] + r["kron_cds_gaps"] == dna[r["complex"]]["nt_cds_d"]
        same = kron.pair_rows(GB, "mus-musculus", "mus-musculus", g, lp)
        assert all(r["kron_cds"] == 0 and r["kron_rrna"] == 0 and r["kron_trna"] == 0 and r["kron_dloop"] == 0 for r in same)
