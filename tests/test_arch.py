"""H38: genome architecture against lifespan (`longevity/arch_measures.py`, `longevity/arch_tests.py`).

The records are real: the mitogenomes of the cybrid hosts and donors in `data/input/cybrid/records.gb` (`mitoage download
cybrid`: human, apes, rat, mouse, hamster; the human record's control region wraps the origin, the hamster's has none) and, for the
tests that need the build, the tables the job makes before the check (`audit records|strand|order|code`, `arch-measures`,
`longevity controls`, `dated-tree`). Every expected value is computed from the record's own text or from the tables at run time by a
path other than the code's. The bar is checked case by case on tables of exactly the shape `run_cells` writes: each case changes the
cells of one lambda model or one set at a time, so that the conditions of the bar are told apart from each other.
"""
import re
import shutil

import numpy as np
import pandas as pd
import pytest
from statsmodels.stats.multitest import multipletests

from mitoage_pipelines import paths
from mitoage_pipelines.audit import order as aorder
from mitoage_pipelines.audit import records as arec
from mitoage_pipelines.build2 import genome as gm
from mitoage_pipelines.longevity import aa_pgls
from mitoage_pipelines.longevity import arch_measures as am
from mitoage_pipelines.longevity import arch_tests as at
from mitoage_pipelines.longevity import rna_tests

RECORDS = paths.CYBRID_GB / "records.gb"
needs_records = pytest.mark.skipif(not RECORDS.is_file(), reason="needs `mitoage download cybrid`")
needs_trf = pytest.mark.skipif(shutil.which("trf") is None, reason="needs the pinned Tandem Repeats Finder on PATH")
COMP = paths.LONGEVITY_OUT / "composition.tsv"
ARCH = paths.ARCH_OUT / "architecture.tsv"
needs_build = pytest.mark.skipif(
    not all(p.is_file() for p in (paths.BUILD / "species.parquet", COMP, paths.DATED_TREE, ARCH, paths.AUDIT_OUT / "order.tsv",
                                  paths.AUDIT_OUT / "order_flags.tsv", paths.AUDIT_OUT / "strand.tsv", paths.AUDIT_OUT / "code.tsv",
                                  paths.AUDIT_OUT / "annotation_vs_build.tsv")),
    reason="needs the build, `longevity controls`, `dated-tree`, the audit's tables and `longevity arch-measures`")
HUMAN, HAMSTER = "NC_012920", "NC_007936"
NAMES = ["NC_012920", "NC_001643", "NC_001644", "NC_011120", "NC_001665", "NC_005089", "NC_007936"]


@pytest.fixture(scope="module")
def records() -> dict[str, str]:
    texts = [t.lstrip() + "\n//\n" for t in RECORDS.read_text().split("\n//\n") if t.strip()]
    out = {re.search(r"^VERSION\s+(\S+)\.\d+", t, re.M).group(1): t for t in texts}
    return {k: out[k] for k in NAMES}


def blocks(text: str) -> tuple[str, list[str], str]:
    head, rest = text.split("\nFEATURES", 1)
    feat, tail = rest.split("\nORIGIN", 1)
    lines = ("FEATURES" + feat).split("\n")
    out, cur = [], []
    for line in lines[1:]:
        if re.match(r"^ {5}\S", line) and cur:
            out.append("\n".join(cur)); cur = []
        cur.append(line)
    out.append("\n".join(cur))
    return head + "\n" + lines[0], out, "\nORIGIN" + tail


def join(head: str, feats: list[str], tail: str) -> str:
    return head + "\n" + "\n".join(feats) + tail


def kind_of(block: str) -> str:
    return re.match(r"^ {5}(\S+)", block).group(1)


def genome(text: str, phylum: str = "Chordata") -> gm.Genome:
    return gm.annotate(text, gm.fallback_table(phylum), phylum == "Chordata")


def arrangement(text: str) -> dict:
    return aorder.arrangement(genome(text))


def without(text: str, keep) -> str:
    head, feats, tail = blocks(text)
    return join(head, [b for b in feats if keep(b)], tail)


def positions(loc: str) -> set[int]:
    """The positions a GenBank location covers, read from its numbers alone (every a..b of a join), for the independent counts."""
    out: set[int] = set()
    for a, b in re.findall(r"(\d+)\.\.(\d+)", loc):
        out |= set(range(int(a), int(b) + 1))
    return out


def location(block: str) -> str:
    return re.sub(r"\s+", "", re.match(r"^ {5}\S+\s+((?:[^/]|\n {21}[^/])*)", block).group(1))


# --- the distance ----------------------------------------------------------------------------------
def naive_distance(a: str, b: str) -> float:
    """The breakpoint distance by its textbook definition for two circular signed orders with the same unique genes: the adjacencies
    of a that b has in neither reading. Written with ordered pairs of signed genes and no multiset."""
    def pairs(o):
        t = o.split(" ")
        return [(t[i], t[(i + 1) % len(t)]) for i in range(len(t))]
    neg = lambda x: ("-" if x[0] == "+" else "+") + x[1:]
    pb = set(pairs(b))
    return float(sum(1 for (x, y) in pairs(a) if (x, y) not in pb and (neg(y), neg(x)) not in pb))


class TestDistance:
    @needs_records
    def test_it_is_the_textbook_breakpoint_distance_where_the_genes_are_the_same(self, records):
        for tier in ("order_full_unoriented", "order_genes_unoriented"):
            o = [arrangement(t)[tier] for t in records.values()]
            t = arrangement(records[HUMAN])[tier].split(" ")
            n = len(t)
            variants = [" ".join(t),
                        " ".join(t[:3] + [t[4], t[3]] + t[5:]),                                                        # a swap
                        " ".join(t[:2] + aorder.flip(t[2:7]) + t[7:]),                                                 # a segment turned over
                        " ".join(t[:2] + t[5:9] + t[2:5] + t[9:]),                                                     # a block moved
                        " ".join(aorder.flip(t[:n // 2]) + t[n // 2:]),                                                # half the circle turned over
                        *o]
            for a in variants:
                for b in variants:
                    assert am.breakpoint_distance(a, b) == naive_distance(a, b) == am.breakpoint_distance(b, a)
            assert max(am.breakpoint_distance(variants[0], v) for v in variants) >= 3

    @needs_records
    def test_a_swap_of_two_neighbouring_genes_is_three_breakpoints_and_a_move_back_is_zero(self, records):
        o = arrangement(records[HUMAN])["order_full_unoriented"].split(" ")
        i = 5
        sw = o[:i] + [o[i + 1], o[i]] + o[i + 2:]
        assert am.breakpoint_distance(" ".join(o), " ".join(sw)) == 3 == naive_distance(" ".join(o), " ".join(sw))
        assert am.breakpoint_distance(" ".join(o), " ".join(o)) == 0

    @needs_records
    def test_it_is_unmoved_by_the_origin_and_by_the_strand_the_record_is_read_from(self, records):
        o = arrangement(records[HUMAN])["order_full"]
        t = o.split(" ")
        rotated = " ".join(t[7:] + t[:7])
        flipped = " ".join(aorder.flip(t))
        assert am.breakpoint_distance(o, rotated) == 0 and am.breakpoint_distance(o, flipped) == 0
        assert am.breakpoint_distance(o, " ".join(aorder.flip(t[7:] + t[:7]))) == 0

    @needs_records
    def test_a_record_deposited_from_the_other_strand_is_at_zero_and_a_real_rearrangement_is_not(self, records):
        g = genome(records[HUMAN])
        rc = arec.reverse_complement(g)
        assert aorder.arrangement(rc)["order_full"] != aorder.arrangement(g)["order_full"]
        assert am.breakpoint_distance(aorder.arrangement(rc)["order_full_unoriented"], aorder.arrangement(g)["order_full_unoriented"]) == 0
        mouse = arrangement(records["NC_005089"])["order_full_unoriented"]
        human = arrangement(records[HUMAN])["order_full_unoriented"]
        assert am.breakpoint_distance(mouse, human) == naive_distance(mouse, human)

    @needs_records
    def test_a_gene_present_in_one_record_only_counts_on_its_side_and_the_distance_stays_symmetric(self, records):
        full = arrangement(records[HUMAN])
        cut = arrangement(without(records[HUMAN], lambda b: not (kind_of(b) == "tRNA" and "Trp" in b)))
        assert "tRNA-Trp" not in cut["order_full"] and "tRNA-Trp" in full["order_full"]
        d = am.breakpoint_distance(full["order_full_unoriented"], cut["order_full_unoriented"])
        assert d == am.breakpoint_distance(cut["order_full_unoriented"], full["order_full_unoriented"]) > 0
        assert d == 1.5                      # the two adjacencies of the lost gene go, the one that closes the gap comes: |{3 changed}| / 2
        assert am.breakpoint_distance(full["order_genes_unoriented"], cut["order_genes_unoriented"]) == 0

    @needs_records
    def test_a_record_with_no_trna_annotated_is_unknown_at_the_full_tier_never_zero(self, records):
        text = without(records[HUMAN], lambda b: kind_of(b) != "tRNA")
        a = arrangement(text)
        assert a["order_full"] is None and a["order_full_unoriented"] is None and a["trnas_annotated"] == 0
        assert a["order_genes_unoriented"] == arrangement(records[HUMAN])["order_genes_unoriented"]
        assert am.breakpoint_distance(a["order_full_unoriented"], arrangement(records[HUMAN])["order_full_unoriented"]) is None
        assert am.breakpoint_distance(None, None) is None and am.breakpoint_distance(float("nan"), "+COX1") is None
        tab = pd.DataFrame({"class": ["C"] * 3, "order_full_unoriented": [a["order_full_unoriented"], "+A +B +C", "+A +B +C"],
                            "order_genes_unoriented": [a["order_genes_unoriented"]] * 3})
        d = am.class_distances(tab, "full")
        assert pd.isna(d.rearr_full.iloc[0]) and list(d.rearr_full.iloc[1:]) == [0, 0]
        assert d.rearr_full.dtype == "Float64"

    def test_the_class_modal_arrangement_is_the_commonest_and_the_first_in_string_order_on_a_tie(self):
        assert am.modal_arrangement(pd.Series(["+B +C", "+A +C", "+B +C", None])) == "+B +C"
        assert am.modal_arrangement(pd.Series(["+B +C", "+A +C", "+B +C", "+A +C"])) == "+A +C"
        assert am.modal_arrangement(pd.Series([None, None], dtype=object)) is None

    @needs_records
    def test_a_genome_with_the_class_consensus_is_at_zero_and_one_record_away_from_it_by_its_distance(self, records):
        human = arrangement(records[HUMAN])["order_full_unoriented"]
        t = human.split(" ")
        moved = " ".join(t[:3] + [t[4], t[3]] + t[5:])
        tab = pd.DataFrame({"class": ["M", "M", "M", "M", None, "N"], "order_full_unoriented": [human, human, human, moved, human, moved]})
        d = am.class_distances(tab, "full")
        assert list(d.rearr_full.iloc[:3]) == [0, 0, 0] and d.rearr_full.iloc[3] == naive_distance(moved, human) == 3
        assert pd.isna(d.rearr_full.iloc[4])                                # no class: unknown, not zero
        assert d.rearr_full.iloc[5] == 0                                    # a class of its own record is its own consensus
        assert list(d.modal_full.iloc[:4]) == [human] * 4

    @needs_build
    def test_on_the_build_the_distance_is_zero_exactly_where_the_arrangement_is_the_classes_commonest(self):
        a = pd.read_csv(ARCH, sep="\t")
        for tier in am.TIERS:
            col, dist = f"order_{tier}_unoriented", f"rearr_{tier}"
            known = a[a[col].notna() & a["class"].notna()]
            modal = known.groupby("class")[col].agg(lambda s: min(s.value_counts()[lambda v: v == v.max()].index))
            assert ((known[dist] == 0) == (known[col] == known["class"].map(modal))).all()
            assert (known[dist] >= 0).all() and a.loc[a[col].isna(), dist].isna().all()
        assert a.loc[a.trnas_annotated == 0, "rearr_full"].isna().all() and a.loc[a.trnas_annotated == 0, "rearr_genes"].notna().all()


# --- the control region and the overlaps ----------------------------------------------------------------
class TestRegions:
    @needs_records
    def test_the_control_region_that_wraps_the_origin_is_its_two_pieces(self, records):
        text = records[HUMAN]
        dloop = [b for b in blocks(text)[1] if kind_of(b) == "D-loop"]
        assert len(dloop) == 1 and "join" in location(dloop[0])
        n_expected = len(positions(location(dloop[0])))
        assert am.control_region(genome(text)) == (1, n_expected)
        assert n_expected == 1122

    @needs_records
    def test_a_record_with_two_control_regions_has_two_and_the_sum_of_their_lengths(self, records):
        text = records["NC_001643"]
        head, feats, tail = blocks(text)
        d = [b for b in feats if kind_of(b) == "D-loop"][0]
        nxt = re.sub(r"\d+\.\.\d+", "100..299", d, count=1)
        assert nxt != d and len(positions(location(nxt))) == 200
        two = join(head, feats + [nxt], tail)
        copies, length = am.control_region(genome(two))
        assert copies == 2 and length == len(positions(location(d))) + 200
        assert am.control_region(genome(text)) == (1, len(positions(location(d))))

    @needs_records
    def test_a_record_with_no_control_region_is_unknown_not_zero_and_keeps_its_count(self, records):
        for text in (records[HAMSTER], without(records[HUMAN], lambda b: kind_of(b) != "D-loop")):
            assert am.control_region(genome(text)) == (0, None)
            row = am.measure_record({"slug": "x", "genbank": text, "phylum": "Chordata"})
            assert row["cr_annotated"] == 0 and row["cr_copies"] is None and row["cr_length"] is None and row["cr_bases"] == 0

    @needs_records
    def test_gene_overlap_is_the_positions_in_two_or_more_genes(self, records):
        text = records[HUMAN]
        feats = [b for b in blocks(text)[1] if kind_of(b) in ("CDS", "rRNA", "tRNA")]
        depth: dict[int, int] = {}
        for b in feats:
            for p in positions(location(b)):
                depth[p] = depth.get(p, 0) + 1
        want = sum(1 for v in depth.values() if v >= 2)
        assert want > 0 and am.gene_overlap(genome(text)) == want

    @needs_records
    def test_a_gene_annotated_twice_in_the_same_place_is_not_an_overlap_and_one_in_a_shifted_place_is(self, records):
        text = records[HUMAN]
        head, feats, tail = blocks(text)
        base = am.gene_overlap(genome(text))
        cox1 = [b for b in feats if kind_of(b) == "CDS" and "COX1" in b][0]
        assert am.gene_overlap(genome(join(head, feats + [cox1], tail))) == base
        a, b = map(int, re.findall(r"\d+", location(cox1))[:2])
        shifted = re.sub(r"\d+\.\.\d+", f"{a + 10}..{b + 10}", cox1, count=1)
        sh = genome(join(head, feats + [shifted], tail))
        assert am.gene_overlap(sh) > base

    @needs_records
    def test_a_control_region_laid_over_a_gene_is_not_a_gene_overlap(self, records):
        text = records["NC_001644"]
        head, feats, tail = blocks(text)
        cox1 = [b for b in feats if kind_of(b) == "CDS" and "COX1" in b][0]
        a, b = map(int, re.findall(r"\d+", location(cox1))[:2])
        d = [x for x in feats if kind_of(x) == "D-loop"][0]
        over = re.sub(r"\d+\.\.\d+", f"{a}..{b}", d, count=1)
        assert am.gene_overlap(genome(join(head, feats + [over], tail))) == am.gene_overlap(genome(text))

    @needs_records
    def test_measure_record_gives_the_audit_arrangement_and_the_independent_counts(self, records):
        text = records["NC_001644"]
        row = am.measure_record({"slug": "x", "genbank": text, "phylum": "Chordata"})
        a = arrangement(text)
        assert row["order_full_unoriented"] == a["order_full_unoriented"] and row["order_genes_unoriented"] == a["order_genes_unoriented"]
        assert row["genome_length"] == len(gm.parse_genbank(text).sequence)
        d = [b for b in blocks(text)[1] if kind_of(b) == "D-loop"][0]
        assert row["cr_length"] == len(positions(location(d))) == row["cr_bases"] and row["cr_copies"] == 1


# --- the repeats --------------------------------------------------------------------------------------
DAT = """Tandem Repeats Finder Program written by:

Gary Benson
Program in Bioinformatics
Boston University
Version 4.09


Sequence: a



Parameters: 2 7 7 80 10 50 500


301 420 10 12.0 10 100 0 240 30 30 20 20 1.97 ACGTTGCAAC ACGTTGCAACACGTTGCAAC
500 530 2 15.5 2 100 0 62 40 0 60 0 1.00 AG AGAGAG

Sequence: b



Parameters: 2 7 7 80 10 50 500

"""


class TestRepeats:
    def test_the_dat_file_is_read_for_every_sequence_and_the_parameters_are_checked(self):
        got = am.parse_dat(DAT, ["a", "b"])
        assert [r["start"] for r in got["a"]] == [301, 500] and got["a"][0]["period"] == 10 and got["a"][0]["consensus"] == "ACGTTGCAAC"
        assert got["b"] == []
        with pytest.raises(RuntimeError, match="no block"):
            am.parse_dat(DAT, ["a", "b", "c"])
        with pytest.raises(RuntimeError, match="parameters"):
            am.parse_dat(DAT.replace("2 7 7 80 10 50 500", "2 7 7 80 10 40 500"), ["a", "b"])

    def test_the_shares_are_the_union_of_the_intervals_and_the_part_inside_the_control_region(self):
        n = 1000
        reps = [{"start": 301, "end": 420}, {"start": 400, "end": 450}, {"start": 990, "end": 1010}]
        cr = np.zeros(n, dtype=bool)
        cr[349:600] = True                              # positions 350..600
        whole, in_cr = am.repeat_shares(reps, n, cr)
        assert whole == pytest.approx((150 + 11) / n)                       # 301..450 and 990..1000: the overlap is counted once
        assert in_cr == pytest.approx(len(set(range(301, 451)) & set(range(350, 601))) / 251)
        assert am.repeat_shares(reps, n, None)[1] is None and am.repeat_shares(reps, n, np.zeros(n, dtype=bool))[1] is None
        assert am.repeat_shares([], n, cr) == (0.0, 0.0)

    @needs_records
    @needs_trf
    def test_tandem_repeats_finder_finds_a_repeat_put_into_a_real_control_region_and_nothing_where_there_is_none(self, records, tmp_path):
        g = genome(records["NC_001644"])
        seq = g.record.sequence.upper()
        mask = am.cr_mask(g)
        k = int(np.where(mask)[0][100])
        rep = "ACGTTGCAAC" * 12
        edited = seq[:k] + rep + seq[k:]
        cr2 = np.concatenate([mask[:k], np.ones(len(rep), dtype=bool), mask[k:]])
        binary = am.find_trf()
        got = am.run_trf(binary, {"plain": seq, "edited": edited}, tmp_path)
        hit = [r for r in got["edited"] if r["start"] <= k + 10 and r["end"] >= k + len(rep) - 10]
        assert hit, got["edited"]
        whole, in_cr = am.repeat_shares(got["edited"], len(edited), cr2)
        cover = set()
        for r in got["edited"]:
            cover |= set(range(r["start"], r["end"] + 1))
        assert whole == pytest.approx(len(cover) / len(edited))
        assert in_cr == pytest.approx(len(cover & set(np.where(cr2)[0] + 1)) / cr2.sum())
        assert in_cr > am.repeat_shares(got["plain"], len(seq), mask)[1] and len(cover) >= len(rep) - 5

    @needs_trf
    def test_a_sequence_without_repeats_gives_none_and_the_binary_is_the_pinned_one(self, tmp_path):
        rng = np.random.default_rng(5)
        s = "".join(rng.choice(list("ACGT"), 400))
        assert am.run_trf(am.find_trf(), {"r": s}, tmp_path)["r"] == []
        assert am.sha256_of(shutil.which("trf")) == am.TRF_SHA256

    def test_another_binary_is_refused(self, tmp_path, monkeypatch):
        fake = tmp_path / "trf"
        fake.write_text("#!/bin/sh\nexit 0\n")
        fake.chmod(0o755)
        monkeypatch.setenv("PATH", f"{tmp_path}:/usr/bin")
        with pytest.raises(SystemExit, match="not the pinned"):
            am.find_trf()


# --- the annotation flag and the control on the build's own tables -------------------------------------
@needs_build
class TestFlagsAndControl:
    def test_a_species_is_flagged_by_each_of_the_four_audit_reasons_and_by_no_other(self):
        f = at.annotation_flags(paths.AUDIT_OUT)
        o = pd.read_csv(paths.AUDIT_OUT / "order.tsv", sep="\t")
        assert set(f.slug) == set(o.slug) and f.slug.is_unique
        flags = pd.read_csv(paths.AUDIT_OUT / "order_flags.tsv", sep="\t")
        strand = pd.read_csv(paths.AUDIT_OUT / "strand.tsv", sep="\t")
        code = pd.read_csv(paths.AUDIT_OUT / "code.tsv", sep="\t")
        build = pd.read_csv(paths.AUDIT_OUT / "annotation_vs_build.tsv", sep="\t", keep_default_na=False)
        assert set(f.slug[f.flag_order]) == set(flags.slug)
        assert set(f.slug[f.flag_strand]) == set(strand.slug[strand.orientation.isin(["reversed", "mixed"])]) and f.flag_strand.any()
        assert set(f.slug[f.flag_code]) == set(code.slug[code.neither_clean == "yes"]) and f.flag_code.any()
        assert set(f.slug[f.flag_build]) == set(build.slug[build.build_annotation == "differs"])
        assert ((f.annot_flag == 1) == f[at.FLAG_COLUMNS].any(axis=1)).all()
        for c in at.FLAG_COLUMNS:                                        # each reason alone flags a record the others do not
            only = f[f[c] & ~f[[x for x in at.FLAG_COLUMNS if x != c]].any(axis=1)]
            assert (only.annot_flag == 1).all()
        assert (f.annot_flag == 0).sum() > 1000

    def test_a_species_missing_from_one_audit_table_is_refused_not_taken_as_unflagged(self, tmp_path):
        for name in ("order.tsv", "order_flags.tsv", "strand.tsv", "code.tsv", "annotation_vs_build.tsv"):
            shutil.copy(paths.AUDIT_OUT / name, tmp_path / name)
        s = pd.read_csv(tmp_path / "strand.tsv", sep="\t")
        s.iloc[1:].to_csv(tmp_path / "strand.tsv", sep="\t", index=False)
        with pytest.raises(ValueError, match="same species"):
            at.annotation_flags(tmp_path)

    def tables(self):
        order = pd.read_csv(paths.AUDIT_OUT / "order.tsv", sep="\t")
        flags = pd.read_csv(paths.AUDIT_OUT / "order_flags.tsv", sep="\t")
        arch = pd.read_csv(ARCH, sep="\t")
        return arch, order, flags

    def test_the_control_holds_on_the_audits_own_tables(self):
        arch, order, flags = self.tables()
        tab, depart, counts = at.control_tables(arch, order, flags)
        assert (tab.exceptions == 0).all() and len(tab) == 9
        assert set(tab.check.str[:2]) == {"G0", "G1", "G2", "G3", "G4"}
        assert (tab[tab.check.str.startswith("G1")].records > 0).all() and (tab[tab.check.str.startswith("G3")].records > 100).all()
        assert counts["full"]["uniform_families"] > 0 and not depart.empty
        for tier in am.TIERS:                    # the families the audit would call uniform, counted from the audit's table alone
            known = order[order[f"order_{tier}"].notna() & order.family.notna()]
            uni = known.groupby("family").filter(lambda d: len(d) >= 3 and d[f"order_{tier}"].nunique() == 1)
            row = tab[tab.check.str.startswith("G3") & (tab.tier == tier)].iloc[0]
            assert row.records == len(uni) and row.families == uni.family.nunique()
            assert counts[tier]["flagged_records"] == int((flags.tier == tier).sum())

    def test_the_arrangements_of_the_measure_must_be_the_audits_g0(self):
        arch, order, flags = self.tables()
        a = arch.copy()
        i = a.index[a.order_full_unoriented.notna()][0]
        a.loc[i, "order_full_unoriented"] = a.loc[i, "order_full_unoriented"] + " +X"
        tab = at.control_tables(a, order, flags)[0].set_index("check")
        assert tab.loc["G0 the arrangements are the audit's", "exceptions"] == 1
        b = arch.iloc[1:]                                                   # a live species the measure lacks is fine; one the audit lacks is not
        assert at.control_tables(b, order, flags)[0].exceptions.eq(0).all()
        c = order[order.slug != arch.slug.iloc[0]]
        assert at.control_tables(arch, c, flags)[0].set_index("check").loc["G0 the arrangements are the audit's", "exceptions"] >= 1

    def test_a_flag_on_a_record_that_sits_at_its_familys_arrangement_fails_g1(self):
        arch, order, flags = self.tables()
        uniform = order[order.status.eq("full")].groupby("family").filter(lambda d: len(d) >= 3 and d.order_full.nunique() == 1).iloc[0]
        fake = {**flags.iloc[0].to_dict(), "tier": "full", "slug": uniform.slug, "kind": "order"}
        tab = at.control_tables(arch, order, pd.concat([flags, pd.DataFrame([fake])], ignore_index=True))[0]
        row = tab[tab.check.str.startswith("G1") & (tab.tier == "full")].iloc[0]
        assert row.exceptions == 1 and (tab.exceptions.sum() == 1)

    def test_a_reversed_deposit_that_is_above_zero_fails_g2_and_one_that_is_a_rearrangement_passes_g1(self):
        arch, order, flags = self.tables()
        flags = flags.copy()
        i = flags.index[(flags.tier == "full") & (flags.kind != "reversed deposit")][0]
        was = flags.at[i, "kind"]
        flags.at[i, "kind"] = "reversed deposit"
        tab = at.control_tables(arch, order, flags)[0].set_index(["check", "tier"])
        assert tab.loc[("G2 a reversed deposit is at zero", "full"), "exceptions"] == 1
        assert was != "reversed deposit"

    def test_a_record_of_a_uniform_family_that_the_distance_puts_away_from_it_fails_g3(self):
        arch, order, flags = self.tables()
        o = order.copy()
        uni = o[o.status.eq("full")].groupby("family").filter(lambda d: len(d) >= 3 and d.order_full.nunique() == 1)
        i = uni.index[0]
        o.loc[i, "order_full_unoriented"] = o.loc[i, "order_full_unoriented"] + " +X"
        tab = at.control_tables(arch, o, flags)[0].set_index(["check", "tier"])
        assert tab.loc[("G3 a record of a uniform family is at zero from its family", "full"), "exceptions"] == 1

    def test_a_class_distance_that_is_not_the_family_distance_where_the_two_references_agree_fails_g4(self):
        arch, order, flags = self.tables()
        a = arch.copy()
        uni = order[order.status.eq("full")].groupby("family").filter(lambda d: len(d) >= 3 and d.order_full.nunique() == 1)
        same_as_class = a[a.slug.isin(uni.slug) & (a.rearr_full == 0)].index[0]
        a.loc[same_as_class, "rearr_full"] = 2.0
        tab = at.control_tables(a, order, flags)[0].set_index(["check", "tier"])
        assert tab.loc[("G4 where the family's arrangement is its class's, class distance equals family distance", "full"), "exceptions"] == 1

    def test_a_family_below_the_audits_minimum_is_neither_flagged_nor_uniform_and_adds_no_exception(self):
        arch, order, flags = self.tables()
        o = order[order.status.eq("full")].copy()
        base = at.control_tables(arch, order, flags)[0].set_index(["check", "tier"])
        free = o[~o.slug.isin(flags.slug)]
        first = free.iloc[:1]
        two = pd.concat([first, free[free.order_full.ne(first.order_full.iloc[0])].iloc[:1]]).slug
        o2 = order.copy()
        o2.loc[o2.slug.isin(two), "family"] = "Pairidae"
        tab = at.control_tables(arch, o2, flags)[0].set_index(["check", "tier"])
        key = ("G3 a record of a uniform family is at zero from its family", "full")
        assert tab.loc[key, "exceptions"] == 0 and (tab.exceptions == 0).all()
        assert tab.loc[key, "records"] <= base.loc[key, "records"]

    def test_the_records_of_families_unlike_their_class_are_listed_with_both_distances(self):
        arch, order, flags = self.tables()
        _, depart, counts = at.control_tables(arch, order, flags)
        up = depart[(depart["in"] == "uniform family") & (depart.tier == "full")]
        assert (up.family_distance == 0).all() and (up.class_distance > 0).all()
        assert counts["full"]["uniform_family_records_above_zero_from_class"] == len(up) and len(up) > 0


# --- the gate ------------------------------------------------------------------------------------------
class TestGate:
    def species_file(self, tmp_path):
        p = tmp_path / "species.parquet"
        p.write_bytes(b"a build")
        return p

    def test_no_test_is_run_without_the_control_or_with_a_failed_one_or_one_of_another_build(self, tmp_path):
        sp = self.species_file(tmp_path)
        with pytest.raises(SystemExit, match="arch-control"):
            at.require_control(tmp_path, tmp_path)
        import json
        (tmp_path / at.CONTROL_FILE).write_text(json.dumps({"reproduced": False, "species_sha256": am.sha256_of(sp), "checks_with_exceptions": 2}))
        with pytest.raises(SystemExit, match="did not reproduce"):
            at.require_control(tmp_path, tmp_path)
        (tmp_path / at.CONTROL_FILE).write_text(json.dumps({"reproduced": True, "species_sha256": "0" * 64}))
        with pytest.raises(SystemExit, match="another build"):
            at.require_control(tmp_path, tmp_path)
        (tmp_path / at.CONTROL_FILE).write_text(json.dumps({"reproduced": True, "species_sha256": am.sha256_of(sp)}))
        assert at.require_control(tmp_path, tmp_path)["reproduced"] is True
        (tmp_path / at.CONTROL_FILE).write_text(json.dumps({"reproduced": "yes", "species_sha256": am.sha256_of(sp)}))
        with pytest.raises(SystemExit, match="did not reproduce"):
            at.require_control(tmp_path, tmp_path)

    def test_the_tests_command_is_refused_without_the_control_and_the_control_command_exits_non_zero_when_it_fails(self, tmp_path):
        from typer.testing import CliRunner
        from mitoage_pipelines.cli import app
        data = tmp_path / "data"
        data.mkdir()
        (data / "species.parquet").write_bytes(b"a build")
        run = CliRunner().invoke(app, ["longevity", "arch-tests", "--build", str(data), "--out", str(tmp_path)])
        assert run.exit_code == 1 and "arch-control" in str(run.exception) and not (tmp_path / "arch_cells.tsv").exists()
        rows = ["slug\tclass\tfamily\torder_full_unoriented\torder_genes_unoriented\tmodal_full\tmodal_genes\trearr_full\trearr_genes"]
        for i in range(3):
            rows.append(f"s{i}\tC\tF\t+A +B\t+A +B\t+A +B\t+A +B\t0\t0")
        (tmp_path / "architecture.tsv").write_text("\n".join(rows) + "\n")
        audit = tmp_path / "audit"
        audit.mkdir()
        orows = ["slug\tclass\tfamily\tstatus\torder_full\torder_genes\torder_full_unoriented\torder_genes_unoriented"]
        for i in range(3):
            orows.append(f"s{i}\tC\tF\tfull\t+A +B\t+A +B\t+A +B\t+A +B")
        (audit / "order.tsv").write_text("\n".join(orows) + "\n")
        (audit / "order_flags.tsv").write_text("tier\tslug\tkind\tfamily\n")
        ok = CliRunner().invoke(app, ["longevity", "arch-control", "--build", str(data), "--out", str(tmp_path), "--audit", str(audit)])
        assert ok.exit_code == 0, ok.output
        (audit / "order_flags.tsv").write_text("tier\tslug\tkind\tfamily\nfull\ts0\torder\tF\n")        # flagged, but at its family's arrangement
        bad = CliRunner().invoke(app, ["longevity", "arch-control", "--build", str(data), "--out", str(tmp_path), "--audit", str(audit)])
        assert bad.exit_code == 1 and "FAILED" in bad.output


# --- the regressions, on the build ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def prepared():
    arch = pd.read_csv(ARCH, sep="\t")
    sp = rna_tests.test_species(paths.BUILD, COMP, paths.DATED_TREE)
    return sp, at.prepare(arch, at.annotation_flags(paths.AUDIT_OUT), sp)


def gls_slope(C, y, X):
    Ci = np.linalg.inv(C)
    beta = np.linalg.solve(X.T @ Ci @ X, X.T @ Ci @ y)
    resid = y - X @ beta
    se = np.sqrt(np.diag(np.linalg.inv(X.T @ Ci @ X)) * (resid @ Ci @ resid / (len(y) - X.shape[1])))
    return beta, se


@needs_build
class TestCells:
    def test_the_species_are_the_default_dataset_with_the_two_sets_by_confidence(self, prepared):
        sp, d = prepared
        assert set(d.slug) == set(sp.slug) and d.slug.is_unique
        assert (d.not_unreliable == d.data_confidence.isin(["high", "caveat", "low"])).all()
        assert (d.high == d.data_confidence.eq("high")).all() and (d.high <= d.not_unreliable).all()
        assert d.not_unreliable.sum() < len(d) and d.annot_flag.isin([0, 1]).all()
        sp2 = pd.read_csv(COMP, sep="\t")
        assert d.slug.isin(sp2.slug[sp2.max_longevity_years > 0]).all()             # every species of the test has a lifespan

    def test_each_slope_is_the_generalised_least_squares_of_the_standardised_measure_with_and_without_the_flag(self, prepared):
        _, d = prepared
        res = at.run_cells(d, paths.DATED_TREE, ["non_coding_frac"], ["Reptilia"])
        cell = d[(d["class"] == "Reptilia") & d.not_unreliable]
        C, order, _ = aa_pgls.class_covariance(paths.DATED_TREE, list(cell.slug))
        cell = cell.set_index("slug").loc[order]
        ok = cell.non_coding_frac.notna().to_numpy()
        i = np.where(ok)[0]
        y = cell.non_coding_frac.to_numpy()[ok]
        z = (y - y.mean()) / y.std(ddof=1)
        ll, lm, fl = cell.ll.to_numpy()[ok], cell.lm.to_numpy()[ok], cell.annot_flag.to_numpy(dtype=float)[ok]
        Cg = C[np.ix_(i, i)]
        for cov, X in ((False, np.column_stack([np.ones(len(z)), ll, lm])), (True, np.column_stack([np.ones(len(z)), ll, lm, fl]))):
            beta, se = gls_slope(Cg, z, X)
            row = res[(res.set == "not_unreliable") & (res.model == "lambda_1") & (res.covariate == cov)].iloc[0]
            assert row.n == len(z) and row.slope == pytest.approx(beta[1], rel=1e-8) and row.se == pytest.approx(se[1], rel=1e-8)
            assert row.df == len(z) - X.shape[1] and row.r == pytest.approx(row.t / np.sqrt(row.t ** 2 + row.df))
        raw = res[(res.set == "not_unreliable") & (res.model == "raw") & ~res.covariate].iloc[0]
        ols = np.linalg.lstsq(np.column_stack([np.ones(len(z)), ll, lm]), z, rcond=None)[0]
        assert raw.slope == pytest.approx(ols[1], rel=1e-6)                          # lambda 0 is the regression that ignores the tree
        ml = res[(res.set == "not_unreliable") & (res.model == "lambda_ml") & ~res.covariate].iloc[0]
        assert 0 <= ml["lambda"] <= 1 and cell.shape[0] >= at.MIN_SPECIES

    def test_the_correction_is_one_holm_and_one_bh_over_every_row_of_the_run_in_the_table_of_the_raw_p(self, prepared):
        _, d = prepared
        res = at.run_cells(d, paths.DATED_TREE, ["non_coding_frac", "gene_overlap_bp", "repeat_genome"], ["Reptilia", "Aves"])
        have = res.p.notna()
        assert {"p", "p_holm_run", "q_run", "family_size"} <= set(res.columns) and res.family_size.eq(have.sum()).all()
        assert np.allclose(res.p_holm_run[have], multipletests(res.p[have], method="holm")[1])
        assert np.allclose(res.q_run[have], multipletests(res.p[have], method="fdr_bh")[1])
        assert res.p_holm_run[~have].isna().all()
        assert set(res.model[have]) == {"lambda_ml", "lambda_1", "raw"} and set(res.covariate[have]) == {False, True}
        assert set(res.set[have]) == {"not_unreliable", "high"} and set(res["class"][have]) == {"Reptilia", "Aves"}
        assert len(res[have]) > 3 * 2 * 2 * 3 + 1                                  # more than within one measure, one group or one model
        by_measure = res[have].groupby("measure").p.transform(lambda s: multipletests(s, method="holm")[1])
        assert (res.p_holm_run[have] >= by_measure - 1e-12).all() and (res.p_holm_run[have] > by_measure + 1e-9).any()

    def test_a_measure_that_is_almost_constant_or_a_group_that_is_small_is_not_assessed_with_its_reason(self, prepared):
        _, d = prepared
        d2 = d.copy()
        m = d2["class"].eq("Reptilia")
        d2.loc[m, "cr_copies"] = 1.0
        d2.loc[d2.index[m][:at.MIN_MINORITY - 1], "cr_copies"] = 2.0
        d2.loc[d2["class"].eq("Aves"), "gene_overlap_bp"] = 5.0
        res = at.run_cells(d2, paths.DATED_TREE, ["cr_copies", "gene_overlap_bp"], ["Reptilia", "Aves"])
        a = res[(res.measure == "cr_copies") & (res["class"] == "Reptilia") & (res.set == "not_unreliable")]
        assert a.p.isna().all() and a.reason.str.contains("off the most common").all()
        b = res[(res.measure == "gene_overlap_bp") & (res["class"] == "Aves")]
        assert b.p.isna().all()
        v = at.judge(res).set_index(["measure", "class"])
        assert v.loc[("cr_copies", "Reptilia"), "verdict"] == "not assessed" and v.loc[("gene_overlap_bp", "Aves"), "verdict"] == "not assessed"
        d4 = d.copy()
        m4 = d4["class"].eq("Reptilia")
        d4.loc[m4, "cr_copies"] = 1.0
        d4.loc[d4.index[m4 & d4.not_unreliable & d4.high][:at.MIN_MINORITY], "cr_copies"] = 2.0
        r4 = at.run_cells(d4, paths.DATED_TREE, ["cr_copies"], ["Reptilia"])
        assert r4[(r4.set == "not_unreliable") & ~r4.covariate].p.notna().all()        # exactly the minimum off the mode is assessed
        d5 = d.copy()
        m5 = d5["class"].eq("Reptilia")
        d5.loc[d5.index[m5][at.MIN_SPECIES - 1:], "cr_length"] = np.nan
        r5 = at.run_cells(d5, paths.DATED_TREE, ["cr_length"], ["Reptilia"])
        assert r5[r5.set == "not_unreliable"].p.isna().all() and r5.reason.str.contains("species with the measure").any()
        d3 = d.copy()
        keep = d3.index[d3["class"].eq("Reptilia") & d3.high][:at.MIN_SPECIES - 1]
        d3.loc[d3["class"].eq("Reptilia") & d3.high & ~d3.index.isin(keep), "high"] = False
        res3 = at.run_cells(d3, paths.DATED_TREE, ["non_coding_frac"], ["Reptilia"])
        h = res3[res3.set == "high"]
        assert h.p.isna().all() and h.reason.str.contains("only 24 species in the set").all()

    def test_the_flag_cannot_enter_a_cell_without_flagged_species_and_its_verdict_is_not_assessed(self, prepared):
        _, d = prepared
        d2 = d.copy()
        d2.loc[d2["class"].eq("Reptilia"), "annot_flag"] = 0
        res = at.run_cells(d2, paths.DATED_TREE, ["non_coding_frac"], ["Reptilia"])
        cov = res[(res.set == "not_unreliable") & res.covariate]
        assert cov.p.isna().all() and cov.reason.str.contains("cannot enter").all()
        assert res[(res.set == "not_unreliable") & ~res.covariate].p.notna().all()
        # the verdict cannot be "counts" there: make the plain cells pass and see
        forced = res.copy()
        forced.loc[forced.p.notna(), ["p", "p_holm_run"]] = [1e-9, 1e-9]
        forced.loc[forced.p.notna(), "slope"] = 0.3
        v = at.judge(forced).iloc[0]
        assert v.verdict == "not assessed" and "flag" in v.reason

    def test_lengths_are_taken_in_log10_and_a_length_of_zero_is_unknown(self, prepared):
        _, d = prepared
        arch = pd.read_csv(ARCH, sep="\t").set_index("slug")
        for c in at.LOGGED:
            want = np.log10(arch.loc[d.slug, c].where(arch.loc[d.slug, c] > 0)).to_numpy()
            assert np.allclose(d[c].to_numpy(dtype=float), want, equal_nan=True) and np.isfinite(d[c]).any()
        for c in set(at.MEASURES) - set(at.LOGGED):
            assert np.allclose(d[c].to_numpy(dtype=float), arch.loc[d.slug, c].to_numpy(dtype=float), equal_nan=True)
        a = arch.reset_index().copy()
        a.loc[a.slug == d.slug.iloc[0], "non_coding_bp"] = 0
        sp = rna_tests.test_species(paths.BUILD, COMP, paths.DATED_TREE)
        assert np.isnan(at.prepare(a, at.annotation_flags(paths.AUDIT_OUT), sp).set_index("slug").at[d.slug.iloc[0], "non_coding_bp"])

    def test_a_test_species_without_an_architecture_or_audit_row_is_refused(self, prepared):
        sp, _ = prepared
        arch = pd.read_csv(ARCH, sep="\t")
        flags = at.annotation_flags(paths.AUDIT_OUT)
        with pytest.raises(ValueError, match="no architecture row"):
            at.prepare(arch[arch.slug != sp.slug.iloc[0]], flags, sp)
        with pytest.raises(ValueError, match="no audit row"):
            at.prepare(arch, flags[flags.slug != sp.slug.iloc[0]], sp)

    def test_a_lifespan_free_species_is_never_in_a_cell(self, prepared):
        sp, d = prepared
        arch = pd.read_csv(ARCH, sep="\t")
        comp = pd.read_csv(COMP, sep="\t")
        no_life = set(arch.slug) - set(comp.slug[comp.max_longevity_years > 0])
        assert no_life and not (no_life & set(d.slug))


@needs_build
@needs_trf
def test_the_table_has_a_row_per_live_species_with_the_distances_the_repeats_and_the_stored_bases(tmp_path):
    import polars as pl
    recs = arec.load(paths.AUDIT)
    sp = pl.read_parquet(paths.BUILD / "species.parquet").to_pandas()
    retired = sp.slug[sp.record_status.eq("retired")]
    in_release = recs[recs.in_release].sort_values("slug")
    pick = pd.concat([in_release.iloc[::70], in_release[in_release.slug.isin(retired)].iloc[:1],
                      in_release[in_release.slug.eq("crocidura-russula")]]).drop_duplicates("slug")
    assert pick.slug.isin(retired).any() and pick.slug.eq("crocidura-russula").any()
    tab, rep = am.table(pick, sp, 2, tmp_path)
    live = [s for s in pick.slug if s not in set(retired)]
    assert list(tab.slug) == sorted(live)                                     # retired species out, sorted by slug
    cro = tab.set_index("slug").loc["crocidura-russula"]
    assert pd.isna(cro.rearr_full) and pd.notna(cro.rearr_genes)
    for tier in am.TIERS:
        col = f"order_{tier}_unoriented"
        for cls, d in tab[tab["class"].notna() & tab[col].notna()].groupby("class"):
            vc = d[col].value_counts()
            modal = sorted(vc.index[vc.values == vc.max()])[0]
            for _, r in d.iterrows():
                assert r[f"rearr_{tier}"] == naive_distance(r[col], modal) or am.breakpoint_distance(r[col], modal) == r[f"rearr_{tier}"]
    s = sp.set_index("slug")
    assert np.allclose(tab.non_coding_frac, (s.loc[tab.slug, "non_coding_bp"] / s.loc[tab.slug, "mtdna_size_bp"]).to_numpy())
    sp2 = sp.copy()
    sp2.loc[sp2.slug == tab.slug.iloc[0], "mtdna_size_bp"] *= 2                  # the stored size, not the record's length, is the share's base
    tab2, _ = am.table(pick, sp2, 2, tmp_path)
    assert tab2.non_coding_frac.iloc[0] == pytest.approx(tab.non_coding_frac.iloc[0] / 2)
    assert set(rep.slug) <= set(tab.slug) and (rep.start >= 1).all() and (rep.start <= rep.end).all()
    for r in tab.itertuples():
        mine = rep[rep.slug == r.slug]
        cover = set()
        for a, b in zip(mine.start, mine.end):
            cover |= set(range(a, b + 1))
        assert r.repeat_genome == pytest.approx(len(cover) / r.genome_length) and max(cover, default=0) <= r.genome_length
    assert (tab.cr_annotated.eq(0) == tab.cr_length.isna()).all() and (tab.cr_annotated.eq(0) == tab.repeat_cr.isna()).all()
    assert tab.cr_copies.dropna().ge(1).all()


# --- the bar, on tables of the shape run_cells writes -------------------------------------------------------
def cells(measure="m", cls="Mammalia", over=None, drop=()):
    """All rows of one measure and group, every condition met; `over`: {(set, covariate, model): {column: value}}; `drop`: keys left out."""
    over = over or {}
    rows = []
    for name in at.SETS:
        for cov in (False, True):
            for model, _ in at.MODELS:
                if (name, cov, model) in drop:
                    continue
                row = {"set": name, "class": cls, "measure": measure, "covariate": cov, "model": model, "n": 300, "n_flagged": 12, "reason": "",
                       "lambda": 0.5, "slope": 0.30, "se": 0.05, "t": 6.0, "df": 297, "p": 1e-6, "r": 0.33,
                       "p_holm_run": 1e-3, "q_run": 1e-3, "family_size": 500}
                if cov:
                    row["slope"] = 0.29
                row.update(over.get((name, cov, model), {}))
                rows.append(row)
    return pd.DataFrame(rows)


def verdict(**kw):
    return at.judge(cells(**kw)).iloc[0]


P, PH = "not_unreliable", "high"
ML, L1 = "lambda_ml", "lambda_1"
FAIL = {"p": 0.2, "p_holm_run": 1.0}
HOLM_FAIL = {"p_holm_run": 0.2}


class TestBar:
    def test_a_measure_that_meets_all_four_conditions_counts(self):
        v = verdict()
        assert v.verdict == "counts" and v.c1_survives_controls and v.c2_holds_on_high and v.c3_survives_run_correction and v.c4_unmoved_by_flag

    def test_both_lambda_forms_must_pass_the_raw_threshold_and_a_pass_in_one_is_a_null(self):
        assert verdict(over={(P, False, ML): FAIL}).verdict == "null"
        assert verdict(over={(P, False, L1): FAIL}).verdict == "null"
        assert verdict(over={(P, False, ML): {"p": 0.05}}).verdict == "null"                       # the threshold is strict

    def test_the_two_lambda_forms_must_agree_in_sign(self):
        assert verdict(over={(P, False, L1): {"slope": -0.3}}).verdict == "null"
        assert verdict(over={(P, False, ML): {"slope": -0.3}}).verdict == "null"

    def test_the_raw_regression_is_only_reported_and_never_decides(self):
        v = verdict(over={(P, False, "raw"): FAIL})
        assert v.verdict == "counts" and v.raw_p == 0.2
        assert verdict(over={(P, False, "raw"): {"slope": -0.3}}).verdict == "counts"

    def test_each_lambda_form_must_hold_on_the_high_confidence_subset_with_the_same_sign(self):
        assert verdict(over={(PH, False, ML): FAIL}).verdict == "not counting"
        assert verdict(over={(PH, False, L1): FAIL}).verdict == "not counting"
        assert verdict(over={(PH, False, ML): {"slope": -0.3}}).verdict == "not counting"
        assert verdict(over={(PH, False, L1): {"slope": -0.3}}).verdict == "not counting"
        assert verdict(over={(PH, False, L1): {"p": 0.05}}).verdict == "not counting"

    def test_each_lambda_form_must_survive_the_correction_over_the_run_on_both_sets(self):
        assert verdict(over={(P, False, ML): HOLM_FAIL}).verdict == "not counting"
        assert verdict(over={(P, False, L1): HOLM_FAIL}).verdict == "not counting"
        assert verdict(over={(PH, False, ML): HOLM_FAIL}).verdict == "not counting"
        assert verdict(over={(PH, False, L1): HOLM_FAIL}).verdict == "not counting"
        assert verdict(over={(P, False, ML): {"p_holm_run": 0.05}}).verdict == "not counting"
        v = verdict(over={(P, False, L1): HOLM_FAIL})
        assert v.c1_survives_controls and v.c2_holds_on_high and not v.c3_survives_run_correction

    def test_a_result_the_flag_moves_is_explained_by_annotation_in_each_way_it_can_move(self):
        for model in (ML, L1):
            assert verdict(over={(P, True, model): {"slope": -0.29}}).verdict == "explained by annotation"
            assert verdict(over={(P, True, model): {"slope": 0.15}}).verdict == "explained by annotation"
            assert verdict(over={(P, True, model): {"slope": 0.45}}).verdict == "explained by annotation"
            assert verdict(over={(P, True, model): HOLM_FAIL}).verdict == "explained by annotation"
        assert verdict(over={(P, True, ML): {"slope": 0.30 * 0.75 - 1e-9}}).verdict == "explained by annotation"
        assert verdict(over={(P, True, ML): {"slope": 0.30 * 0.75 + 1e-9}}).verdict == "counts"        # 25 % is within the tolerance
        assert verdict(over={(P, True, L1): {"slope": 0.30 * 1.25 + 1e-9}}).verdict == "explained by annotation"
        v = verdict(over={(P, True, L1): {"slope": 0.0001}})
        assert v.c1_survives_controls and v.c2_holds_on_high and v.c3_survives_run_correction and not v.c4_unmoved_by_flag

    def test_a_change_of_exactly_a_quarter_of_the_slope_is_not_a_move_and_a_change_of_sign_is(self):
        both = {(P, False, ML): {"slope": 4.0}, (P, False, L1): {"slope": 4.0}}
        assert verdict(over={**both, (P, True, ML): {"slope": 3.0}, (P, True, L1): {"slope": 5.0}}).verdict == "counts"
        assert verdict(over={**both, (P, True, ML): {"slope": -3.0}}).verdict == "explained by annotation"
        assert verdict(over={**both, (P, True, L1): {"slope": -5.0}}).verdict == "explained by annotation"

    def test_a_flag_that_moves_a_result_that_failed_an_earlier_condition_does_not_make_it_explained_by_annotation(self):
        assert verdict(over={(P, True, ML): {"slope": -0.29}, (PH, False, ML): FAIL}).verdict == "not counting"
        assert verdict(over={(P, True, ML): {"slope": -0.29}, (P, False, ML): FAIL}).verdict == "null"
        assert verdict(over={(P, True, ML): {"slope": -0.29}, (P, False, L1): HOLM_FAIL}).verdict == "not counting"

    def test_a_cell_the_flag_cannot_enter_or_with_no_high_subset_is_not_assessed_and_says_why(self):
        v = verdict(drop=[(P, True, ML), (P, True, L1), (P, True, "raw")])
        assert v.verdict == "not assessed" and "flag" in v.reason and pd.isna(v.c4_unmoved_by_flag)
        v = verdict(drop=[(P, True, L1)])
        assert v.verdict == "not assessed" and pd.isna(v.c4_unmoved_by_flag)
        v = verdict(drop=[(PH, False, ML), (PH, False, L1), (PH, False, "raw")])
        assert v.verdict == "not assessed" and "high" in v.reason
        gone = verdict(drop=[(PH, False, ML), (PH, False, L1)], over={(P, False, ML): FAIL})
        assert gone.verdict == "null"                                              # a failed first condition is a null whatever else is missing
        none = at.judge(cells(over={(P, False, ML): {"p": np.nan}})).iloc[0]
        assert none.verdict == "not assessed"

    def test_judge_refuses_a_table_whose_correction_was_not_computed_with_the_p_values(self):
        c = cells()
        for col in ("p_holm_run", "q_run"):
            with pytest.raises(ValueError, match="whole-run correction"):
                at.judge(c.drop(columns=[col]))

    def test_a_measure_counts_where_one_of_its_groups_counts_and_the_summary_names_it(self):
        both = pd.concat([cells(measure="rearr_full", cls="Aves"), cells(measure="rearr_full", cls="Mammalia", over={(P, False, ML): FAIL}),
                          cells(measure="cr_length", cls="Aves", over={(P, False, ML): FAIL})], ignore_index=True)
        s = at.summarise(at.judge(both)).set_index("measure")
        assert bool(s.loc["rearr_full", "measure_counts"]) and s.loc["rearr_full", "counts_in"] == "Aves"
        assert s.loc["rearr_full", "n_counts"] == 1 and s.loc["rearr_full", "n_null"] == 1
        assert not s.loc["cr_length", "measure_counts"] and s.loc["cr_length", "n_null"] == 1
        assert s.loc["gene_overlap_bp", "groups"] == 0 and not s.loc["gene_overlap_bp", "measure_counts"]

    def test_the_verdicts_are_per_measure_and_group_and_do_not_mix_cells_of_other_measures(self):
        both = pd.concat([cells(measure="a"), cells(measure="b", over={(P, False, ML): FAIL})], ignore_index=True)
        v = at.judge(both).set_index("measure")
        assert v.loc["a", "verdict"] == "counts" and v.loc["b", "verdict"] == "null"
        mixed = pd.concat([cells(measure="a", cls="Aves"), cells(measure="a", cls="Mammalia", over={(P, False, L1): FAIL})], ignore_index=True)
        v = at.judge(mixed).set_index("class")
        assert v.loc["Aves", "verdict"] == "counts" and v.loc["Mammalia", "verdict"] == "null"

    def test_a_correction_added_to_the_table_afterwards_is_the_same_as_the_one_computed_with_it(self):
        rng = np.random.default_rng(9)
        raw = cells()
        raw["p"] = rng.uniform(1e-6, 0.3, len(raw))
        out = at.correct_run(raw.drop(columns=["p_holm_run", "q_run", "family_size"]))
        assert np.allclose(np.sort(out.p_holm_run), np.sort(multipletests(raw.p, method="holm")[1]))
        assert (out.p_holm_run >= out.p - 1e-12).all() and out.family_size.eq(len(raw)).all()
