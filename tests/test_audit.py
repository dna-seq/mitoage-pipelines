"""H13, the audit of the build's records, on the real records of RefSeq release 237 and the downloaded build. Expected
values are computed here from the records themselves: a record read from the other strand, started elsewhere on
its circle, given a stop codon, stripped of an annotation, copied under another name or given another order."""
import copy
import dataclasses
from collections import Counter

import numpy as np
import pandas as pd
import pytest

from mitoage_pipelines import paths
from mitoage_pipelines.audit import code, common, order, records, similar, strand
from mitoage_pipelines.build2 import genome as gm

RELEASE = paths.refseq_release("237")
CONTIG = paths.B2_GENBANK / "contig"
TAXDUMP = paths.taxdump("current-2026-11-13")
pytestmark = pytest.mark.skipif(not (paths.BUILD / "species.parquet").is_file() or not RELEASE.is_file(),
                                reason="needs `mitoage download build` and `mitoage download refseq --from-snapshot`")
needs_taxdump = pytest.mark.skipif(not (TAXDUMP / "nodes.dmp").is_file(), reason="needs `mitoage download taxonomy --from-snapshot`")
HUMAN, FLY = "homo-sapiens", "drosophila-melanogaster"


@pytest.fixture(scope="module")
def recs() -> pd.DataFrame:
    return records.collect(paths.BUILD, RELEASE, CONTIG).set_index("slug", drop=False)


def genome(recs, slug):
    return records.genome_of(recs.loc[slug])


def with_features(g, keep=None, change=None):
    """A copy of the record whose annotated features are filtered (keep) or edited (change: feature -> feature)."""
    g = copy.deepcopy(g)
    feats = [f for f in g.features if keep is None or keep(f)]
    if change:
        feats = [change(f) for f in feats]
    return dataclasses.replace(g, features=feats)


def gc_sign(seq: str) -> str:
    b = gm.count_bases(seq)
    return common.sign(gm.skew(b["g"], b["c"]))


# --- the records ------------------------------------------------------------------------------
def test_every_species_has_a_row_and_its_own_record(recs):
    sp = pd.read_parquet(paths.BUILD / "species.parquet")
    assert set(recs.slug) == set(sp.slug) and len(recs) == len(sp)
    live = recs[recs.in_release]
    for r in live.itertuples(index=False):
        assert records.genome_of(r).record.accession_version == r.refseq
    gone = recs[~recs.in_release]
    assert (gone.reason != "").all() and gone.genbank.isna().all()
    assert set(sp[sp.record_status != "retired"].slug) <= set(live.slug)       # every live species has its record
    for c in records.TAXONOMY_COLUMNS:                 # an empty label in the build is a missing one, not a group named ""
        blank = set(sp.slug[sp[c].fillna("") == ""])
        assert set(recs.slug[recs[c].isna()]) == blank and not (recs[c] == "").any()


def test_the_builds_annotation_is_compared_not_taken(recs):
    """Each record's features as RefSeq deposits them against the build's: same for an unpatched record, and a
    feature only one side has is named, on its side."""
    assert set(recs.build_annotation) <= {"same", "differs", "no record"}
    assert (recs.build_annotation == "no record").sum() == (~recs.in_release).sum()
    assert recs.loc[HUMAN, "build_annotation"] == "same"
    assert (recs[recs.build_annotation == "differs"].build_annotation_detail != "").all()
    g = genome(recs, HUMAN)
    built = records.build_features(paths.BUILD)[HUMAN]
    assert records.compare_with_build(g, built) == ("same", "")
    nd6 = next(f for f in g.features if f.name == "ND6")
    st, de = records.compare_with_build(g, built - {(nd6.name, nd6.type, nd6.start, nd6.end, nd6.strand)})
    assert st == "differs" and de == f"only in the record: ND6 CDS {nd6.start}..{nd6.end} (-)"
    st, de = records.compare_with_build(with_features(g, keep=lambda f: f.name != "ND6"), built)
    assert st == "differs" and de.startswith("only in the build: ND6 CDS")
    assert records.compare_with_build(None, built) == ("no record", "")


# --- strand -----------------------------------------------------------------------------------
def test_reverse_complement_flips_the_class_and_the_skew(recs):
    vert = recs[recs.in_release & recs.vertebrate.eq(True)]
    mirror = {"conventional": "reversed", "reversed": "conventional", "mixed": "mixed", "unknown": "unknown"}
    seen = Counter()
    for slug in list(vert.slug[:: max(1, len(vert) // 40)]) + [HUMAN]:
        g = genome(recs, slug)
        a = strand.classify(g, True)
        rc = records.reverse_complement(g)
        b = strand.classify(rc, True)
        assert b["orientation"] == mirror[a["orientation"]], slug
        assert strand.classify(records.reverse_complement(rc), True)["orientation"] == a["orientation"]
        assert (a["copies_conventional"], a["copies_opposite"]) == (b["copies_opposite"], b["copies_conventional"])
        s, t = gc_sign(g.record.sequence), gc_sign(rc.record.sequence)
        assert {s, t} == {"+", "-"} or s == t == "0"
        seen[a["orientation"]] += 1
    assert strand.classify(genome(recs, HUMAN), True)["orientation"] == "conventional"
    assert gc_sign(genome(recs, HUMAN).record.sequence) == "-"          # the human L-strand deposit: G-poor


def test_skew_beside_the_class_is_the_builds(recs):
    sk = strand.skews(paths.BUILD)
    g = genome(recs, HUMAN)
    want = gm.make_stats(0, "total_mtDNA", gm.count_bases(g.record.sequence), "genbank")
    assert sk.loc[HUMAN, "gc_skew"] == pytest.approx(want["gc_skew"]) and sk.loc[HUMAN, "at_skew"] == pytest.approx(want["at_skew"])


def test_one_gene_on_the_other_strand_is_mixed(recs):
    g = genome(recs, HUMAN)
    flip = lambda f: dataclasses.replace(f, strand=-f.strand) if f.name == "COX2" else f
    c = strand.classify(with_features(g, change=flip), True)
    assert c["orientation"] == "mixed" and c["opposite_genes"] == "COX2" and c["copies_opposite"] == 1
    assert "majority is conventional" in c["reason"]


def test_missing_annotation_is_unknown_not_conventional(recs):
    g = genome(recs, HUMAN)
    assert strand.classify(with_features(g, keep=lambda f: f.type != "CDS"), True)["orientation"] == "unknown"
    assert strand.classify(None, None)["orientation"] == "unknown"
    present = [x for x in gm.PROTEIN_GENES]
    six = set(present[:6]); seven = set(present[:7])
    assert strand.classify(with_features(g, keep=lambda f: f.type != "CDS" or f.name in six), True)["orientation"] == "unknown"
    c = strand.classify(with_features(g, keep=lambda f: f.type != "CDS" or f.name in seven), True)
    assert c["orientation"] == "conventional" and "not annotated" in c["reason"]
    no_nd6 = strand.classify(with_features(g, keep=lambda f: f.name != "ND6"), True)
    assert no_nd6["orientation"] == "conventional" and no_nd6["missing_genes"] == "ND6"


def test_outside_vertebrata_is_unknown(recs):
    assert recs.loc[FLY, "vertebrate"] == False and recs.loc[HUMAN, "vertebrate"] == True   # noqa: E712
    g = genome(recs, FLY)
    for x in (g, records.reverse_complement(g)):
        c = strand.classify(x, False)
        assert c["orientation"] == "unknown" and "outside Vertebrata" in c["reason"] and c["genes_annotated"] == 13


# --- gene order --------------------------------------------------------------------------------
def free_position(g, after: int = 1000) -> int:
    """A base no named feature and no CDS, rRNA or tRNA covers: a place to start the circle."""
    n = len(g.record.sequence)
    covered = np.zeros(n + 2, dtype=bool)
    for segs in [f.segments for f in g.features] + [f.segments for f in g.record.features if f.key in ("CDS", "rRNA", "tRNA")]:
        for s in segs:
            covered[s.start:s.end + 1] = True
    k = next(i for i in range(after, n) if not covered[i] and not covered[i + 1])
    return k


def test_rotation_does_not_change_the_arrangement(recs):
    for slug in (HUMAN, FLY, "gallus-gallus"):
        g = genome(recs, slug)
        k = free_position(g)
        r = records.rotate(g, k)
        assert r.record.sequence == g.record.sequence[k:] + g.record.sequence[:k]
        a, b = order.arrangement(g), order.arrangement(r)
        assert a == b, slug
        assert order.tokens(r) != order.tokens(g) or k == 0      # the raw order did move
    with pytest.raises(ValueError):
        f = next(f for f in genome(recs, HUMAN).features if f.name == "COX1")
        records.rotate(genome(recs, HUMAN), f.start + 10)


def test_reverse_complement_keeps_the_unoriented_arrangement(recs):
    for slug in (HUMAN, FLY):
        g = genome(recs, slug)
        a, b = order.arrangement(g), order.arrangement(records.reverse_complement(g))
        assert a["order_full_unoriented"] == b["order_full_unoriented"] and a["order_genes_unoriented"] == b["order_genes_unoriented"]
        assert a["order_full"] != b["order_full"]
        assert b["order_genes"] == order.canonical(order.flip(order.tokens(g, False)))


@pytest.fixture(scope="module")
def bovids(recs):
    """Three bovids that share one arrangement at both tiers (the largest such group, its first three by slug)."""
    d = recs[recs.in_release & recs.family.eq("Bovidae")]
    rows = [(s, order.arrangement(genome(recs, s))) for s in d.slug]
    groups = Counter((a["order_full"], a["order_genes"]) for _, a in rows if a["order_full"])
    top = groups.most_common(1)[0][0]
    pick = sorted(s for s, a in rows if (a["order_full"], a["order_genes"]) == top)[:3]
    assert len(pick) == 3
    return pick


def family_table(recs, slugs, odd=None):
    rows = []
    for i, s in enumerate(slugs):
        r = recs.loc[s]
        g = records.genome_of(r)
        if odd is not None and i == len(slugs) - 1:
            g = odd(g)
        rows.append(order.row_of({"slug": s, "refseq": r.refseq, "class": r["class"], "order": r["order"], "family": r.family}, g))
    return pd.DataFrame(rows)


def test_a_family_that_shares_one_arrangement_flags_nothing(recs, bovids):
    assert order.flags(family_table(recs, bovids)).empty


def test_a_family_of_two_is_too_small_to_share_one(recs, bovids):
    """Two records that differ: neither is the odd one out. With a third that shares one's arrangement, the other is."""
    assert order.flags(family_table(recs, bovids[:2], records.reverse_complement)).empty
    assert len(order.flags(family_table(recs, bovids, records.reverse_complement))) == 2


@pytest.mark.parametrize("edit,tiers,kind", [
    ("reverse", {"full", "genes"}, "reversed deposit"),
    ("strand", {"full", "genes"}, "strand only"),
    ("drop_nd6", {"full", "genes"}, "gene content"),
    ("swap_nd1_nd2", {"full", "genes"}, "order"),
    ("swap_trna", {"full"}, "tRNAs only"),
])
def test_an_edited_record_is_flagged_with_what_differs(recs, bovids, edit, tiers, kind):
    def swap(a, b):
        return lambda f: dataclasses.replace(f, name=b) if f.name == a else dataclasses.replace(f, name=a) if f.name == b else f
    trnas = [f.name for f in genome(recs, bovids[-1]).features if f.type == "tRNA"]
    edits = {"reverse": records.reverse_complement,
             "strand": lambda g: with_features(g, change=lambda f: dataclasses.replace(f, strand=-f.strand) if f.name == "ATP8" else f),
             "drop_nd6": lambda g: with_features(g, keep=lambda f: f.name != "ND6"),
             "swap_nd1_nd2": lambda g: with_features(g, change=swap("ND1", "ND2")),
             "swap_trna": lambda g: with_features(g, change=swap(trnas[0], trnas[5]))}
    fl = order.flags(family_table(recs, bovids, edits[edit]))
    assert set(fl.tier) == tiers and (fl.slug == bovids[-1]).all() and (fl.kind == kind).all()
    assert (fl.family_records == 3).all()
    if kind != "reversed deposit":
        assert (fl.verdict == order.CANNOT_TELL).all()
    if edit == "drop_nd6":
        assert (fl.detail == "absent here: ND6").all()


def test_a_trna_less_record_has_no_full_arrangement(recs, bovids):
    g = genome(recs, HUMAN)
    bare = with_features(g, keep=lambda f: f.type != "tRNA")
    a, b = order.arrangement(g), order.arrangement(bare)
    assert b["order_full"] is None and b["trnas_annotated"] == 0 and a["order_genes"] == b["order_genes"]
    row = order.row_of({"slug": HUMAN}, bare)
    assert row["status"] == "genes only" and "no tRNA" in row["reason"]
    tab = family_table(recs, bovids, lambda x: with_features(x, keep=lambda f: f.type != "tRNA"))
    assert order.flags(tab).empty                                  # two at the full tier: under the family minimum
    cnt = order.counts(tab)
    assert cnt.set_index(["tier", "level"]).loc[("full", "family"), "records"] == 2


# --- genetic code -----------------------------------------------------------------------------
def translate(c: dict, table: dict) -> str:
    nt = c["nt"]
    return "".join(table.get(nt[i:i + 3].replace("T", "U"), "X") for i in range(0, len(nt) - 2, 3))


def test_translation_under_the_declared_code_is_refseqs(recs):
    """Our reading of each CDS under its declared table gives RefSeq's own /translation (the first residue aside, which
    initiation makes M), and the internal stops counted are the * inside it."""
    live = recs[recs.in_release]
    checked = 0
    for slug in list(live.slug[:: max(1, len(live) // 60)]) + [HUMAN, FLY]:
        g = genome(recs, slug)
        cds = code.cds_of(g)
        dec, _ = code.declared(cds)
        table = code.code_table(dec)
        if table is None:
            continue
        raw = [f for f in g.record.features if f.key == "CDS" and gm.feature_name(f) and gm.feature_name(f)[0] in gm.PROTEIN_GENES]
        for c, f in zip(cds, raw):
            ref = f.qualifiers.get("translation", [""])[0]
            ours = translate(c, table)
            if not ref or f.qualifiers.get("transl_except") or "X" in ours:
                continue
            body = ours[:-1] if ours.endswith("*") and len(ours) == len(ref) + 1 else ours
            assert body[1:] == ref[1:], (slug, c["gene"])
            assert code.internal_stops(c, table)[0] == ours[:-1].count("*")
            checked += 1
    assert checked > 500


def test_the_standard_code_counts_every_internal_tga(recs):
    g = genome(recs, HUMAN)
    std = code.code_table(1)
    for c in code.cds_of(g):
        codons = [c["nt"][i:i + 3] for i in range(0, len(c["nt"]) - 2, 3)][:-1]
        want = sum(x in ("TGA", "TAA", "TAG") for x in codons)
        assert code.internal_stops(c, std)[0] == want
        assert code.internal_stops(c, code.code_table(2))[0] == sum(x in ("AGA", "AGG", "TAA", "TAG") for x in codons)


def test_an_inserted_stop_is_counted_once_and_survives_reverse_complement(recs):
    g = genome(recs, HUMAN)
    f = next(f for f in g.features if f.name == "COX1")
    p = f.start - 1 + 300                                             # codon 101 of a + strand gene, codon_start 1
    seq = g.record.sequence[:p] + "TAA" + g.record.sequence[p + 3:]
    mut = dataclasses.replace(g, record=dataclasses.replace(g.record, sequence=seq))
    before = {x["gene"]: x for x in code.assess(g, 2)[1]}
    after = {x["gene"]: x for x in code.assess(mut, 2)[1]}
    assert after["COX1"]["stops_declared"] == before["COX1"]["stops_declared"] + 1 == 1
    assert after["COX1"]["first_stop_declared"] == 101
    assert all(after[x]["stops_declared"] == before[x]["stops_declared"] for x in after if x != "COX1")
    rc = {x["gene"]: x for x in code.assess(records.reverse_complement(mut), 2)[1]}
    assert {x: rc[x]["stops_declared"] for x in rc} == {x: after[x]["stops_declared"] for x in after}
    row = code.assess(mut, 2)[0]
    assert row["clean_declared"] == "no" and row["neither_clean"] == "yes" and row["genes_with_stops_declared"] == "COX1:1"


def test_codon_start_sets_the_frame(recs):
    """COX1 annotated one base earlier with codon_start 2 reads the same codons and the same stops."""
    g = genome(recs, HUMAN)
    shifted = copy.deepcopy(g)
    f = next(f for f in shifted.record.features if f.key == "CDS" and gm.feature_name(f)[0] == "COX1")
    assert f.segments[0].strand == 1 and f.qualifiers.get("codon_start", ["1"]) == ["1"]
    f.segments[0] = gm.Segment(f.segments[0].start - 1, f.segments[0].end, 1)
    f.qualifiers["codon_start"] = ["2"]
    a = next(c for c in code.cds_of(g) if c["gene"] == "COX1")
    b = next(c for c in code.cds_of(shifted) if c["gene"] == "COX1")
    assert b["nt"] == a["nt"] and b["pos"] == a["pos"]
    assert code.internal_stops(b, code.code_table(1)) == code.internal_stops(a, code.code_table(1))


def test_an_unknown_code_is_unknown_not_clean(recs):
    g = genome(recs, HUMAN)
    assert code.code_table(99) is None and code.code_table(None) is None
    row, genes = code.assess(g, 99)
    assert row["clean_expected"] == "unknown" and row["stops_expected"] is None and row["code_match"] == "no"
    assert row["clean_declared"] == "yes" and row["neither_clean"] == "no"
    assert "expected table 99 not known" in row["reason"]
    bare = copy.deepcopy(g)
    for f in bare.record.features:
        f.qualifiers.pop("transl_table", None)
    row, _ = code.assess(bare, None)
    assert (row["declared_table"], row["code_match"], row["clean_declared"], row["neither_clean"]) == (None, "unknown", "unknown", "unknown")
    none = dataclasses.replace(g, record=dataclasses.replace(g.record, features=[f for f in g.record.features if f.key != "CDS"]))
    row, genes = code.assess(none, 2)
    assert genes == [] and row["stops_declared"] is None and row["clean_expected"] == "unknown" and "no protein gene" in row["reason"]


def test_neither_clean_has_three_answers():
    Y, N, U = "yes", "no", "unknown"
    assert code.neither_clean(N, N) == Y
    assert code.neither_clean(Y, N) == N and code.neither_clean(N, Y) == N and code.neither_clean(Y, U) == N
    assert code.neither_clean(N, U) == U and code.neither_clean(U, U) == U


@needs_taxdump
def test_expected_code_comes_from_the_taxonomy_with_merges():
    mc = code.MitoCodes(TAXDUMP)
    assert mc.of(9606)[0] == 2 and mc.of(7227)[0] == 5 and mc.of(7955)[0] == 2
    old, new = next((o, n) for o, n in sorted(mc.merged.items()) if mc.mgc.get(n))
    got, note = mc.of(old)
    assert got == mc.mgc[new] and f"merged into {new}" in note
    assert mc.of(None)[0] is None and mc.of(10 ** 9)[0] is None


# --- duplicates and nearest neighbours ----------------------------------------------------------
@pytest.fixture(scope="module")
def alignments():
    return similar.load_alignments(paths.BUILD, "codon"), similar.load_alignments(paths.BUILD, "aa")


def test_identity_counts_equal_a_count_by_hand(alignments):
    nt, _ = alignments
    slugs = sorted(nt["COX1"])[::400][:6]
    seqs = similar.concatenate(nt, slugs)
    m, k = similar.identity_counts(seqs, similar.NT)
    for i in range(len(slugs)):
        for j in range(len(slugs)):
            both = [(a, b) for a, b in zip(seqs[i], seqs[j]) if a in "ACGT" and b in "ACGT"]
            assert k[i, j] == len(both) and m[i, j] == sum(a == b for a, b in both)
    assert (m == m.T).all() and (k == k.T).all()


def test_two_identical_records_pair_and_others_do_not(recs, alignments):
    nt, _ = alignments
    a, b, c = HUMAN, "gallus-gallus", FLY
    copied = {g: {**aln, "copy-of-human": aln[a]} if a in aln else aln for g, aln in nt.items()}
    slugs = sorted([a, b, c, "copy-of-human"])
    tax = recs.set_index("slug")[["refseq", "scientific_name", "class", "order", "family", "genus"]].loc[[a, b, c]]
    tax.loc["copy-of-human"] = ["NC_X", "Copy copy", "Mammalia", "Rodentia", "Muridae", "Mus"]
    m, k = similar.identity_counts(similar.concatenate(copied, slugs), similar.NT)
    pr = similar.pairs(tax.loc[slugs], m, k, {})
    assert len(pr) == 1
    p = pr.iloc[0]
    assert {p.slug_a, p.slug_b} == {a, "copy-of-human"} and p.identity == 1.0 and p.differences == 0
    assert p.near_identical == "yes" and p.relation == "same class, different order" and p.same_sequence == "no"
    same = similar.pairs(tax.loc[slugs], m, k, {b: "h", c: "h"})
    assert set(map(tuple, same[["slug_a", "slug_b"]].values)) == {tuple(sorted((a, "copy-of-human"))), tuple(sorted((b, c)))}
    assert same.set_index(["slug_a", "slug_b"]).same_sequence.to_dict() == {tuple(sorted((b, c))): "yes", tuple(sorted((a, "copy-of-human"))): "no"}


def test_a_record_given_another_order_sits_outside_it(recs, alignments, bovids):
    _, aa = alignments
    b1, b2 = bovids[:2]
    slugs = sorted([b1, b2, HUMAN])
    tax = recs.set_index("slug")[["refseq", "scientific_name", "class", "order", "family", "genus"]].loc[slugs].copy()
    m, k = similar.identity_counts(similar.concatenate(aa, slugs), similar.AA)
    i = {s: n for n, s in enumerate(slugs)}
    assert m[i[b1], i[b2]] / k[i[b1], i[b2]] > m[i[b1], i[HUMAN]] / k[i[b1], i[HUMAN]]
    nn = similar.nearest(tax, m, k).set_index("slug")
    assert nn.loc[b1, "out_of_order"] == "no" and nn.loc[b1, "nn_slug"] == b2
    assert nn.loc[HUMAN, "out_of_order"] == "unknown"                  # no other primate here: unknown, not in place
    tax.loc[b1, "order"] = tax.loc[HUMAN, "order"]
    nn = similar.nearest(tax, m, k).set_index("slug")
    assert nn.loc[b1, "out_of_order"] == "yes" and nn.loc[b1, "margin_order"] > 0 and nn.loc[b1, "best_in_order_slug"] == HUMAN
    assert nn.loc[b1, "out_of_class"] == "no"


# --- tables and warnings -----------------------------------------------------------------------
def test_tables_are_sorted_and_rewritten_byte_for_byte(tmp_path):
    df = pd.DataFrame({"slug": ["b", "a", "c"], "x": [0.1234567891, None, 3.0]})
    p1 = common.write_tsv(df, tmp_path / "1.tsv", ["slug"])
    p2 = common.write_tsv(df.sample(frac=1, random_state=1), tmp_path / "2.tsv", ["slug"])
    assert p1.read_bytes() == p2.read_bytes()
    assert p1.read_text().splitlines()[1:] == ["a\t", "b\t0.123457", "c\t3"]
    with pytest.raises(ValueError):
        common.write_tsv(pd.concat([df, df]), tmp_path / "3.tsv", ["slug"])


def test_warnings_are_counted_by_reason():
    w = common.Warnings("step")
    for s in ("x", "y"):
        w.add("same reason", s)
    w.add("another")
    assert w.lines() == ["step: another: 1", "step: same reason: 2 (first: x)"]
