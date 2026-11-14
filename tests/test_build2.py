"""The build-2 rules on small real cases, with expected values computed in the test: the base-label fix of
build 1, join() segments, the lifespan precedence and its tie-break, name matching, id minting, the
change classes, and the agreement of record_status with the changelog."""
import pytest

from mitoage_pipelines.build2 import assemble, build1, diff, genome, lifespan, validate
from mitoage_pipelines.build2.lifespan import Obs

SEQ = "ATGAAACCCGGGTTTAAAGGGCCCTTTAAACCCGGGTAA"     # 39 nt


def record(features: str, seq: str = SEQ) -> str:
    """A minimal GenBank record with these feature lines."""
    lines = [seq[i:i + 60] for i in range(0, len(seq), 60)]
    origin = "".join(f"{i * 60 + 1:>9} {' '.join(l[j:j + 10] for j in range(0, len(l), 10)).lower()}\n" for i, l in enumerate(lines))
    return (f"LOCUS       NC_000001              {len(seq)} bp    DNA     circular VRT 01-JAN-2020\n"
            "DEFINITION  Testus testus mitochondrion, complete genome.\nACCESSION   NC_000001\nVERSION     NC_000001.1\n"
            "SOURCE      mitochondrion Testus testus\n  ORGANISM  Testus testus\n            Eukaryota; Metazoa; Chordata.\n"
            "FEATURES             Location/Qualifiers\n"
            f"     source          1..{len(seq)}\n                     /db_xref=\"taxon:12345\"\n"
            f"{features}ORIGIN\n{origin}//\n")


# --- build 1's base labels --------------------------------------------------------------------
def test_build1_base_labels_are_complemented_and_fixed():
    seq = "AAAAACCGTT"                                    # on the RefSeq strand: A 5, C 2, G 1, T 2
    true = genome.count_bases(seq)
    assert (true["a"], true["c"], true["g"], true["t"]) == (5, 2, 1, 2)
    # build 1 wrote the complement's counts under each label: column "G" = C, "C" = G, "A" = T, "T" = A
    row = build1.Composition(label_g=true["c"], label_c=true["g"], label_a=true["t"], label_t=true["a"], other=0,
                             length=10, gc_percent=30.0, non_coding=None, non_coding_no_dloop=None)
    fixed = build1.fixed_counts(row)
    assert {k: fixed[k] for k in "acgt"} == {k: true[k] for k in "acgt"}
    st = genome.make_stats(1, "total_mtDNA", fixed, "mitoage")
    assert st["gc_skew"] == round((1 - 2) / 3, 5) and st["at_skew"] == round((5 - 2) / 7, 5)
    assert st["gc_percent"] == 30.0                       # GC% does not depend on the labels
    as_labelled = genome.make_stats(1, "total_mtDNA", {"a": row.label_a, "c": row.label_c, "g": row.label_g, "t": row.label_t, "other": 0}, "mitoage")
    assert as_labelled["gc_skew"] == -st["gc_skew"] and as_labelled["at_skew"] == -st["at_skew"]   # the old sign


def test_halves_round_up_as_in_javascript():
    assert genome.js_round(2.5, 0) == 3 and genome.js_round(-2.5, 0) == -2 and genome.js_round(0.125, 2) == 0.13
    assert round(2.5) == 2                               # Python's own rounding would differ


# --- join() features ---------------------------------------------------------------------------
JOIN = ('     CDS             join(1..9,11..21)\n                     /gene="ND3"\n                     /codon_start=1\n'
        '                     /transl_table=2\n'
        '     tRNA            complement(30..36)\n                     /product="tRNA-Leu"\n                     /anticodon=(pos:complement(33..35),aa:Leu,seq:taa)\n')


def test_join_counts_annotated_segments_only():
    g = genome.annotate(record(JOIN), 2)
    nd3 = next(f for f in g.features if f.name == "ND3")
    assert (nd3.start, nd3.end) == (1, 21)                        # the span, as a feature row
    counts, length = genome.compute_region(g, "ND3")
    assert length == 9 + 11 == 20                                 # not 21: base 10 is skipped
    want = genome.count_bases(SEQ[0:9] + SEQ[10:21])
    assert counts == want
    codons = genome.gene_codons(g, "ND3")
    seg = (SEQ[0:9] + SEQ[10:21]).replace("T", "U")
    assert sum(codons.values()) == len(seg) // 3 and codons[seg[0:3]] >= 1
    # what build 1 counted: the full span; the diff classifies it as the full-span rule
    full = genome.count_bases(SEQ[0:21])
    b1 = {"length_bp": 21, "other": 0, **{k: full[k] for k in "acgt"}}
    c = diff.classify_region({"evidence": g, "evidence_text": record(JOIN), "primary": None}, "ND3", b1)
    assert c["diff_class"] == "rule_full_span" and c["scope"] == "rule"


def test_complement_join_reads_segments_in_reverse():
    segs, partial = genome.parse_location("complement(join(1..5,10..12))")
    assert [(s.start, s.end, s.strand) for s in segs] == [(10, 12, -1), (1, 5, -1)] and not partial
    assert genome.extract("AACCGTTTTGCA", segs) == genome.reverse_complement("GCA") + genome.reverse_complement("AACCG")


def test_trna_isoacceptor_from_anticodon():
    g = genome.annotate(record(JOIN), 2)
    assert [f.name for f in g.features if f.type == "tRNA"] == ["tRNA-Leu(UUR)"]    # anticodon TAA reads UUR


# --- lifespan precedence and matching ----------------------------------------------------------
def obs(source, value, record_type="observed", quality="acceptable", rid=None, aquatic=None):
    return Obs(source, "v", rid or f"{source}:{value}:{record_type}", "Testus testus", None, value, record_type, quality,
               None, None, None, None, "lic", aquatic=aquatic)


def test_precedence_anage_first_then_the_taxon_source():
    fb = [obs("FishBase", 30.0), obs("FishBase", 50.0, "estimated", "low")]
    assert lifespan.choose(fb + [obs("AnAge", 12.0)], "fish").source == "AnAge"           # AnAge always wins
    assert lifespan.choose(fb, "fish").value_years == 30.0          # a reported value beats a larger estimate
    assert lifespan.choose([obs("FishBase", 50.0, "estimated", "low")], "fish").value_years == 50.0
    assert lifespan.choose([obs("FishBase", 40.0, quality="questionable"), obs("FishBase", 20.0)], "fish").value_years == 20.0
    assert lifespan.choose(fb, "bird") is None and lifespan.choose(fb, "other_vertebrate") is None
    assert lifespan.choose([obs("USGS_BBL", 9.5, "wild")], "bird").value_years == 9.5
    slb = [obs("SeaLifeBase", 70.0, aquatic=False), obs("SeaLifeBase", 15.0, aquatic=True)]
    assert lifespan.choose(slb, "invertebrate").value_years == 15.0                       # aquatic species only


def test_tie_prefers_a_stated_specimen_origin():
    popchar = obs("FishBase", 25.0, "observed", rid="popchar:1")
    wild = obs("FishBase", 25.0, "wild", rid="species.LongevityWild:7")
    captive = obs("FishBase", 25.0, "captive", rid="species.LongevityCaptive:7")
    assert lifespan.choose([popchar, wild, captive], "fish") is wild     # first with an origin, not the first listed
    assert lifespan.choose([wild, popchar], "fish") is wild
    assert lifespan.choose([popchar, obs("FishBase", 25.0, "observed", rid="popchar:2")], "fish").source_record_id == "popchar:1"


class Tax:
    def __init__(self, sci, syn): self.sci, self.syn = sci, syn


def test_matching_order_and_ambiguity():
    tax = Tax({"aa bb": [1], "cc dd": [2, 3], "ee ff": [4]}, {"old name": [1]})
    o = [Obs("AnAge", "v", "1", "Aa bb", None, 1.0, "wild", None, None, None, None, None, "l"),
         Obs("AnAge", "v", "2", "Old name", None, 1.0, "wild", None, None, None, None, None, "l"),
         Obs("AnAge", "v", "3", "Cc dd", None, 1.0, "wild", None, None, None, None, None, "l"),
         Obs("AnAge", "v", "4", "Ee ff gg", None, 1.0, "wild", None, None, None, None, None, "l"),
         Obs("USGS_BBL", "v", "5", "Xx (yy) zz", None, 1.0, "wild", None, None, None, None, None, "l")]
    m = lifespan.match_names(o, tax, {})
    get = lambda r: m[lifespan.rec_key(r.source, r.source_name, r.speccode)]
    assert (get(o[0]).method, get(o[0]).status) == ("ncbi_scientific_name", "ok_conflict_primary")   # beats the synonym
    assert (get(o[1]).method, get(o[1]).status) == ("ncbi_synonym", "conflict_secondary")
    assert get(o[2]).status == "ambiguous" and get(o[2]).species_taxid is None       # two species: never guessed
    assert (get(o[3]).method, get(o[3]).species_taxid) == ("binomial_of_trinomial", 4)
    assert get(o[4]).status == "unmatched"
    assert [lifespan.usable(get(x)) for x in o] == [True, False, False, True, False]


def test_fishbase_synonym_table_used_after_ncbi_names():
    tax = Tax({"new genus": [9]}, {})
    o = Obs("FishBase", "v", "1", "Old genus", None, 3.0, "observed", None, None, None, None, None, "l", speccode=77)
    syn = {"FishBase": {77: [("New genus", "accepted name"), ("Bad genus", "misapplied name")]}}
    m = lifespan.match_names([o], tax, syn)[lifespan.rec_key("FishBase", "Old genus", 77)]
    assert (m.method, m.species_taxid, m.status) == ("fishbase_synonym", 9, "ok")


# --- ids, names, groups ------------------------------------------------------------------------
def test_ids_and_slugs_of_new_species():
    assert assemble.mint_id(9606) == 10_009_606
    assert assemble.slugify("Pœcilia  reticulata") == "p-cilia-reticulata"
    assert assemble.slugify("Rhamphochromis esox (Lake Malawi)") == "rhamphochromis-esox-lake-malawi"
    assert assemble.slugify("Ćichla") == "cichla"


def test_name_change_classes():
    assert assemble.name_change("Carduelis spinus", "Spinus spinus") == "genus_transfer"
    assert assemble.name_change("Mustela vison", "Neovison vison") == "genus_transfer"
    assert assemble.name_change("Felis silvestris", "Felis silvestris") is None
    assert assemble.name_change("Aotus lemurinus", "Aotus lemurina") == "epithet_spelling_or_gender"
    assert assemble.name_change("Bos grunniens", "Poephagus grunniens") == "genus_transfer"
    assert assemble.name_change("Carassius auratus", "Carassius langsdorfii") == "different_species_name"


def test_taxon_groups_and_class_scheme():
    fish = "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Actinopterygii; Neopterygii"
    lamprey = "Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Cyclostomata; Hyperoartia; Petromyzontiformes"
    assert lifespan.taxon_group(fish) == "fish" and lifespan.taxon_group(lamprey) == "fish"
    assert lifespan.taxon_group("Eukaryota; Metazoa; Mollusca; Bivalvia") == "invertebrate"
    assert lifespan.harmonised_class(lamprey, "Hyperoartia") == "Cephalaspidomorphi"     # build 1's label
    assert lifespan.harmonised_class("Eukaryota; Metazoa; Chordata; Vertebrata; Sauropsida; Testudines", "") == "Reptilia"


# --- record_status and changelog ---------------------------------------------------------------
def test_record_status_agrees_with_changelog():
    ch = lambda key, table, field, ct: {"row_key": key, "table": table, "field": field, "change_type": ct}
    species = [{"id": 1, "record_status": "original"}, {"id": 2, "record_status": "updated"},
               {"id": 3, "record_status": "renamed"}, {"id": 4, "record_status": "retired"},
               {"id": 10_000_005, "record_status": "new"}]
    log = [ch("2:ND1", "region_stats", "a", "correction"), ch("3", "species", "scientific_name", "renamed"),
           ch("4", "species", "record_status", "retired"), ch("10000005", "species", "(row)", "added"), ch("*", "region_stats", "a", "correction")]
    assert validate.problems(species, log) == []
    assert validate.problems(species, log + [ch("1:COX1", "codon_usage", "UUU", "correction")]) == ["1: original but 1 changelog rows"]
    assert validate.problems(species, log[1:]) == ["2: updated without a changelog row"]
    assert validate.problems([{"id": 10_000_006, "record_status": "updated"}], [ch("10000006", "species", "x", "y")]) == \
        ["10000006: a minted id on a build 1 species"]


def test_window_search_is_circular():
    seq = "GGGGAAAATTCC"
    target = genome.count_bases("CCGG")                    # only across the origin: positions 11, 12, 1, 2
    assert diff.window_starts(seq, 4, target) == [11]
