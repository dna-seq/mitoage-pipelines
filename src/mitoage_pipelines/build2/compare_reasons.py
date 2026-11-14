"""Reason classes of the disagreements between a reproduction and the published build 2.

species_in_one_build_only        the row belongs to a species present in only one of the two builds (the
                                 inputs differ: another taxonomy dump, release or lifespan table)
input_label                      the value names an input's version (RefSeq release, taxonomy dump), and the
                                 reproduction read another version: the text is equal once the names are swapped
same_version_record_differs      the species keeps the same RefSeq accession.version in both builds but its
                                 genome rows differ, so the record text read differs (NCBI re-annotated it in
                                 place); its changelog rows, provenance and species fields follow
unexplained                      none of the above: read the row in disagreements.tsv
"""
import re

from . import pins

GENOME_TABLES = ("region_stats", "codon_usage", "genome_features")


def species_of_key(table: str, key: str) -> int | None:
    """The species id a disagreement's key names (None for rule rows and review items)."""
    if table == "review_list": return None
    parts = key.split("|")
    head = parts[1] if table == "changelog" and len(parts) > 1 else parts[0]
    m = re.match(r"^(\d+)", head or "")
    return int(m.group(1)) if m else None


def classify(table, field, published, reproduced, row) -> str:
    return "unexplained"


def post_classifier(dis: list[dict], one_side: set[int], same_refseq: set[int], labels: dict):
    """The second pass, which sees every disagreement: returns post(d) -> reason."""
    record_differs = {s for d in dis if d["table"] in GENOME_TABLES and (s := species_of_key(d["table"], d["key"])) is not None} & same_refseq
    swaps = [(pins.REFSEQ_VERSION, labels.get("refseq")), (pins.TAXONOMY_VERSION, labels.get("taxonomy"))]
    swaps = [(a, b) for a, b in swaps if b and a != b]

    def post(d: dict) -> str:
        sid = species_of_key(d["table"], d["key"])
        if sid is not None and sid in one_side: return "species_in_one_build_only"
        p, r = d["published"], d["reproduced"]
        if p is not None and r is not None and swaps:
            q = p
            for a, b in swaps: q = q.replace(a, b)
            if q == r: return "input_label"
        if sid is not None and sid in record_differs: return "same_version_record_differs"
        return "unexplained"
    return post
