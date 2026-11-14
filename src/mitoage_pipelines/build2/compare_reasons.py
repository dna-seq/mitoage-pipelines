"""Reason classes of the disagreements between a reproduction and the published build 2.

species_in_one_build_only   the row belongs to a species present in only one of the two builds (its
                            inputs differ: another taxonomy dump, release or lifespan table)
taxonomy_label / refseq_label  the value names the input's version, and the input is another version
unexplained                 none of the above: read the row in disagreements.tsv
"""
import re


def species_of_key(table: str, key: str) -> int | None:
    """The species id a disagreement's key names (None for rule rows and review items)."""
    parts = key.split("|")
    head = parts[1] if table == "changelog" else parts[0]
    m = re.match(r"^(\d+)", head or "")
    return int(m.group(1)) if m and table != "review_list" else None


def make(one_side: set[int]):
    def classify(table, field, published, reproduced, row) -> str:
        return "unexplained"

    def post(d: dict) -> str:
        sid = species_of_key(d["table"], d["key"])
        if sid is not None and sid in one_side: return "species_in_one_build_only"
        p, r = d["published"] or "", d["reproduced"] or ""
        if p.startswith("new_taxdump") and r.startswith("new_taxdump"): return "taxonomy_label"
        if "release" in p and "release" in r and d["field"] in ("upstream_version", "reason", "retired_reason", "source_version"):
            return "refseq_label"
        if d["table"] == "changelog" and d["field"] == "(row)" and sid is None and "taxdump" in d["key"]: return "taxonomy_label"
        return "unexplained"
    return classify, post


def classify(table, field, published, reproduced, row) -> str:
    return "unexplained"
