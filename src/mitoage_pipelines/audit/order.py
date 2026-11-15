"""Gene order. Each record's arrangement of the 13 protein genes, the 2 rRNAs and its tRNAs as a signed circular
order (+ on the record's + strand, - on the other), named as build 2 names them, in order of their start on the
record, every copy kept. The control region is left out (its annotation varies more than the genes').

Canonical by rotation: the order is started at COX1 when the record has exactly one COX1, else at the rotation that
is first in string order. A second form is canonical under rotation and reversal (the record read from the other
strand: the order reversed and every sign flipped): it tells a reversed deposit from a rearrangement.

Two tiers: "full" (genes, rRNAs and tRNAs; unknown for a record with no tRNA annotated) and "genes" (the 13 and the
2 rRNAs only). Arrangements are counted per class, order and family at both. A record is flagged when its family has
at least 3 records with a known arrangement at that tier, its own arrangement occurs once there, and every other
record of the family shares one arrangement. What differs is said, never whether it is real: a rearrangement and an
annotation error look the same in a record.
"""
from collections import Counter
from pathlib import Path

import pandas as pd

from ..build2.genome import PROTEIN_GENES
from . import records as rec
from .common import UNKNOWN, Warnings, write_tsv

GENES = PROTEIN_GENES + ("rRNA-12S", "rRNA-16S")
MIN_FAMILY = 3
CANNOT_TELL = "a real rearrangement or an annotation error: cannot tell from the record"


def tokens(g, with_trna: bool = True) -> list[str]:
    """The signed order of the record's genes (and tRNAs), by start position on the record."""
    keep = {"CDS", "rRNA", "tRNA"} if with_trna else {"CDS", "rRNA"}
    return [("+" if f.strand == 1 else "-") + f.name for f in g.features if f.type in keep and (f.type == "tRNA" or f.name in GENES)]


def flip(tok: list[str]) -> list[str]:
    """The same order read from the other strand."""
    return [("-" if t[0] == "+" else "+") + t[1:] for t in reversed(tok)]


def canonical(tok: list[str]) -> str:
    if not tok:
        return ""
    cox1 = [i for i, t in enumerate(tok) if t[1:] == "COX1"]
    if len(cox1) == 1:
        i = cox1[0]
        return " ".join(tok[i:] + tok[:i])
    return min(" ".join(tok[i:] + tok[:i]) for i in range(len(tok)))


def unoriented(tok: list[str]) -> str:
    return min(canonical(tok), canonical(flip(tok)))


def arrangement(g) -> dict:
    full = tokens(g, True)
    genes = tokens(g, False)
    n_trna = sum(1 for t in full if t[1:].startswith("tRNA"))
    return {"genes_present": len({t[1:] for t in genes}), "trnas_annotated": n_trna,
            "order_genes": canonical(genes) if genes else None, "order_genes_unoriented": unoriented(genes) if genes else None,
            "order_full": canonical(full) if n_trna else None, "order_full_unoriented": unoriented(full) if n_trna else None}


def _names(order: str) -> list[str]:
    return [t[1:] for t in order.split(" ")]


def diagnose(own: dict, other: dict, tier: str) -> tuple[str, str]:
    """What separates a record's arrangement from the one the rest of its family shares (rows of `table`)."""
    a, b = own[f"order_{tier}"], other[f"order_{tier}"]
    if own[f"order_{tier}_unoriented"] == other[f"order_{tier}_unoriented"]:
        return "reversed deposit", "the family's arrangement read from the other strand: the record is the other strand"
    ca, cb = Counter(_names(a)), Counter(_names(b))
    if ca != cb:
        lost = sorted((cb - ca).elements()); extra = sorted((ca - cb).elements())
        return "gene content", "; ".join(x for x in (f"absent here: {', '.join(lost)}" if lost else "",
                                                       f"extra here: {', '.join(extra)}" if extra else "") if x)
    if canonical(["+" + x for x in _names(a)]) == canonical(["+" + x for x in _names(b)]):
        sa = {t[1:]: t[0] for t in a.split(" ")}; sb = {t[1:]: t[0] for t in b.split(" ")}
        diff = sorted(x for x in sa if sa[x] != sb.get(x))
        return "strand only", f"same order, other strand for: {', '.join(diff)}"
    if tier == "full" and own["order_genes"] == other["order_genes"]:
        return "tRNAs only", "the 13 genes and 2 rRNAs in the family's order; tRNAs placed or named differently"
    return "order", "the genes in another order"


def row_of(base: dict, g) -> dict:
    """The arrangement row of one record (g None: no record)."""
    if g is None:
        return {**base, "status": UNKNOWN, "reason": "no record", "genes_present": None, "trnas_annotated": None,
                "order_genes": None, "order_genes_unoriented": None, "order_full": None, "order_full_unoriented": None}
    a = arrangement(g)
    if a["order_genes"] is None:
        status, reason = UNKNOWN, "no protein gene or rRNA annotated"
    elif a["order_full"] is None:
        status, reason = "genes only", "no tRNA annotated: the full arrangement is unknown"
    else:
        status, reason = "full", ""
    return {**base, "status": status, "reason": reason, **a}


def table(df: pd.DataFrame) -> pd.DataFrame:
    warn = Warnings("audit order")
    rows = []
    for r in df.itertuples(index=False):
        base = {"slug": r.slug, "refseq": r.refseq, "class": r.class_, "order": r.order, "family": r.family}
        row = row_of(base, rec.genome_of(r) if r.in_release else None)
        if not r.in_release:
            row["reason"] = r.reason
        if row["reason"]:
            warn.add(row["reason"], r.slug)
        rows.append(row)
    warn.emit()
    return pd.DataFrame(rows)


def counts(tab: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tier in ("full", "genes"):
        col = f"order_{tier}"
        known = tab[tab[col].notna()]
        for level in ("class", "order", "family"):
            for name, d in known.groupby(known[level].fillna("(none)"), sort=True):
                c = d[col].value_counts()
                rows.append({"tier": tier, "level": level, "group": name, "records": len(d), "arrangements": len(c),
                             "arrangements_unoriented": d[f"{col}_unoriented"].nunique(),
                             "modal_share": round(int(c.iloc[0]) / len(d), 4)})
        c = known[col].value_counts()
        rows.append({"tier": tier, "level": "all", "group": "all", "records": len(known), "arrangements": len(c),
                     "arrangements_unoriented": known[f"{col}_unoriented"].nunique(), "modal_share": round(int(c.iloc[0]) / len(known), 4) if len(c) else None})
    return pd.DataFrame(rows)


def flags(tab: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tier in ("full", "genes"):
        col = f"order_{tier}"
        known = tab[tab[col].notna() & tab.family.notna()]
        for fam, d in known.groupby("family", sort=True):
            if len(d) < MIN_FAMILY:
                continue
            c = d[col].value_counts()
            if len(c) != 2 or int(c.min()) != 1:
                continue
            odd = d[d[col] == c.index[c.values == 1][0]].iloc[0]
            rest = d[d[col] != odd[col]].iloc[0]
            kind, detail = diagnose(odd, rest, tier)
            rows.append({"tier": tier, "slug": odd.slug, "refseq": odd.refseq, "class": odd["class"], "order": odd["order"], "family": fam,
                         "family_records": len(d), "kind": kind, "detail": detail,
                         "verdict": "a deposit of the other strand, not a rearrangement" if kind == "reversed deposit" else CANNOT_TELL,
                         "arrangement": odd[col], "family_arrangement": rest[col]})
    cols = ["tier", "slug", "refseq", "class", "order", "family", "family_records", "kind", "detail", "verdict", "arrangement", "family_arrangement"]
    return pd.DataFrame(rows, columns=cols)


def arrangements(tab: pd.DataFrame) -> pd.DataFrame:
    """Every distinct arrangement at each tier, with its records per class."""
    rows = []
    for tier in ("full", "genes"):
        col = f"order_{tier}"
        known = tab[tab[col].notna()]
        for a, d in known.groupby(col, sort=True):
            cls = d["class"].fillna("(none)").value_counts()
            rows.append({"tier": tier, "arrangement": a, "records": len(d),
                         "classes": "; ".join(f"{k} {v}" for k, v in sorted(cls.items())),
                         "families": d.family.nunique()})
    return pd.DataFrame(rows)


def main(interim: Path, out: Path) -> pd.DataFrame:
    df = rec.load(interim).rename(columns={"class": "class_"})
    tab = table(df)
    tab = tab.merge(df[["slug", "build_annotation"]], on="slug", how="left", validate="one_to_one")
    write_tsv(tab, out / "order.tsv", ["slug"])
    write_tsv(counts(tab), out / "order_counts.tsv", ["tier", "level", "group"])
    write_tsv(arrangements(tab), out / "order_arrangements.tsv", ["tier", "arrangement"])
    fl = flags(tab)
    write_tsv(fl, out / "order_flags.tsv", ["tier", "family", "slug"])
    print(tab.status.value_counts().sort_index().to_string())
    print(f"flagged: {len(fl)} ({', '.join(f'{t} {n}' for t, n in fl.tier.value_counts().sort_index().items())})")
    return tab
