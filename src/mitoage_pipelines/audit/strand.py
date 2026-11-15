"""Which strand each record deposits, against the vertebrate convention: 12 of the 13 protein genes on one strand
and ND6 on the other. A record that follows it with the 12 on its + strand is conventional; with the 12 on its -
strand, reversed (the deposit is the other strand); with some copies each way, mixed; and unknown when the
convention cannot be read from it (no record, outside Vertebrata, fewer than 7 of the 13 genes annotated).

Beside it, the record's GC/AT balance over the whole sequence (gc_skew = (G-C)/(G+C), at_skew = (A-T)/(A+T), on
the + strand as deposited, from the build's region_stats), and the build's own orientation column for comparison.
"""
from pathlib import Path

import pandas as pd
import polars as pl

from ..build2.genome import PROTEIN_GENES
from . import records as rec
from .common import UNKNOWN, Warnings, sign, write_tsv

MIN_GENES = 7          # of the 13: fewer say too little about a pattern of two strands


def conventional_strand(gene: str) -> int:
    return -1 if gene == "ND6" else 1


def classify(g, vertebrate) -> dict:
    """orientation (conventional, reversed, mixed, unknown) and the reason, from the annotated protein genes."""
    cds = [f for f in g.features if f.type == "CDS"] if g is not None else []
    present = [x for x in PROTEIN_GENES if any(f.name == x for f in cds)]
    agree = [f.name for f in cds if f.strand == conventional_strand(f.name)]
    opposite = [f.name for f in cds if f.strand != conventional_strand(f.name)]
    missing = [x for x in PROTEIN_GENES if x not in present]
    pattern = " ".join(f"{x}:{''.join('+' if f.strand == 1 else '-' for f in cds if f.name == x) or '.'}" for x in PROTEIN_GENES)
    out = {"genes_annotated": len(present), "copies_conventional": len(agree), "copies_opposite": len(opposite),
           "opposite_genes": ",".join(x for x in PROTEIN_GENES if x in opposite), "missing_genes": ",".join(missing),
           "pattern": pattern}
    note = f"; not annotated: {', '.join(missing)}" if missing else ""
    if g is None:
        return {**out, "orientation": UNKNOWN, "reason": "no record"}
    if not vertebrate:
        return {**out, "orientation": UNKNOWN, "reason": "outside Vertebrata: the convention is the vertebrates'"}
    if not cds:
        return {**out, "orientation": UNKNOWN, "reason": "no protein gene annotated"}
    if len(present) < MIN_GENES:
        return {**out, "orientation": UNKNOWN, "reason": f"only {len(present)} of the 13 protein genes annotated (at least {MIN_GENES} needed)"}
    if not opposite:
        return {**out, "orientation": "conventional", "reason": f"every copy of the {len(present)} annotated genes on its conventional strand{note}"}
    if not agree:
        return {**out, "orientation": "reversed", "reason": f"every copy of the {len(present)} annotated genes on the other strand{note}"}
    most = "conventional" if len(agree) > len(opposite) else "reversed" if len(opposite) > len(agree) else "neither"
    return {**out, "orientation": "mixed",
            "reason": f"{len(opposite)} of {len(cds)} copies against the convention ({out['opposite_genes']}); the majority is {most}{note}"}


def skews(build: Path) -> pd.DataFrame:
    rs = pl.read_parquet(build / "region_stats.parquet").filter(pl.col("region") == "total_mtDNA")
    return rs.select("slug", "gc_skew", "at_skew").to_pandas().set_index("slug")


def summary(tab: pd.DataFrame) -> pd.DataFrame:
    """Per class (and over all records, and over all vertebrate records): how many records of each orientation have
    a negative or positive GC and AT skew."""
    rows = []
    groups = [("all", tab), ("all vertebrates", tab[tab.vertebrate.eq(True)])]
    groups += [(c, d) for c, d in tab.groupby(tab["class"].fillna("(none)"), sort=True)]
    for name, d in groups:
        for o, e in d.groupby("orientation", sort=True):
            rows.append({"group": name, "orientation": o, "records": len(e),
                         "gc_skew_negative": int((e.gc_skew_sign == "-").sum()), "gc_skew_positive": int((e.gc_skew_sign == "+").sum()),
                         "gc_skew_zero_or_unknown": int((~e.gc_skew_sign.isin(["-", "+"])).sum()),
                         "at_skew_positive": int((e.at_skew_sign == "+").sum()), "at_skew_negative": int((e.at_skew_sign == "-").sum()),
                         "at_skew_zero_or_unknown": int((~e.at_skew_sign.isin(["-", "+"])).sum())})
    return pd.DataFrame(rows)


def table(df: pd.DataFrame, build: Path) -> pd.DataFrame:
    sk = skews(build)
    built = pl.read_parquet(build / "mtdna_structure.parquet").select("slug", "orientation").to_pandas().set_index("slug").orientation
    warn = Warnings("audit strand")
    rows = []
    for r in df.itertuples(index=False):
        g = rec.genome_of(r) if r.in_release else None
        c = classify(g, r.vertebrate if r.in_release else None)
        if not r.in_release:
            c["reason"] = r.reason
        warn.add(f"{c['orientation']}: {c['reason'].split(';')[0].split(' (')[0]}" if c["orientation"] in (UNKNOWN, "mixed") else c["orientation"], r.slug)
        gs = sk.gc_skew.get(r.slug); at = sk.at_skew.get(r.slug)
        rows.append({"slug": r.slug, "refseq": r.refseq, "class": r.class_, "order": r.order, "family": r.family,
                     "vertebrate": r.vertebrate, **c, "gc_skew": gs, "at_skew": at, "gc_skew_sign": sign(gs),
                     "at_skew_sign": sign(at), "build_orientation": built.get(r.slug, UNKNOWN)})
    warn.emit()
    cols = ["slug", "refseq", "class", "order", "family", "vertebrate", "orientation", "reason", "genes_annotated",
            "copies_conventional", "copies_opposite", "opposite_genes", "missing_genes", "gc_skew", "at_skew", "gc_skew_sign",
            "at_skew_sign", "build_orientation", "pattern"]
    return pd.DataFrame(rows)[cols]


def main(interim: Path, build: Path, out: Path) -> pd.DataFrame:
    df = rec.load(interim).rename(columns={"class": "class_"})
    tab = table(df, build)
    tab = tab.merge(df[["slug", "build_annotation"]], on="slug", how="left", validate="one_to_one")
    write_tsv(tab, out / "strand.tsv", ["slug"])
    write_tsv(summary(tab), out / "strand_summary.tsv", ["group", "orientation"])
    print(tab.orientation.value_counts().sort_index().to_string())
    return tab
