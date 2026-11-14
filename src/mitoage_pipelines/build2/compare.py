"""The reproduction against the published build 2: per table, the row sets (by key) and, per field, how many
values are identical, equal within a tolerance (floats: 1e-9 absolute), or differ. Every disagreement is
listed with a reason class. This is the only step that reads the published build 2 tables.
"""
import json
from pathlib import Path

import polars as pl

KEYS = {
    "species": ["id"],
    "region_stats": ["species_id", "region"],
    "codon_usage": ["species_id", "region"],
    "genome_features": ["species_id", "feature_key"],
    "lifespan_observations": ["species_id", "source", "source_record_id"],
    "changelog": ["table", "row_key", "field", "old_value", "_n"],
    "review_list": ["review_type", "source", "source_name", "scientific_name", "_n"],
}
TOL = 1e-9


def read(folder: Path, table: str) -> pl.DataFrame:
    if table == "review_list":
        p = folder / "review_list.csv"
        if not p.is_file(): raise FileNotFoundError(f"no {p}")
        df = pl.read_csv(p, infer_schema=False)
    else:
        p = folder / f"{table}.parquet"
        if not p.is_file(): raise FileNotFoundError(f"no {p}")
        df = pl.read_parquet(p)
    if "_n" in KEYS[table]:      # several rows with one key: number them in file order
        k = [c for c in KEYS[table] if c != "_n"]
        df = df.with_columns(pl.int_range(pl.len()).over([pl.col(c).fill_null("\u0000") for c in k]).alias("_n"))
    return df


def as_text(v) -> str | None:
    return None if v is None else str(v)


def compare_table(table: str, pub: pl.DataFrame, rep: pl.DataFrame, reason) -> tuple[list[dict], list[dict], dict]:
    keys = KEYS[table]
    fields = [c for c in pub.columns if c not in keys and c in rep.columns]
    j = pub.join(rep, on=keys, how="full", suffix="__rep", coalesce=True, nulls_equal=True)
    both = j.filter(pl.col(f"{fields[0]}__rep").is_not_null() | pl.col(fields[0]).is_not_null()) if fields else j
    pub_keys = set(pub.select(keys).iter_rows()); rep_keys = set(rep.select(keys).iter_rows())
    sets = {"table": table, "published": pub.height, "reproduced": rep.height, "in_both": len(pub_keys & rep_keys),
            "only_published": len(pub_keys - rep_keys), "only_reproduced": len(rep_keys - pub_keys)}
    inner = pub.join(rep, on=keys, how="inner", suffix="__rep", nulls_equal=True)
    summary, dis = [], []
    for f in fields:
        a, b = inner[f], inner[f"{f}__rep"]
        same = (a == b).fill_null(False) | (a.is_null() & b.is_null())
        tol = pl.Series([False] * inner.height)
        if a.dtype.is_float() and b.dtype.is_float():
            tol = (~same) & ((a - b).abs() <= TOL).fill_null(False)
        diff = ~(same | tol)
        summary.append({"table": table, "field": f, "compared": inner.height, "identical": int(same.sum()),
                        "within_tolerance": int(tol.sum()), "differs": int(diff.sum())})
        for r in inner.filter(diff).select(keys + [f, f"{f}__rep"]).iter_rows(named=True):
            key = "|".join(as_text(r[k]) or "" for k in keys if k != "_n")
            dis.append({"table": table, "key": key, "field": f, "published": as_text(r[f]), "reproduced": as_text(r[f"{f}__rep"]),
                        "reason": reason(table, f, r[f], r[f"{f}__rep"], r)})
    for side, ks in (("only_published", pub_keys - rep_keys), ("only_reproduced", rep_keys - pub_keys)):
        for k in sorted(ks, key=lambda x: tuple("" if v is None else str(v) for v in x)):
            dis.append({"table": table, "key": "|".join("" if v is None else str(v) for v, n in zip(k, keys) if n != "_n"),
                        "field": "(row)", "published": "present" if side == "only_published" else None,
                        "reproduced": "present" if side == "only_reproduced" else None, "reason": reason(table, "(row)", side, None, None)})
    return summary, dis, sets


def compare_alignments(published: Path, reproduced: Path) -> dict:
    """Per alignment file: the SHA-256 of the uncompressed FASTA in each manifest; and the gene lists."""
    pa, ra = (json.loads((d / "alignments.json").read_text()) for d in (published, reproduced))
    pf, rf = {f["file"]: f for f in pa["files"]}, {f["file"]: f for f in ra["files"]}
    same = sorted(k for k in pf if k in rf and pf[k]["sha256"] == rf[k]["sha256"])
    return {"files_published": len(pf), "files_reproduced": len(rf), "identical_sha256": len(same),
            "differ": sorted(set(pf) ^ set(rf) | {k for k in pf if k in rf and pf[k]["sha256"] != rf[k]["sha256"]}),
            **{k: pa[k] == ra[k] for k in ("absent", "duplicates", "codonNotes", "aligner", "source")}}


def tree_splits(node: dict, out: set, ages: dict) -> frozenset:
    if "leaf" in node and not node["children"]:
        return frozenset([node["leaf"]])
    below = frozenset().union(*(tree_splits(c, out, ages) for c in node["children"]))
    out.add(below)
    return below


def compare_tree(published: Path, reproduced: Path) -> dict:
    """phylogeny.json: identical bytes, or else the leaves, the clades (splits) and the branch lengths that differ."""
    import hashlib, sys
    sys.setrecursionlimit(1_000_000)
    pb, rb = published.read_bytes(), reproduced.read_bytes()
    res = {"sha256_published": hashlib.sha256(pb).hexdigest(), "sha256_reproduced": hashlib.sha256(rb).hexdigest()}
    res["identical_bytes"] = res["sha256_published"] == res["sha256_reproduced"]
    p, r = json.loads(pb), json.loads(rb)
    sp, sr = set(), set()
    lp, lr = tree_splits(p["root"], sp, {}), tree_splits(r["root"], sr, {})
    res.update(version_published=p["version"], version_reproduced=r["version"], leaves_published=len(lp), leaves_reproduced=len(lr),
               leaves_in_both=len(lp & lr), clades_published=len(sp), clades_reproduced=len(sr), clades_in_both=len(sp & sr))
    return res


def main(published: Path, reproduced: Path, out: Path, reason=None, site_tree: Path | None = None) -> dict:
    from . import compare_reasons
    psp, rsp = read(published, "species"), read(reproduced, "species")
    one_side = set(psp["id"].to_list()) ^ set(rsp["id"].to_list())
    same_refseq = set(psp.join(rsp, on="id", suffix="_r").filter(pl.col("refseq") == pl.col("refseq_r"))["id"].to_list())
    info = reproduced / "build-info.json"
    labels = json.loads(info.read_text()) if info.is_file() else {}
    summary, dis, sets = [], [], []
    for t in KEYS:
        s, d, st = compare_table(t, read(published, t), read(reproduced, t), reason or compare_reasons.classify)
        summary += s; dis += d; sets.append(st)
    post = compare_reasons.post_classifier(dis, one_side, same_refseq, labels)
    for d in dis:
        if d["reason"] == "unexplained": d["reason"] = post(d)
    out.mkdir(parents=True, exist_ok=True)
    sm = pl.DataFrame(summary); sm.write_csv(out / "fields.tsv", separator="\t")
    pl.DataFrame(sets).write_csv(out / "rows.tsv", separator="\t")
    dd = pl.DataFrame(dis, schema={"table": pl.String, "key": pl.String, "field": pl.String, "published": pl.String,
                                   "reproduced": pl.String, "reason": pl.String})
    dd.write_csv(out / "disagreements.tsv", separator="\t")
    by_reason = dd.group_by("table", "field", "reason").len().sort("table", "field", "reason")
    by_reason.write_csv(out / "disagreements_by_reason.tsv", separator="\t")
    tot = {"values_compared": int(sm["compared"].sum()), "identical": int(sm["identical"].sum()),
           "within_tolerance": int(sm["within_tolerance"].sum()), "differs": int(sm["differs"].sum())}
    extra = {}
    if (published / "alignments" / "alignments.json").is_file() and (reproduced / "alignments" / "alignments.json").is_file():
        extra["alignments"] = compare_alignments(published / "alignments", reproduced / "alignments")
    if site_tree and site_tree.is_file() and (reproduced / "phylogeny.json").is_file():
        extra["phylogeny"] = compare_tree(site_tree, reproduced / "phylogeny.json")
    (out / "summary.json").write_text(json.dumps({"rows": sets, "values": tot, **extra}, indent=1) + "\n")
    for k, v in extra.items(): print(k, v)
    print(pl.DataFrame(sets))
    print(tot)
    if dd.height:
        with pl.Config(tbl_rows=60, fmt_str_lengths=60):
            print(by_reason)
    return {"rows": sets, "values": tot, "fields": summary, **extra}
