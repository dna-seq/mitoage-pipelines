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


def main(published: Path, reproduced: Path, out: Path, reason=None) -> dict:
    from .compare_reasons import classify
    reason = reason or classify
    summary, dis, sets = [], [], []
    for t in KEYS:
        s, d, st = compare_table(t, read(published, t), read(reproduced, t), reason)
        summary += s; dis += d; sets.append(st)
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
    (out / "summary.json").write_text(json.dumps({"rows": sets, "values": tot}, indent=1) + "\n")
    print(pl.DataFrame(sets))
    print(tot)
    if dd.height:
        with pl.Config(tbl_rows=60, fmt_str_lengths=60):
            print(by_reason)
    return {"rows": sets, "values": tot, "fields": summary}
