"""What explains each difference between a build 1 value (base labels fixed) and the row build 2 emits
for the same species (site repository: pipeline/steps/diff.ts). The evidence is upstream data only: the
RefSeq record the row was computed from, the INSDC record its COMMENT names (when downloaded), and the
record's gene features. Classes, tried in this order for a single region:

  rule_full_span (correction, rule)            build 1 equals the full start..end span of a join() feature
  upstream_reannotation_same_version           build 1 equals the INSDC record: RefSeq re-annotated since 2015
  upstream_reannotation_same_version_probable  build 1 equals the span of the gene feature
  coordinates_differ_origin_unresolved         build 1 equals a window of the sequence <= 60 bp away
  build1_disagrees_with_cited_record           none of the above
Aggregates (total_pc, total_rRNA) and total_pc codons inherit their members' classes (derived_from_...);
a changed accession makes every difference upstream_version_change or refseq_successor.
"""
import re

import numpy as np

from .build1 import fixed_counts
from .genome import CODONS, PROTEIN_GENES, compute_region, count_bases, feature_keys, gene_codons, make_stats
from .lifespan import js_str

RF = ("length_bp", "a", "c", "g", "t", "other", "gc_percent")
SINGLE = PROTEIN_GENES + ("rRNA-12S", "rRNA-16S", "D_Loop")
AGG_MEMBERS = {"total_pc": list(PROTEIN_GENES), "total_rRNA": ["rRNA-12S", "rRNA-16S"], "total_tRNA": []}
GENE_CLS_ORDER = ["gap_fill_from_refseq", "upstream_reannotation_same_version", "upstream_reannotation_same_version_probable",
                  "rule_full_span", "coordinates_differ_origin_unresolved", "build1_disagrees_with_cited_record"]
CORRECTION_CLASSES = {"gap_fill_from_refseq", "rule_full_span", "coordinates_differ_origin_unresolved", "build1_disagrees_with_cited_record"}
REASONS = {
    "duplicate": "Build-1 CSV has two rows for this region; this row does not match the cited RefSeq annotation and is dropped (the sibling row matches)",
    "full_span": "Build 1 counted the full start..end span of a segmented (join) feature, including the untranslated frameshift nucleotide(s); build 2 counts only annotated segments",
    "probable": "Build-1 value equals the span of the RefSeq gene feature; the CDS/rRNA coordinates were changed (trimmed/extended) after build 1 while the gene feature kept the old span",
    "disagrees": "Build-1 value is not reproducible from the cited RefSeq sequence and annotation (no contiguous window with these counts near the feature)",
    "gap_fill": "Region missing from the build-1 CSV export but annotated in the cited RefSeq record; computed from the sequence (D-loop may come from misc_feature control region)",
}


def cls(change_type, diff_class, reason, scope) -> dict:
    return {"change_type": change_type, "diff_class": diff_class, "reason": reason, "scope": scope}


def vals(s: dict) -> dict:
    return {k: s[k] for k in RF}


def vals_csv(sid: int, region: str, row) -> dict:
    return vals(make_stats(sid, region, fixed_counts(row), "mitoage"))


def same_acgt(x: dict, y: dict) -> bool:
    return x["a"] == y["a"] and x["c"] == y["c"] and x["g"] == y["g"] and x["t"] == y["t"]


def gene_spans(text: str) -> list[tuple[int, int]]:
    parts = text.split("\nFEATURES")
    feats = parts[1].split("\nORIGIN")[0] if len(parts) > 1 else ""
    out = []
    for m in re.finditer(r"^ {5}gene\s+(\S+)", feats, flags=re.M):
        nums = [int(x) for x in re.findall(r"\d+", m.group(1))]
        if len(nums) >= 2: out.append((min(nums), max(nums)))
    return out


def window_starts(seq: str, L: int, target: dict) -> list[int]:
    """1-based starts of the circular windows of length L with exactly these A, C, G, T counts."""
    n = len(seq)
    if L <= 0 or L > n: return []
    s = np.frombuffer((seq + seq[:L]).encode(), dtype=np.uint8)
    ok = np.ones(n, dtype=bool)
    for b, k in ((b"A", "a"), (b"C", "c"), (b"G", "g"), (b"T", "t")):
        pre = np.concatenate([[0], np.cumsum(s == b[0])])
        ok &= (pre[L:L + n] - pre[:n]) == target[k]
    return (np.nonzero(ok)[0] + 1).tolist()


def classify_region(inp: dict, region: str, b1: dict) -> dict:
    g = inp["evidence"]; seq = g.record.sequence
    f = None
    if region in SINGLE:
        f = next((x for x in g.features if (x.type == "D-loop" if region == "D_Loop" else x.name == region)), None)
    if f and f.start <= f.end and same_acgt(count_bases(seq[f.start - 1:f.end]), b1):
        return cls("correction", "rule_full_span", REASONS["full_span"], "rule")
    if inp["primary"]:
        p = compute_region(inp["primary"][1], region)
        if p and same_acgt(p[0], b1):
            return cls("value_update", "upstream_reannotation_same_version",
                       f"Build-1 value equals the source INSDC record {inp['primary'][0]}; the RefSeq record (same accession.version) was re-annotated after build 1", "row")
    if f:
        if f.start <= f.end:
            for gs, ge in gene_spans(inp["evidence_text"]):
                ov = min(ge, f.end) - max(gs, f.start)
                if ov > 0 and not (gs == f.start and ge == f.end) and same_acgt(count_bases(seq[gs - 1:ge]), b1):
                    return cls("value_update", "upstream_reannotation_same_version_probable", REASONS["probable"], "row")
        L = b1["length_bp"] - b1["other"]
        if b1["other"] == 0 and L > 0:
            hits = window_starts(seq, L, b1)
            if hits:
                off = hits[0] - f.start
                for h in hits:
                    if abs(h - f.start) < abs(off): off = h - f.start
                if abs(off) <= 60:
                    cur = compute_region(g, region)
                    dlen = b1["length_bp"] - (cur[1] if cur else 0)
                    return cls("correction", "coordinates_differ_origin_unresolved",
                               f"Build-1 value reproduces a window of the cited sequence offset {off} bp from the current feature start (length diff {dlen} bp); 2015 annotation not retrievable, so upstream re-annotation vs build-1 extraction error is unresolved", "row")
    return cls("correction", "build1_disagrees_with_cited_record", REASONS["disagrees"], "row")


def row_text(v: dict) -> str:
    return f"length_bp={v['length_bp']};a={v['a']};c={v['c']};g={v['g']};t={v['t']}"


def diff_genome(inp: dict) -> dict:
    """inp: id, csv_composition {region: [rows]}, csv_codons {region: counts}, emitted_stats, emitted_codons,
    evidence (Genome), evidence_text, primary (accession, Genome) or None, acc (from, to, kind) or None,
    baseline_features, emitted_features. Returns changes, aggregates and removed row keys."""
    sid = inp["id"]
    changes, aggregates, removed = [], [], []
    emitted = {s["region"]: s for s in inp["emitted_stats"]}
    acc = inp["acc"]
    def acc_class(what):
        reason = (f"RefSeq {acc[0]} was suppressed; build 2 uses the live RefSeq mitogenome of the same species {acc[1]} ({what})"
                  if acc[2] == "refseq_successor" else f"RefSeq sequence version changed {acc[0]} -> {acc[1]} ({what})")
        return cls("value_update", acc[2], reason, "row")
    def change(table, key, field, old, new, c):
        changes.append({"table": table, "row_key": key, "field": field, "old_value": old, "new_value": new, **c})

    mismatch, dup, gap = {}, set(), set()
    for region, rows in inp["csv_composition"].items():
        cur = emitted.get(region); key = f"{sid}:{region}"
        if cur is None:
            removed.append(f"region_stats:{key}")
            b1 = vals_csv(sid, region, rows[0])
            c = acc_class("region not annotated in the new record; build 1 row kept as retired") if acc else cls(
                "value_update", "region_absent_in_current_record",
                "Region of the build-1 CSV is not annotated in the current RefSeq record; build 1 row kept as retired", "row")
            change("region_stats", key, "(row)", row_text(b1), None, c); continue
        cv = vals(cur)
        if len(rows) > 1:
            dup.add(region)
            for r in rows:
                b1 = vals_csv(sid, region, r)
                if same_acgt(b1, cv): continue
                change("region_stats", key, "(row)", row_text(b1), "dropped (duplicate)",
                       cls("correction", "duplicate_csv_row", REASONS["duplicate"], "rule"))
            continue
        b1 = vals_csv(sid, region, rows[0])
        if same_acgt(b1, cv) and b1["other"] == cv["other"]: continue
        mismatch[region] = {"b1": b1, "cur": cv, "c": acc_class("values recomputed from the new record") if acc else classify_region(inp, region, b1)}
    for s in inp["emitted_stats"]:
        if s["region"] in inp["csv_composition"]: continue
        gap.add(s["region"])
        c = acc_class("region absent from build 1, annotated in the new record") if acc else cls("correction", "gap_fill_from_refseq", REASONS["gap_fill"], "rule")
        change("region_stats", f"{sid}:{s['region']}", "(row)", None, f"length_bp={s['length_bp']};gc_percent={js_str(s['gc_percent'])}", c)
    if not acc:
        for region, m in mismatch.items():
            members = AGG_MEMBERS.get(region)
            if members is None or m["c"]["diff_class"] == "upstream_reannotation_same_version": continue
            sub = []
            for r in members:
                mm = mismatch.get(r)
                if mm and r not in dup: sub.append((r, mm["c"]["diff_class"], mm["b1"]["length_bp"] - mm["cur"]["length_bp"]))
                elif r in gap: sub.append((r, "gap_fill_from_refseq", -emitted[r]["length_bp"]))
            if not sub: continue
            classes = sorted(dict.fromkeys(x[1] for x in sub), key=lambda c: GENE_CLS_ORDER.index(c) if c in GENE_CLS_ORDER else -1)
            dsum = sum(x[2] for x in sub)
            m["c"] = cls("correction" if any(c in CORRECTION_CLASSES for c in classes) else "value_update",
                         f"derived_from_{'+'.join(classes)}",
                         f"Aggregate of member regions with differences ({', '.join(x[0] for x in sub)}; their length differences sum to {dsum} bp, aggregate differs by {m['b1']['length_bp'] - m['cur']['length_bp']} bp)",
                         "rule" if all(c in ("rule_full_span", "gap_fill_from_refseq") for c in classes) else "row")
    for region, m in mismatch.items():
        if region in AGG_MEMBERS: aggregates.append((region, m["c"]))
        for f in RF:
            if m["b1"][f] == m["cur"][f]: continue
            change("region_stats", f"{sid}:{region}", f, js_str(m["b1"][f]), js_str(m["cur"][f]), m["c"])

    # codon usage
    emitted_codons = {u["region"]: u["counts"] for u in inp["emitted_codons"]}
    csv_codons = inp["csv_codons"]
    codon_mismatch = {}
    for region, counts in csv_codons.items():
        cur = emitted_codons.get(region)
        if cur is None: continue
        codon_mismatch[region] = any(counts.get(c, 0) != cur.get(c, 0) for c in CODONS)

    def codon_class(region):
        if acc: return acc_class("codons recounted from the new record")
        b1 = csv_codons[region]
        if inp["primary"] and region != "total_pc":
            p = gene_codons(inp["primary"][1], region)
            if p and all(p.get(c, 0) == b1.get(c, 0) for c in CODONS):
                return cls("value_update", "upstream_reannotation_same_version",
                           "Build-1 codon counts equal the source INSDC record; the RefSeq CDS was re-annotated after build 1 (same accession.version)", "row")
        if region == "total_pc":
            genes = [g for g in PROTEIN_GENES if g in csv_codons]
            is_sum = bool(genes) and all(b1.get(c, 0) == sum(csv_codons[g].get(c, 0) for g in genes) for c in CODONS)
            diff_genes = [g for g in genes if codon_mismatch.get(g)]
            gap_genes = [g for g in PROTEIN_GENES if g not in csv_codons and g in emitted_codons]
            if is_sum and (diff_genes or gap_genes):
                cl = sorted(set([codon_class(g)["diff_class"] for g in diff_genes] + (["gap_fill_from_refseq"] if gap_genes else [])))
                gaps_txt = f" (gene rows missing from build 1: {', '.join(gap_genes)})" if gap_genes else ""
                return cls("value_update" if all(c.startswith("upstream") for c in cl) else "correction", f"derived_from_{'+'.join(cl)}",
                           f"Build-1 total_pc equals the sum of its gene rows; differs because of gene-level differences in {', '.join(diff_genes + gap_genes)}{gaps_txt}",
                           "rule" if cl == ["gap_fill_from_refseq"] else "row")
            if genes and not is_sum:
                return cls("correction", "build1_total_pc_not_sum_of_genes",
                           "Build-1 total_pc codon counts are not the sum of the build-1 gene rows and do not match the cited record", "row")
        g = mismatch.get(region)
        if g and region not in dup and g["c"]["diff_class"] != "rule_full_span":
            return {**g["c"], "reason": f"Same classification as the base-composition difference of this gene: {g['c']['reason']}"}
        return cls("correction", "build1_disagrees_with_cited_record",
                   "Build-1 codon counts are not reproducible from the cited RefSeq CDS (base composition of the gene matches or is not diagnostic)", "row")

    for region, differs in codon_mismatch.items():
        if not differs: continue
        c = codon_class(region); b1 = csv_codons[region]; cur = emitted_codons[region]
        for cd in CODONS:
            if b1.get(cd, 0) == cur.get(cd, 0): continue
            change("codon_usage", f"{sid}:{region}", cd, str(b1.get(cd, 0)), str(cur.get(cd, 0)), c)
    for region, cur in emitted_codons.items():
        if region in csv_codons: continue
        c = acc_class("codon row absent from build 1, present in the new record") if acc else cls(
            "correction", "gap_fill_from_refseq", "Codon row missing from build-1 CSV; counted from the cited RefSeq CDS annotation", "rule")
        change("codon_usage", f"{sid}:{region}", "(row)", None, f"codons={sum(cur.get(c, 0) for c in CODONS)}", c)
    for region in csv_codons:
        if region in emitted_codons: continue
        removed.append(f"codon_usage:{sid}:{region}")
        c = acc_class("gene not annotated in the new record; build 1 row kept as retired") if acc else cls(
            "value_update", "region_absent_in_current_record",
            "Gene of the build-1 CSV is not annotated in the current RefSeq record; build 1 row kept as retired", "row")
        change("codon_usage", f"{sid}:{region}", "(row)", "present", None, c)

    # genome features, keyed name#n in genome order
    def keyed(fs):
        return dict(zip(feature_keys(fs), fs))
    before, after = keyed(inp["baseline_features"]), keyed(inp["emitted_features"])
    fc = acc_class("feature table of the new record") if acc else cls(
        "value_update", "upstream_reannotation_same_version",
        "Feature coordinates differ between the record cached for build 1 and the same accession.version in the pinned RefSeq release", "row")
    span = lambda f: f"{f['start']}..{f['end']}({f['strand']})"
    for k in sorted(set(before) | set(after)):
        o, n = before.get(k), after.get(k); key = f"{sid}:{k}"
        if o and not n:
            removed.append(f"genome_features:{key}"); change("genome_features", key, "(row)", span(o), None, fc)
        elif n and not o:
            change("genome_features", key, "(row)", None, span(n), fc)
        else:
            for f in ("start", "end", "strand", "type"):
                if o[f] != n[f]: change("genome_features", key, f, str(o[f]), str(n[f]), fc)
    return {"changes": changes, "aggregates": aggregates, "removed": removed}
