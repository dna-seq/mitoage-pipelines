"""Genetic code. The translation table a record declares (the /transl_table its protein genes state, the most
common one when they disagree) against the one NCBI Taxonomy assigns to its lineage (the mitochondrial genetic
code of the record's taxon in nodes.dmp, the taxon made current through merged.dmp).

Internal stop codons in each copy of the 13 protein genes: the coding sequence in its own sense from codon_start,
complete codons only, every codon but the last complete one, under the declared code and under the expected one.
A codon with an ambiguity code is neither a stop nor a sense codon (counted apart); a codon a /transl_except names
with an amino acid is not a stop. A record is flagged when both codes can be read and neither gives a clean
translation. A code this module does not know, or a record without one, is unknown: never a clean translation.
"""
import re
from collections import Counter
from pathlib import Path

import pandas as pd
import polars as pl

from ..build2 import genome as gm
from ..build2 import taxonomy as tx
from . import records as rec
from .common import NO, UNKNOWN, YES, Warnings, write_tsv

# NCBI's mitochondrial codes that genome.py does not carry: trematode (21), Rhabdopleuridae (24), Cephalodiscidae (33)
EXTRA = {21: {"AGA": "S", "AGG": "S", "AUA": "M", "UGA": "W", "AAA": "N"},
         24: {"AGA": "S", "AGG": "K", "UGA": "W"},
         33: {"UAA": "Y", "UGA": "W", "AGA": "S", "AGG": "K"}}


def code_table(n) -> dict[str, str] | None:
    """Codon (RNA letters) -> amino acid or *, for an NCBI table number; None for one this audit does not know."""
    if n is None or pd.isna(n):
        return None
    n = int(n)
    if n in gm.CODE_DIFFS:
        return gm.genetic_code(n)
    if n in EXTRA:
        return {c: EXTRA[n].get(c, gm.STANDARD[i]) for i, c in enumerate(gm.CODONS)}
    return None


def sense_positions(segments) -> list[int]:
    """Genomic position of every base of a feature, in its own sense."""
    out = []
    for s in segments:
        out.extend(range(s.start, s.end + 1) if s.strand == 1 else range(s.end, s.start - 1, -1))
    return out


def excepted(f: gm.GbFeature) -> set[int]:
    """Genomic positions of the codons a /transl_except gives an amino acid (not TERM)."""
    out = set()
    for q in f.qualifiers.get("transl_except", []):
        m = re.search(r"pos:(.+?),\s*aa:(\w+)", q)
        if not m or m.group(2).upper() == "TERM":
            continue
        segs, _ = gm.parse_location(m.group(1))
        out.update(sense_positions(segs))
    return out


def cds_of(g: gm.Genome) -> list[dict]:
    """Every CDS named as one of the 13 genes, in record order, with its coding sequence."""
    seq = g.record.sequence
    copies: Counter = Counter()
    out = []
    for f in g.record.features:
        if f.key != "CDS" or not f.segments:
            continue
        n = gm.feature_name(f)
        if not n or n[0] not in gm.PROTEIN_GENES:
            continue
        copies[n[0]] += 1
        try:
            cs = int(f.qualifiers.get("codon_start", ["1"])[0]) or 1
        except ValueError:
            cs = 1
        try:
            t = int(f.qualifiers.get("transl_table", [""])[0])
        except ValueError:
            t = None
        out.append({"gene": n[0], "copy": copies[n[0]], "nt": gm.extract(seq, f.segments)[cs - 1:].upper(),
                    "pos": sense_positions(f.segments)[cs - 1:], "except": excepted(f), "table": t})
    return out


def internal_stops(c: dict, table: dict[str, str] | None) -> tuple[int | None, int | None, int]:
    """(internal stops, codon number of the first, codons with an ambiguity code); stops None under an unknown code."""
    nt = c["nt"]
    n = len(nt) // 3
    stops, first, ambiguous = 0, None, 0
    for i in range(n - 1):
        codon = nt[3 * i:3 * i + 3]
        if set(codon) - set("ACGT"):
            ambiguous += 1
            continue
        if c["except"] and set(c["pos"][3 * i:3 * i + 3]) & c["except"]:
            continue
        if table is not None and table[codon.replace("T", "U")] == "*":
            stops += 1
            if first is None:
                first = i + 1
    return (None, None, ambiguous) if table is None else (stops, first, ambiguous)


def declared(cds: list[dict]) -> tuple[int | None, str]:
    """The table most protein genes state (the lower number on a tie), and every stated table with its count."""
    c = Counter(x["table"] for x in cds if x["table"] is not None)
    if not c:
        return None, ""
    best = sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
    return best, ",".join(f"{t}x{n}" for t, n in sorted(c.items()))


def answer(stops) -> str:
    return UNKNOWN if stops is None else YES if stops == 0 else NO


def neither_clean(clean_declared: str, clean_expected: str) -> str:
    """yes when both codes were read and neither is clean; no when either is clean; unknown otherwise."""
    if YES in (clean_declared, clean_expected):
        return NO
    if clean_declared == NO and clean_expected == NO:
        return YES
    return UNKNOWN


class MitoCodes:
    """The mitochondrial genetic code NCBI Taxonomy gives each taxon (nodes.dmp field 9, inherited values filled in)."""

    def __init__(self, folder: Path):
        nodes = tx._dmp(folder / "nodes.dmp", [0, 8], ["taxid", "mgc"])
        self.mgc = dict(zip(nodes["taxid"].cast(pl.Int64).to_list(), nodes["mgc"].cast(pl.Int64).to_list()))
        merged = tx._dmp(folder / "merged.dmp", [0, 1], ["old", "new"])
        self.merged = dict(zip(merged["old"].cast(pl.Int64).to_list(), merged["new"].cast(pl.Int64).to_list()))

    def of(self, taxid) -> tuple[int | None, str]:
        if taxid is None or pd.isna(taxid):
            return None, "the record names no taxon"
        t = int(taxid)
        cur = self.merged.get(t, t)
        note = f"taxid {t}" + (f" (merged into {cur})" if cur != t else "")
        m = self.mgc.get(cur)
        if m is None:
            return None, f"{note} is not in the taxonomy dump"
        if m == 0:
            return None, f"{note} has no mitochondrial code in the taxonomy dump"
        return m, note


def assess(g: gm.Genome, expected: int | None) -> tuple[dict, list[dict]]:
    cds = cds_of(g)
    dec, stated = declared(cds)
    tdec, texp = code_table(dec), code_table(expected)
    genes = []
    for c in cds:
        sd, fd, amb = internal_stops(c, tdec)
        se, fe, _ = internal_stops(c, texp)
        genes.append({"gene": c["gene"], "copy": c["copy"], "nt_length": len(c["nt"]), "whole_codons": YES if len(c["nt"]) % 3 == 0 else NO,
                      "codons": len(c["nt"]) // 3, "ambiguous_codons": amb, "stops_declared": sd, "first_stop_declared": fd,
                      "stops_expected": se, "first_stop_expected": fe})
    reasons = []
    if not cds:
        reasons.append("no protein gene annotated")
    if dec is None and cds:
        reasons.append("no transl_table stated")
    elif dec is not None and tdec is None:
        reasons.append(f"declared table {dec} not known to the audit")
    if expected is not None and texp is None:
        reasons.append(f"expected table {expected} not known to the audit")
    sum_or_none = lambda k: None if not cds or any(x[k] is None for x in genes) else sum(x[k] for x in genes)
    sd, se = sum_or_none("stops_declared"), sum_or_none("stops_expected")
    with_stops = lambda k: ",".join(f"{x['gene']}{'#' + str(x['copy']) if x['copy'] > 1 else ''}:{x[k]}" for x in genes if x[k])
    cd, ce = answer(sd), answer(se)
    row = {"declared_table": dec, "stated_tables": stated,
           "code_match": UNKNOWN if dec is None or expected is None else YES if dec == expected else NO,
           "stops_declared": sd, "stops_expected": se, "genes_with_stops_declared": with_stops("stops_declared"),
           "genes_with_stops_expected": with_stops("stops_expected"), "clean_declared": cd, "clean_expected": ce,
           "neither_clean": neither_clean(cd, ce), "reason": "; ".join(reasons)}
    return row, genes


def main(interim: Path, taxdump: Path, out: Path) -> pd.DataFrame:
    df = rec.load(interim).rename(columns={"class": "class_"})
    codes = MitoCodes(taxdump)
    warn = Warnings("audit code")
    rows, gene_rows = [], []
    for r in df.itertuples(index=False):
        base = {"slug": r.slug, "refseq": r.refseq, "class": r.class_, "order": r.order, "family": r.family,
                "build_transl_table": r.build_transl_table}
        exp, note = codes.of(r.record_taxid if r.in_release else None)
        if not r.in_release:
            warn.add("no record", r.slug)
            rows.append({**base, "expected_table": None, "expected_from": "", "code_match": UNKNOWN, "clean_declared": UNKNOWN,
                         "clean_expected": UNKNOWN, "neither_clean": UNKNOWN, "reason": r.reason})
            continue
        if exp is None:
            warn.add(note.split(" ", 2)[-1] if note.startswith("taxid") else note, r.slug)
        row, genes = assess(rec.genome_of(r), exp)
        if row["reason"]:
            for x in row["reason"].split("; "):
                warn.add(x, r.slug)
        rows.append({**base, "expected_table": exp, "expected_from": note, **row})
        gene_rows += [{"slug": r.slug, **x} for x in genes]
    warn.emit()
    cols = ["slug", "refseq", "class", "order", "family", "declared_table", "stated_tables", "build_transl_table", "expected_table",
            "expected_from", "code_match", "stops_declared", "stops_expected", "clean_declared", "clean_expected", "neither_clean",
            "genes_with_stops_declared", "genes_with_stops_expected", "reason"]
    tab = pd.DataFrame(rows).reindex(columns=cols)
    for c in ("declared_table", "build_transl_table", "expected_table", "stops_declared", "stops_expected"):
        tab[c] = tab[c].astype("Int64")
    genes = pd.DataFrame(gene_rows)
    for c in ("stops_declared", "first_stop_declared", "stops_expected", "first_stop_expected"):
        genes[c] = genes[c].astype("Int64")
    tab = tab.merge(df[["slug", "build_annotation"]], on="slug", how="left", validate="one_to_one")
    write_tsv(tab, out / "code.tsv", ["slug"])
    write_tsv(genes, out / "code_genes.tsv", ["slug", "gene", "copy"])
    print(f"code match: {tab.code_match.value_counts().sort_index().to_dict()}; neither code clean: "
          f"{tab.neither_clean.value_counts().sort_index().to_dict()}")
    return tab
