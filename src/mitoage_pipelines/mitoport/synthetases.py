"""The nuclear aminoacyl-tRNA synthetases that charge the mitochondrial tRNAs, resolved from NCBI at run time.

H37 part C's entry requires the human symbols to be read from NCBI, not from the entry, and the genes it cannot
resolve to be named as "not assessed". So nothing here holds a list of symbols. The resolution is:

1. Every live human gene NCBI annotates with a tRNA ligase activity (the esearch term below). This is wider than
   the mitochondrial set on purpose: a search for the mitochondrial ones by their own words misses them. The
   narrower terms were tried and are recorded here because the next reader will try them too: the description
   "tRNA synthetase 2, mitochondrial" is not indexed in `[Title]` and returns nothing, "mitochondrial
   aminoacyl-tRNA synthetase" returns 6 of them, and "aminoacyl-tRNA ligase activity" returns 38 and drops HARS2.
2. Each gene's record is read for its official description and its GO terms, by category.
3. A gene is `mitochondrial` when NCBI's own description ends in ", mitochondrial"; `dual` when it does not but the
   gene carries a mitochondrial cellular component *and* a GO process naming a mitochondrial aminoacylation; and
   `ambiguous` when it carries the component with no such process. Only the first two are used. An ambiguous gene
   is written out with its evidence and counted as not assessed, never quietly kept or quietly dropped: a
   mitochondrial component alone is also what a cytosolic synthetase with one stray annotation looks like.

The amino acid is taken from the description's prefix, not from the GO process terms: NCBI's process annotation for
DARS2 (aspartyl-tRNA synthetase 2) reads "mitochondrial asparaginyl-tRNA aminoacylation", so the process terms name
the wrong amino acid for at least one gene of the set, while the description and the symbol agree.

An amino acid left with no gene is a result, not a gap to fill by hand. Two are expected and are reported as not
assessed rather than guessed: glutamine, which has no direct mitochondrial synthetase (the mitochondrial Gln-tRNA
is charged by the GatCAB amidotransferase, which is not a synthetase and is not in the search), and any amino acid
whose only candidate came out `ambiguous`.

Writes synthetases.tsv: amino_acid, code, symbol, gene_id, status, description, components, processes.
"""
import re
from pathlib import Path
from xml.etree import ElementTree as ET

from ..download import EUTILS, _request

# the search, fixed here and not narrowed to the mitochondrial ones: see the docstring
TERM = '"Homo sapiens"[orgn] AND "tRNA ligase activity"[All Fields] AND alive[prop]'
RETMAX = 500
# the adjectival form NCBI's descriptions use -> the amino acid. A naming table, not a gene list: every symbol,
# every gene id and every localisation comes from the records themselves.
ACYL = {
    "alanyl": ("alanine", "Ala"), "arginyl": ("arginine", "Arg"), "asparaginyl": ("asparagine", "Asn"),
    "aspartyl": ("aspartate", "Asp"), "cysteinyl": ("cysteine", "Cys"), "glutaminyl": ("glutamine", "Gln"),
    "glutamyl": ("glutamate", "Glu"), "glycyl": ("glycine", "Gly"), "histidyl": ("histidine", "His"),
    "isoleucyl": ("isoleucine", "Ile"), "leucyl": ("leucine", "Leu"), "lysyl": ("lysine", "Lys"),
    "methionyl": ("methionine", "Met"), "phenylalanyl": ("phenylalanine", "Phe"), "prolyl": ("proline", "Pro"),
    "seryl": ("serine", "Ser"), "threonyl": ("threonine", "Thr"), "tryptophanyl": ("tryptophan", "Trp"),
    "tyrosyl": ("tyrosine", "Tyr"), "valyl": ("valine", "Val"),
}
MITO_COMPONENTS = ("mitochondrion", "mitochondrial matrix", "mitochondrial membrane", "mitochondrial inner membrane")
COLUMNS = ["amino_acid", "code", "symbol", "gene_id", "status", "description", "components", "processes"]
USED = ("mitochondrial", "dual")          # the statuses part C tests; everything else is not assessed


# --- what NCBI answers: parsed, with no network -----------------------------------------------------------------------


def amino_acid(description: str) -> tuple[str, str] | None:
    """(amino acid, three-letter code) named by the prefix of an NCBI gene description, or None when the
    description names none or names more than one (`glutamyl-prolyl-tRNA synthetase 1`, a fused cytosolic pair)."""
    m = re.match(r"^([a-z]+yl)[- ]tRNA synthetase\b", description.strip(), flags=re.IGNORECASE)
    if m is None:
        return None
    return ACYL.get(m.group(1).lower())


def go_terms(gene: ET.Element) -> dict[str, list[str]]:
    """{"Function"|"Process"|"Component": the gene's GO terms of that category, each once, sorted}."""
    out: dict[str, set[str]] = {}
    for c in gene.iter("Gene-commentary"):
        if (c.findtext("Gene-commentary_heading") or "") != "GeneOntology":
            continue
        comment = c.find("Gene-commentary_comment")
        if comment is None:
            continue
        for sub in comment:
            label = sub.findtext("Gene-commentary_label") or ""
            terms = {a.text for a in sub.iter("Other-source_anchor") if a.text}
            out.setdefault(label, set()).update(terms)
    return {k: sorted(v) for k, v in sorted(out.items())}


def status(description: str, components: list[str], processes: list[str]) -> str:
    """`mitochondrial`, `dual`, `ambiguous` or `other`, by the rule in the docstring. A gene with no mitochondrial
    cellular component at all is `other` whatever its description says."""
    mito_component = any(t in MITO_COMPONENTS for t in components)
    if description.strip().lower().endswith(", mitochondrial"):
        return "mitochondrial"
    if not mito_component:
        return "other"
    mito_process = any(p.lower().startswith("mitochondrial ") and "aminoacylation" in p.lower() for p in processes)
    return "dual" if mito_process else "ambiguous"


def parse(xml: bytes) -> list[dict]:
    """One row per gene record in an efetch answer, in the answer's order, whatever its status."""
    rows = []
    for g in ET.fromstring(xml).iter("Entrezgene"):
        symbol = g.findtext(".//Gene-ref/Gene-ref_locus") or ""
        gene_id = g.findtext(".//Gene-track/Gene-track_geneid") or ""
        description = g.findtext(".//Gene-ref/Gene-ref_desc") or ""
        go = go_terms(g)
        components, processes = go.get("Component", []), go.get("Process", [])
        aa = amino_acid(description)
        rows.append({
            "amino_acid": aa[0] if aa else "", "code": aa[1] if aa else "", "symbol": symbol, "gene_id": gene_id,
            "status": status(description, components, processes), "description": description,
            "components": "; ".join(components), "processes": "; ".join(processes),
        })
    return rows


def resolve(rows: list[dict]) -> tuple[dict[str, dict], list[dict]]:
    """({code: the one gene part C uses for that amino acid}, the rows it does not use).

    Refused when two used genes claim one amino acid, or when a used gene names no amino acid: either is a set
    this code cannot resolve, and the entry's answer to that is to stop, not to choose.
    """
    used, rest = {}, []
    for r in rows:
        if r["status"] not in USED:
            rest.append(r); continue
        if not r["code"]:
            raise ValueError(f"{r['symbol']} ({r['status']}) names no amino acid in its description: {r['description']!r}")
        if r["code"] in used:
            raise ValueError(f"two genes for {r['code']}: {used[r['code']]['symbol']} and {r['symbol']}")
        used[r["code"]] = r
    return used, rest


def not_assessed(used: dict[str, dict], rest: list[dict]) -> list[tuple[str, str, str]]:
    """(code, amino acid, why) for each of the 20 amino acids part C has no gene for, sorted by code."""
    out = []
    for acyl, (name, code) in sorted(ACYL.items(), key=lambda kv: kv[1][1]):
        if code in used:
            continue
        others = [r for r in rest if r["code"] == code]
        if not others:
            why = "no human gene with a tRNA ligase activity names this amino acid in its description"
        else:
            why = "; ".join(f"{r['symbol']} is {r['status']}" for r in sorted(others, key=lambda r: r["symbol"]))
        out.append((code, name, why))
    return out


# --- the run ------------------------------------------------------------------------------------------------------


def search() -> list[str]:
    """The gene ids the search term returns, in NCBI's order."""
    import json
    from urllib.parse import quote
    url = f"{EUTILS}/esearch.fcgi?db=gene&retmode=json&retmax={RETMAX}&term={quote(TERM)}"
    answer = json.loads(_request(url))["esearchresult"]
    ids = answer["idlist"]
    if int(answer["count"]) > len(ids):
        raise ValueError(f"the search returned {answer['count']} genes and only {len(ids)} were read: raise RETMAX")
    return ids


def fetch(ids: list[str]) -> bytes:
    """The full gene records of `ids`, as efetch XML."""
    return _request(f"{EUTILS}/efetch.fcgi?db=gene&retmode=xml&id={','.join(ids)}", accept="*/*")


def main(out: Path, xml: Path | None = None) -> None:
    """Write the table to `out`. `xml` is an efetch answer already on disk, used instead of the network when given
    (a run of a recorded answer); when given and absent, the answer is fetched and written there, so a run can be
    repeated on the same bytes. The amino acids with no gene are printed, not dropped."""
    out.parent.mkdir(parents=True, exist_ok=True)
    if xml is not None and xml.is_file():
        body = xml.read_bytes()
        print(f"the records of {xml} ({len(body)} bytes), not the network")
    else:
        ids = search()
        print(f"{len(ids)} human genes with a tRNA ligase activity")
        body = fetch(ids)
        if xml is not None:
            xml.parent.mkdir(parents=True, exist_ok=True); xml.write_bytes(body)
    rows = parse(body)
    used, rest = resolve(rows)
    missing = not_assessed(used, rest)
    lines = ["\t".join(COLUMNS)]
    for code in sorted(used):
        lines.append("\t".join(str(used[code][c]) for c in COLUMNS))
    for code, name, why in missing:
        lines.append("\t".join([name, code, "", "", "not assessed", why, "", ""]))
    out.write_text("\n".join(lines) + "\n")
    kinds = {}
    for r in rows:
        kinds[r["status"]] = kinds.get(r["status"], 0) + 1
    print(f"{len(rows)} records read: " + ", ".join(f"{v} {k}" for k, v in sorted(kinds.items())))
    print(f"{len(used)} of the 20 amino acids have a gene: " + " ".join(f"{c}={used[c]['symbol']}" for c in sorted(used)))
    for code, name, why in missing:
        print(f"not assessed: {code} ({name}) — {why}")
    print(f"wrote {out}")
