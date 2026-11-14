"""The RefSeq mitochondrion release -> one complete mitogenome per NCBI species (site repository:
pipeline/sources/refseq-release.ts, pipeline/steps/refseq.ts).

Complete: an NC_ record whose definition is not partial, nearly complete or one of several molecules
(chromosome, minicircle, segment, part, fragment, minichromosome, molecule), and says complete genome /
complete mitochondrial genome, complete sequence, or is a circular organelle assembly. Metazoa only (by
the record's lineage). Several complete records of one species: newest LOCUS date first, then full length
(not shorter than 75 % of the median of complete records of its NCBI class), then the higher
accession.version. A build 1 species whose record left RefSeq takes a live complete record of the same
species (one annotated at species rank before a subspecies one), never one another species cites.
"""
import gzip
import re
from dataclasses import dataclass
from functools import cmp_to_key
from pathlib import Path
from statistics import median

MONTHS = {m: f"{i + 1:02d}" for i, m in enumerate("JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC".split())}
MULTI = re.compile(r"chromosome|minicircle|segment|\bpart\b|fragment|minichromosome|molecule")
STATUS = re.compile(r"(REVIEWED|PROVISIONAL|VALIDATED|INFERRED|PREDICTED|MODEL|WGS) REFSEQ")


@dataclass
class Header:
    version: str
    accession: str
    length: int
    topology: str
    update_date: str
    definition: str
    organism: str
    lineage: str
    taxid: int | None
    has_sequence: bool


def locus_date(s: str) -> str:
    m = re.match(r"^(\d{2})-([A-Z]{3})-(\d{4})$", s)
    if not m or m.group(2) not in MONTHS: raise ValueError(f"bad LOCUS date {s}")
    return f"{m.group(3)}-{MONTHS[m.group(2)]}-{m.group(1)}"


def parse_header(rec: str) -> Header:
    end = rec.find("\nFEATURES")
    lines = (rec[:end] if end >= 0 else rec).split("\n")
    locus = lines[0].strip().split()
    if locus[0] != "LOCUS": raise ValueError(f"record does not start with LOCUS: {lines[0]}")
    ln = re.search(r"(\d+)\s+bp", lines[0])
    definition = version = organism = lineage = ""
    i = 0
    while i < len(lines):
        l = lines[i]
        if l.startswith("DEFINITION"):
            definition = l[12:].strip()
            while i + 1 < len(lines) and lines[i + 1].startswith(" " * 12):
                i += 1; definition += " " + lines[i].strip()
        elif l.startswith("VERSION"):
            parts = l[12:].strip().split(); version = parts[0] if parts else ""
        elif l.startswith("  ORGANISM"):
            organism = l[12:].strip(); parts = []
            while i + 1 < len(lines) and lines[i + 1].startswith(" " * 12):
                i += 1; parts.append(lines[i].strip())
            lineage = re.sub(r"\.$", "", " ".join(parts))
        i += 1
    tx = re.search(r'/db_xref="taxon:(\d+)"', rec)
    return Header(version, version.split(".")[0], int(ln.group(1)) if ln else 0,
                  "circular" if re.search(r"\bcircular\b", lines[0]) else "linear", locus_date(locus[-1]),
                  definition, organism, lineage, int(tx.group(1)) if tx else None, "\nORIGIN" in rec)


def completeness(h: Header) -> str:
    d = h.definition.lower()
    if not h.accession.startswith("NC_"): return "excluded_non_NC_contig"
    if "partial" in d or "nearly complete" in d: return "excluded_partial"
    if MULTI.search(d): return "excluded_multipartite"
    if "complete genome" in d or "complete mitochondrial genome" in d: return "complete_genome"
    if "complete sequence" in d: return "complete_sequence"
    if "genome assembly, organelle: mitochondrion" in d:
        return "complete_assembly_organelle" if h.topology == "circular" else "excluded_assembly_linear"
    return "excluded_unclassified"


def is_complete(c: str) -> bool:
    return c.startswith("complete")


def scan(path: Path, keep=None):
    """Every record of a gzipped GenBank file: yields (header, text if keep(header) else None)."""
    if not path.is_file():
        raise FileNotFoundError(f"no {path}: run `mitoage download refseq` first")
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        buf = ""
        while True:
            chunk = f.read(1 << 24)
            if not chunk: break
            buf += chunk
            start = 0
            while (i := buf.find("\n//\n", start)) >= 0:
                rec = buf[start:i] + "\n//\n"; start = i + 4
                h = parse_header(rec)
                yield h, (rec if keep is not None and keep(h) else None)
            buf = buf[start:]
    if buf.strip():
        raise ValueError(f"{path} ends with an unterminated record ({len(buf)} characters)")


@dataclass
class ReleaseRecord:
    header: Header
    completeness: str
    taxid_current: int | None
    species_taxid: int | None
    species_class: str | None
    length_qc: str = "ok"


def classify(headers: list[Header], tax) -> list[ReleaseRecord]:
    recs = []
    for h in headers:
        if not re.search(r"\bMetazoa\b", h.lineage): continue
        cur = None if h.taxid is None else tax.current(h.taxid)
        sp = None if cur is None else tax.species_of(cur)
        cls = None if sp is None else ((tax.lineage(sp) or {}).get("class") or "")
        recs.append(ReleaseRecord(h, completeness(h), cur, sp, cls))
    lengths: dict[str, list[int]] = {}
    for r in recs:
        if r.species_class is None or not is_complete(r.completeness): continue
        lengths.setdefault(r.species_class, []).append(r.header.length)
    med = {c: median(sorted(l)) for c, l in lengths.items()}
    for r in recs:
        m = None if r.species_class is None else med.get(r.species_class)
        if m is not None and r.header.length < 0.75 * m: r.length_qc = "short_for_class"
    return recs


def compare_records(a: ReleaseRecord, b: ReleaseRecord) -> int:
    """Newest LOCUS date, then full length, then the higher accession.version (string order, as the site's code)."""
    if a.header.update_date != b.header.update_date: return 1 if a.header.update_date < b.header.update_date else -1
    if a.length_qc != b.length_qc: return -1 if a.length_qc == "ok" else 1
    if a.header.version != b.header.version: return 1 if a.header.version < b.header.version else -1
    return 0


def one_per_species(recs: list[ReleaseRecord]) -> dict[int, list[ReleaseRecord]]:
    """species taxid -> its complete records, best first (the chosen one is [0])."""
    by: dict[int, list[ReleaseRecord]] = {}
    for r in recs:
        if is_complete(r.completeness) and r.species_taxid is not None:
            by.setdefault(r.species_taxid, []).append(r)
    return {t: sorted(by[t], key=cmp_to_key(compare_records)) for t in sorted(by)}


def successor_for(records: list[ReleaseRecord] | None, taken: set[str]) -> ReleaseRecord | None:
    if not records: return None
    free = [r for r in records if r.header.version not in taken]
    species_level = [r for r in free if r.taxid_current == r.species_taxid]
    pool = sorted(species_level or free, key=cmp_to_key(compare_records))
    return pool[0] if pool else None
