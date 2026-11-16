"""The sequencing technology behind each mitogenome of the build (the control H8 could not run).

A RefSeq mitochondrial record rarely states how it was sequenced; the INSDC record it was copied
from sometimes does, in the structured comment GenBank calls Assembly-Data:

    ##Assembly-Data-START##
    Assembly Method       :: MitoHiFi v. 2.2
    Sequencing Technology :: PacBio Sequel
    ##Assembly-Data-END##

`fetch` downloads, by E-utilities, the RefSeq record of every species of the build and then the
INSDC record each one names in its COMMENT ("derived from X", "identical to X"). `main` reads both
and writes technology.tsv, one row per species:

    slug, refseq, source, technology, stated_in, technology_class

`technology` is the text as the record gives it; `stated_in` says which record gave it (the INSDC
source first, since it is the submitter's own; the RefSeq record otherwise); `technology_class` is
long_read, short_read, sanger or unstated. A record that states nothing is `unstated`: it is never
guessed from its date or its journal, and never counted as short-read.

A record naming both a long-read and a short-read platform (a hybrid or polished assembly) is
long_read: the repeats and control regions the control is about were spanned by long reads.
"""
import re
from pathlib import Path

import polars as pl

from ..build2.assemble import primary_accession
from ..build2.download import efetch

# word-bounded, lower case; "ont" alone would match inside other words
LONG_READ = re.compile(r"pacbio|pac bio|pacific biosciences|nanopore|\bont\b|hifi|\bccs\b|\bclr\b|minion|promethion|gridion|sequel|revio|\brs ?ii\b|\bsmrt\b")
SHORT_READ = re.compile(r"illumina|solexa|hiseq|miseq|novaseq|nextseq|\b454\b|roche|ion ?torrent|ion ?proton|ion pgm|\bsolid\b|bgi|dnbseq|mgi|genome analyzer|x ?ten\b|complete genomics|\bflx\b")
SANGER = re.compile(r"sanger|dideoxy|\babi\b|abi ?3\d{3}|\b3730|\b3130|\b3100|\b377\b|capillary|megabace|li-?cor|primer walking")

CLASSES = ("long_read", "short_read", "sanger", "unstated")


def stated_technology(record: str) -> str | None:
    """The `Sequencing Technology ::` value of a GenBank record's COMMENT, with line wraps joined; None when absent."""
    m = re.search(r"\nCOMMENT(.*?)\nFEATURES", record, flags=re.S)
    if not m:
        return None
    t = re.search(r"Sequencing Technology[ \t]*::[ \t]*(.*?)(?=\n\s*[A-Za-z][^\n]*::|\n\s*##|\Z)", m.group(1), flags=re.S)
    if not t:
        return None
    value = re.sub(r"\s+", " ", t.group(1)).strip()
    return value or None


def classify(technology: str | None) -> str:
    """long_read when any long-read platform is named; else short_read, else sanger; unstated for nothing or for a text naming no known platform."""
    if not technology:
        return "unstated"
    t = technology.lower()
    if LONG_READ.search(t):
        return "long_read"
    if SHORT_READ.search(t):
        return "short_read"
    if SANGER.search(t):
        return "sanger"
    return "unstated"


def _record(folder: Path, accession: str) -> str | None:
    """The stored record of an accession, with or without a version in the name asked for."""
    exact = folder / f"{accession}.gb"
    if exact.is_file():
        return exact.read_text()
    hits = sorted(folder.glob(f"{accession.split('.')[0]}.*.gb"))
    return hits[-1].read_text() if hits else None


def fetch(data: Path, folder: Path) -> None:
    """The RefSeq record of each species of the build into folder/refseq, then the INSDC records they name into folder/source."""
    species = pl.read_parquet(data / "species.parquet", columns=["slug", "refseq"]).drop_nulls("refseq")
    versions = species["refseq"].to_list()
    missing = efetch(versions, folder / "refseq")
    if missing:
        print(f"NCBI returned no record for {len(missing)} of {len(versions)} RefSeq versions: {', '.join(missing[:10])}")
    sources = []
    for v in versions:
        rec = _record(folder / "refseq", v)
        src = primary_accession(rec) if rec else None
        if src:
            sources.append(src)
    print(f"{len(sources)} of {len(versions)} RefSeq records name an INSDC source")
    missing = efetch(sources, folder / "source")
    if missing:
        print(f"NCBI returned no record for {len(missing)} of {len(set(sources))} INSDC sources: {', '.join(missing[:10])}")


def table(data: Path, folder: Path) -> pl.DataFrame:
    species = pl.read_parquet(data / "species.parquet", columns=["slug", "refseq"])
    rows = []
    for slug, refseq in species.iter_rows():
        rec = _record(folder / "refseq", refseq) if refseq else None
        src = primary_accession(rec) if rec else None
        src_rec = _record(folder / "source", src) if src else None
        from_source = stated_technology(src_rec) if src_rec else None
        from_refseq = stated_technology(rec) if rec else None
        tech = from_source or from_refseq
        rows.append({
            "slug": slug, "refseq": refseq, "source": src,
            "refseq_record": rec is not None, "source_record": src_rec is not None,
            "technology": tech,
            "stated_in": "source" if from_source else "refseq" if from_refseq else None,
            "technology_class": classify(tech),
        })
    return pl.DataFrame(rows, schema_overrides={"technology": pl.String, "stated_in": pl.String, "source": pl.String}).sort("slug")


def main(data: Path, folder: Path, out: Path) -> pl.DataFrame:
    """Writes technology.tsv and prints what could and could not be read."""
    if not (folder / "refseq").is_dir():
        raise SystemExit(f"{folder / 'refseq'} is missing: run `mitoage download technology` first")
    t = table(data, folder)
    out.mkdir(parents=True, exist_ok=True)
    t.write_csv(out / "technology.tsv", separator="\t")
    n = len(t)
    print(f"{n} species; RefSeq record read for {t['refseq_record'].sum()}, an INSDC source named by "
          f"{t['source'].is_not_null().sum()}, its record read for {t['source_record'].sum()}")
    counts = dict(t.group_by("technology_class").len().iter_rows())
    for c in CLASSES:
        print(f"  {c}: {counts.get(c, 0)}")
    odd = t.filter(pl.col("technology").is_not_null() & (pl.col("technology_class") == "unstated"))
    if len(odd):
        print(f"  {len(odd)} records state a technology that names no known platform (kept as unstated):")
        for text, k in odd.group_by("technology").len().sort("len", descending=True).head(15).iter_rows():
            print(f"    {k} x {text}")
    print(f"written {out / 'technology.tsv'}")
    return t
