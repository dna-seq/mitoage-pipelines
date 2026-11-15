"""The RefSeq record of every species of the build, read from the pinned RefSeq release (237, the one build 2
read; `mitoage download refseq --from-snapshot` checks its md5), with the sequence of CONTIG-only records from
`download contig`. Every species of the build gets a row, retired ones included: a record that is not in the
release is a row with `in_release` false and the reason, which every later step reports as not assessed.

The audit reads each record's annotation as RefSeq deposits it. Where the build's own features (genome_features)
differ from it, as for an annotation the build patched, the record's row says so (`build_annotation`: same, differs
with the features on either side, or no record) and every per-record table carries that column: the audit does not
take the build's patch in place of the record, and does not hide that the two disagree.

The inputs of the whole audit are hashed here, once, into inputs.tsv: the build's tables and alignments, the
release, the CONTIG-only records and the taxonomy dump that `audit code` reads.
"""
import hashlib
from pathlib import Path

import pandas as pd
import polars as pl

from ..build2 import genome as gm
from ..build2 import pins
from .common import Warnings, write_tsv

TAXONOMY_COLUMNS = ["class", "order", "family", "genus"]


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def inputs_table(build: Path, release: Path, contig: Path, taxdump: Path) -> pd.DataFrame:
    rows = []
    for name in ("species.parquet", "region_stats.parquet", "mtdna_structure.parquet", "genome_features.parquet"):
        rows.append({"input": f"build/{name}", "source": "longevity-genie/mitoage (download build --revision)",
                     "sha256": sha256_of(build / name)})
    for p in sorted((build / "alignments").glob("*.fasta.gz")):
        if p.name.startswith("rRNA"):
            continue
        rows.append({"input": f"build/alignments/{p.name}", "source": "longevity-genie/mitoage (download build --revision)",
                     "sha256": sha256_of(p)})
    from ..build2.download import md5_of
    md5 = md5_of(release)
    pinned = "pinned md5 matches" if md5 == pins.REFSEQ_MD5 else f"md5 {md5} is not the pinned {pins.REFSEQ_MD5}"
    rows.append({"input": f"refseq/{release.parent.name}/{release.name}", "source": f"NCBI RefSeq mitochondrion ({pinned})",
                 "sha256": sha256_of(release)})
    for p in sorted(contig.glob("*.gb")):
        rows.append({"input": f"genbank/contig/{p.name}", "source": "NCBI efetch gbwithparts (download contig)", "sha256": sha256_of(p)})
    for name in ("nodes.dmp", "merged.dmp"):
        rows.append({"input": f"taxdump/{taxdump.name}/{name}", "source": "NCBI Taxonomy new_taxdump (download taxonomy)",
                     "sha256": sha256_of(taxdump / name)})
    return pd.DataFrame(rows)


def collect(build: Path, release: Path, contig: Path, refuse_unpinned: bool = True) -> pd.DataFrame:
    """One row per species of the build, with its record's GenBank text when the release has it."""
    from ..build2.download import md5_of
    from ..build2.refseq import scan
    if not (build / "species.parquet").is_file():
        raise FileNotFoundError(f"no {build / 'species.parquet'}: run `mitoage download build` first")
    if refuse_unpinned and md5_of(release) != pins.REFSEQ_MD5:
        raise ValueError(f"{release} is not RefSeq release {pins.REFSEQ_RELEASE} (md5 differs from the pin); "
                         "`mitoage download refseq --from-snapshot`, or pass --any-release to audit another release")
    sp = pl.read_parquet(build / "species.parquet").select(
        "id", "slug", "scientific_name", "refseq", "record_status", "phylum", *TAXONOMY_COLUMNS, "transl_table").to_pandas()
    want = set(sp.refseq.dropna())
    headers, texts = {}, {}
    for h, t in scan(release, keep=lambda h: h.version in want):
        if t is not None:
            headers[h.version] = h
            texts[h.version] = t
    warn = Warnings("audit records")
    rows = []
    for d in sp.to_dict("records"):
        s = pd.Series(d)
        row = {"slug": s.slug, "species_id": s.id, "scientific_name": s.scientific_name, "refseq": s.refseq,
               "record_status": s.record_status, "phylum": s.phylum, "build_transl_table": s.transl_table,
               **{c: d[c] or None for c in TAXONOMY_COLUMNS}}        # the build writes a missing order as "": missing, not a name
        h = headers.get(s.refseq)
        if h is None:
            warn.add(f"not in the release (status {s.record_status})", s.refseq)
            rows.append({**row, "in_release": False, "reason": "the record is not in the RefSeq release read",
                         "vertebrate": None, "lineage": None, "record_taxid": None, "length": None, "genbank": None})
            continue
        text = texts[s.refseq]
        if not h.has_sequence:
            p = contig / f"{s.refseq}.gb"
            if not p.is_file():
                raise FileNotFoundError(f"{s.refseq} has only a CONTIG line in the release: `mitoage download contig --from-snapshot`")
            text = p.read_text()
            warn.add("sequence read from the CONTIG-only download", s.refseq)
        rows.append({**row, "in_release": True, "reason": "", "vertebrate": "Vertebrata" in h.lineage.split("; "),
                     "lineage": h.lineage, "record_taxid": h.taxid, "length": h.length, "genbank": text})
    df = pd.DataFrame(rows)
    built = build_features(build)
    status, detail = [], []
    for r in df.itertuples(index=False):
        st, de = compare_with_build(genome_of(r) if r.in_release else None, built.get(r.slug))
        status.append(st); detail.append(de)
        if st != "same":
            warn.add(f"the build's annotation: {st}", f"{r.slug}: {de}" if de else r.slug)
    df["build_annotation"], df["build_annotation_detail"] = status, detail
    warn.emit()
    return df


def _fmt(x: tuple) -> str:
    name, typ, start, end, strand = x
    return f"{name} {typ} {start}..{end} ({'+' if strand == 1 else '-'})"


def build_features(build: Path) -> dict[str, set]:
    """slug -> the build's features as (name, type, start, end, strand)."""
    gf = pl.read_parquet(build / "genome_features.parquet").select("slug", "name", "type", "start", "end", "strand")
    out: dict[str, set] = {}
    for slug, name, typ, start, end, strand in gf.iter_rows():
        out.setdefault(slug, set()).add((name, typ, int(start), int(end), int(strand)))
    return out


def compare_with_build(g, built: set | None) -> tuple[str, str]:
    """same, differs (with what only the build and only the record carry), or no record."""
    if g is None:
        return "no record", ""
    mine = {(f.name, f.type, f.start, f.end, f.strand) for f in g.features}
    built = built or set()
    if mine == built:
        return "same", ""
    parts = []
    if built - mine:
        parts.append("only in the build: " + ", ".join(_fmt(x) for x in sorted(built - mine)))
    if mine - built:
        parts.append("only in the record: " + ", ".join(_fmt(x) for x in sorted(mine - built)))
    return "differs", "; ".join(parts)


def main(build: Path, release: Path, contig: Path, taxdump: Path, interim: Path, out: Path, refuse_unpinned: bool = True) -> pd.DataFrame:
    df = collect(build, release, contig, refuse_unpinned)
    interim.mkdir(parents=True, exist_ok=True)
    df.sort_values("slug", kind="mergesort").to_parquet(interim / "records.parquet", index=False)
    write_tsv(inputs_table(build, release, contig, taxdump), out / "inputs.tsv", ["input"])
    write_tsv(df[["slug", "refseq", "build_annotation", "build_annotation_detail"]], out / "annotation_vs_build.tsv", ["slug"])
    print(f"{len(df)} species, {int(df.in_release.sum())} records in the release -> {interim / 'records.parquet'}; "
          f"input hashes -> {out / 'inputs.tsv'}")
    return df


def load(interim: Path) -> pd.DataFrame:
    p = interim / "records.parquet"
    if not p.is_file():
        raise FileNotFoundError(f"no {p}: run `mitoage audit records` first")
    return pd.read_parquet(p)


def genome_of(row) -> gm.Genome:
    """The record annotated with build 2's rules (feature names, Leu/Ser isoacceptors by gene order in chordates)."""
    return gm.annotate(row.genbank, gm.fallback_table(row.phylum), row.phylum == "Chordata")


# --- the same record written another way: the tests' controls -----------------------------------
def _map(g: gm.Genome, seq: str, seg_map) -> gm.Genome:
    n = len(seq)
    feats = []
    for f in g.features:
        segs = [seg_map(s) for s in f.segments]
        start, end, strand = gm.feature_span(segs, n)
        feats.append(gm.Feature(f.name, f.type, start, end, strand, segs, f.origin, f.codon_start, f.via))
    feats.sort(key=lambda f: (f.start, f.end))
    raw = [gm.GbFeature(f.key, f.location, [seg_map(s) for s in f.segments], f.partial, f.qualifiers) for f in g.record.features]
    r = g.record
    rec = gm.GbRecord(r.accession_version, r.length, r.circular, r.definition, r.organism, raw, seq)
    return gm.Genome(rec, feats, g.transl_table, g.taxid, list(g.unrecognised))


def reverse_complement(g: gm.Genome) -> gm.Genome:
    """The record as it would be deposited from the other strand: the sequence reverse-complemented and every
    feature moved to the mirrored coordinates and the other strand, its segments kept in reading order."""
    n = len(g.record.sequence)
    return _map(g, gm.reverse_complement(g.record.sequence), lambda s: gm.Segment(n - s.end + 1, n - s.start + 1, -s.strand))


def rotate(g: gm.Genome, k: int) -> gm.Genome:
    """The circular record started at base k + 1 instead of base 1. No segment of a named feature or of a CDS, rRNA
    or tRNA may contain the new origin (the source feature and gene features are not read by the audit)."""
    seq = g.record.sequence
    n = len(seq)
    read = [(f.name, f.segments) for f in g.features] + [(f.key, f.segments) for f in g.record.features if f.key in ("CDS", "rRNA", "tRNA")]
    for name, segs in read:
        for s in segs:
            if s.start <= k < s.end:
                raise ValueError(f"the new origin {k + 1} falls inside {name} {s.start}..{s.end}")
    shift = lambda p: (p - 1 - k) % n + 1
    return _map(g, seq[k:] + seq[:k], lambda s: gm.Segment(shift(s.start), shift(s.end), s.strand))
