"""H37 part A, step 1: the 22 mitochondrial tRNA genes and the two rRNA genes, one alignment per gene, with what
each species carries of them.

Extraction. Each live species' RefSeq record (release 237, md5-checked, the one build 2 read; the CONTIG-only
records from `download contig`) is annotated by `build2.genome.annotate`, the rules that made the build's
`genome_features`: tRNA-Leu and tRNA-Ser copies whose isoacceptor the annotation does not give are placed by the
vertebrate gene order in chordates only, and one it cannot place keeps the plain name `tRNA-Leu` or `tRNA-Ser`
(`unresolved`). A gene is the first copy of its name in genome order, in its own sense (a minus-strand feature
reverse-complemented). A gene the record does not annotate, an `unresolved` Leu or Ser copy and a duplicated
gene are counted per gene in `gene_census.tsv`, never dropped silently.

Alignment. One alignment per tRNA gene, by the rule of the build's nucleotide alignments (`build2/align.py`):
MAFFT FFT-NS-2 (`--retree 2 --maxiterate 0 --quiet`), species sorted by slug, v7.526 pinned; the SHA-256 of
each aligned FASTA and the aligner's version are written to `alignments.tsv`. The two rRNA genes are NOT
realigned: the build's `rRNA-12S.nt.fasta.gz` and `rRNA-16S.nt.fasta.gz` were made by the same rule from the
same records, and the sequences extracted here are checked against them (`rrna_check.tsv`).

Per gene and species (`rna_features.tsv`): length, GC content and the two skews as the build's region
statistics compute them (`genome.make_stats` on the gene's bases, in gene sense), and the shares of A, C, G and
T among the resolved bases, the analogue of the fourfold shares H1 reads. The gene's informative alignment
columns (`informative_columns`) are counted for the whole alignment.
"""
import hashlib
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import polars as pl

from .. import external
from ..build2 import align as b2align
from ..build2 import genome as gm
from ..build2 import pins
from ..seqio import read_fasta

TRNA_GENES = ["tRNA-Ala", "tRNA-Arg", "tRNA-Asn", "tRNA-Asp", "tRNA-Cys", "tRNA-Gln", "tRNA-Glu", "tRNA-Gly",
              "tRNA-His", "tRNA-Ile", "tRNA-Leu(CUN)", "tRNA-Leu(UUR)", "tRNA-Lys", "tRNA-Met", "tRNA-Phe",
              "tRNA-Pro", "tRNA-Ser(AGY)", "tRNA-Ser(UCN)", "tRNA-Thr", "tRNA-Trp", "tRNA-Tyr", "tRNA-Val"]
RRNA_GENES = list(b2align.RRNA_GENES)
RNA_GENES = TRNA_GENES + RRNA_GENES
PLAIN_ISO = {"tRNA-Leu": "Leu", "tRNA-Ser": "Ser"}        # the names an `unresolved` copy keeps
BASES = "ACGT"
MIN_ALIGNED = 3                                            # fewer sequences than this are not aligned
INFORMATIVE_OCCUPANCY = 0.5                                # a column needs at least this share of non-gap rows
WEAK_COLUMNS = 20                                          # fewer informative columns than this: a null is uninformative


def file_name(gene: str) -> str:
    """tRNA-Leu(CUN) -> tRNA-Leu_CUN: a gene's name as a file name."""
    return gene.replace("(", "_").replace(")", "")


# --- the records -----------------------------------------------------------------------------------
def load_records(build: Path, release: Path, contig: Path) -> list[dict]:
    """The RefSeq record of every live species of the build: slug, refseq, phylum, the GenBank text. The release must be
    the pinned one; a record that is not in it, or has only a CONTIG line with no download, is refused: a species
    never silently loses its record here."""
    from ..build2.download import md5_of
    from ..build2.refseq import scan
    if md5_of(release) != pins.REFSEQ_MD5:
        raise ValueError(f"{release} is not RefSeq release {pins.REFSEQ_RELEASE} (md5 differs from the pin)")
    sp = pl.read_parquet(build / "species.parquet").select("slug", "refseq", "phylum", "record_status").to_pandas()
    sp = sp[sp.record_status.ne("retired")]
    want = set(sp.refseq.dropna())
    headers, texts = {}, {}
    for h, t in scan(release, keep=lambda h: h.version in want):
        if t is not None:
            headers[h.version], texts[h.version] = h, t
    out, missing = [], []
    for r in sp.sort_values("slug").itertuples(index=False):
        h = headers.get(r.refseq)
        if h is None:
            missing.append(f"{r.slug} ({r.refseq})")
            continue
        text = texts[r.refseq]
        if not h.has_sequence:
            p = contig / f"{r.refseq}.gb"
            if not p.is_file():
                raise FileNotFoundError(f"{r.refseq} has only a CONTIG line in the release: `mitoage download contig --from-snapshot`")
            text = p.read_text()
        out.append({"slug": r.slug, "refseq": r.refseq, "phylum": r.phylum, "text": text})
    if missing:
        raise ValueError(f"{len(missing)} live species have a record that is not in the release: {', '.join(missing[:5])}")
    return out


# --- one record ------------------------------------------------------------------------------------
def extract_record(text: str, phylum: str) -> dict:
    """The 22 tRNA and 2 rRNA genes of one record, as `genome.annotate` names them.

    -> {"sequences": gene -> first copy in genome order, gene sense, upper case;
        "copies": gene -> number of features of that name;
        "unresolved": {"Leu": n, "Ser": n}: tRNA copies of that amino acid whose isoacceptor could not be resolved}"""
    g = gm.annotate(text, gm.fallback_table(phylum), phylum == "Chordata")
    seqs, copies, unresolved = {}, {}, {}
    for f in g.features:
        if f.type not in ("tRNA", "rRNA"):
            continue
        if f.name in PLAIN_ISO:
            aa = PLAIN_ISO[f.name]
            unresolved[aa] = unresolved.get(aa, 0) + 1
            continue
        if f.name not in RNA_GENES:                              # genome.py names only these, so this is a guard
            raise ValueError(f"{f.name} is a tRNA or rRNA feature that is none of the 24 genes")
        copies[f.name] = copies.get(f.name, 0) + 1
        if f.name not in seqs:                                   # the first copy in genome order
            seqs[f.name] = gm.extract(g.record.sequence, f.segments).upper()
    return {"sequences": seqs, "copies": copies, "unresolved": unresolved}


def census(records: dict[str, dict]) -> pd.DataFrame:
    """Per gene, over the species: how many carry it and why the others do not. `records`: slug -> extract_record().

    n_species          species considered
    n_gene             species with the gene (the first copy, if duplicated)
    n_not_annotated    species whose record has no feature of that name and no unresolved copy of its amino acid
    n_unresolved       species without the gene whose record has an unresolved copy of its amino acid (Leu, Ser)
    n_duplicated       species with the gene annotated more than once (the first copy is used)
    n_unresolved_copies  unresolved copies of the gene's amino acid over all species (Leu and Ser genes only)
    """
    rows = []
    for gene in RNA_GENES:
        aa = gene.split("(")[0][5:] if gene.startswith("tRNA") and "(" in gene else None
        has = [s for s, r in records.items() if gene in r["sequences"]]
        lacking = [s for s in records if s not in set(has)]
        unres = [s for s in lacking if aa and records[s]["unresolved"].get(aa, 0) > 0]
        rows.append({"gene": gene, "n_species": len(records), "n_gene": len(has),
                     "n_not_annotated": len(lacking) - len(unres), "n_unresolved": len(unres),
                     "n_duplicated": sum(1 for r in records.values() if r["copies"].get(gene, 0) > 1),
                     "n_unresolved_copies": sum(r["unresolved"].get(aa, 0) for r in records.values()) if aa else 0})
    out = pd.DataFrame(rows)
    assert (out.n_gene + out.n_not_annotated + out.n_unresolved == out.n_species).all()
    return out


# --- composition and columns ------------------------------------------------------------------------
def gene_features(seq: str) -> dict:
    """Length, GC content, GC skew, AT skew (the build's region statistics, on the gene's bases in gene sense) and the
    shares of A, C, G and T among the resolved bases. Gaps are not bases."""
    seq = seq.upper().replace("-", "")
    b = gm.count_bases(seq)
    st = gm.make_stats(0, "gene", b, "genbank")
    n = b["a"] + b["c"] + b["g"] + b["t"]
    return {"length": len(seq), "n_other": b["other"], "gc_percent": st["gc_percent"], "gc_skew": st["gc_skew"],
            "at_skew": st["at_skew"], **{f"{x}_share": (b[x.lower()] / n if n else None) for x in BASES}}


def informative_columns(aligned: dict[str, str], subset: set[str] | None = None) -> int:
    """Columns of an alignment, over the rows in `subset` (all when None), that are occupied (at least half of the rows
    have a base) and parsimony-informative (two different bases each in at least two rows)."""
    rows = [s.upper() for k, s in aligned.items() if subset is None or k in subset]
    n = 0
    for col in zip(*rows):
        counts: dict[str, int] = {}
        for ch in col:
            if ch in BASES:
                counts[ch] = counts.get(ch, 0) + 1
        if sum(counts.values()) < INFORMATIVE_OCCUPANCY * len(rows):
            continue
        if sum(1 for v in counts.values() if v >= 2) >= 2:
            n += 1
    return n


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


# --- alignment -------------------------------------------------------------------------------------
def _align_one(args) -> tuple[str, str]:
    binary, gene, seqs = args
    rows = [(slug, "", seq) for slug, seq in sorted(seqs.items())]
    if len(rows) < MIN_ALIGNED:
        return gene, ""
    aligned = b2align.align(binary, rows, 1)
    return gene, b2align.fasta(aligned)


def align_genes(sequences: dict[str, dict[str, str]], workers: int) -> tuple[dict[str, str], str]:
    """gene -> aligned FASTA text, one MAFFT process (one thread) per gene, and MAFFT's version. Genes with fewer than
    MIN_ALIGNED sequences give an empty text."""
    binary = external.find("mafft")
    version = b2align.mafft_version(binary)
    jobs = [(binary, g, s) for g, s in sequences.items()]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        done = dict(pool.map(_align_one, jobs))
    return {g: done[g] for g in sequences}, version


def read_alignment(path: Path) -> dict[str, str]:
    return {k: v.upper() for k, v in read_fasta(path).items()}


# --- the step --------------------------------------------------------------------------------------
def main(build: Path, release: Path, contig: Path, out: Path, workers: int | None = None) -> None:
    """Writes into `out`: <gene>.nt.fasta (the tRNAs, aligned), gene_census.tsv, unresolved.tsv, alignments.tsv,
    rrna_check.tsv, rna_features.tsv (one row per gene and species)."""
    out.mkdir(parents=True, exist_ok=True)
    workers = workers or len(os.sched_getaffinity(0))
    recs = load_records(build, release, contig)
    print(f"{len(recs)} live species with a record in the release")
    per = {r["slug"]: extract_record(r["text"], r["phylum"]) for r in recs}
    cen = census(per)
    cen.to_csv(out / "gene_census.tsv", sep="\t", index=False)
    pd.DataFrame([{"slug": s, "aa": aa, "unresolved_copies": n} for s, r in per.items() for aa, n in sorted(r["unresolved"].items())],
                 columns=["slug", "aa", "unresolved_copies"]).to_csv(out / "unresolved.tsv", sep="\t", index=False)
    print(cen.to_string(index=False))

    tr = {g: {s: r["sequences"][g] for s, r in per.items() if g in r["sequences"]} for g in TRNA_GENES}
    aligned, version = align_genes(tr, workers)
    rows, feats = [], []
    for g in TRNA_GENES:
        text = aligned[g]
        (out / f"{file_name(g)}.nt.fasta").write_text(text)
        ali = b2align.parse_fasta(text) if text else []
        d = {k: s.upper() for k, _, s in ali}
        rows.append({"gene": g, "source": "aligned here", "aligner": f"MAFFT {version} " + " ".join(b2align.MAFFT_ARGS),
                     "sequences": len(d), "columns": len(ali[0][2]) if ali else 0,
                     "informative_columns": informative_columns(d), "sha256": sha256_text(text) if text else ""})
        feats += [{"gene": g, "slug": s, **gene_features(seq)} for s, seq in tr[g].items()]
    check = []
    for g in RRNA_GENES:
        path = build / "alignments" / f"{g}.nt.fasta.gz"
        d = read_alignment(path)
        mine = {s: r["sequences"][g] for s, r in per.items() if g in r["sequences"]}
        same = sum(1 for s in mine if s in d and d[s].replace("-", "") == mine[s])
        check.append({"gene": g, "in_record": len(mine), "in_build_alignment": len(d), "same_sequence": same,
                      "only_in_record": len(set(mine) - set(d)), "only_in_build": len(set(d) - set(mine)),
                      "different_sequence": len(set(mine) & set(d)) - same})
        raw = (build / "alignments" / f"{g}.nt.fasta.gz")
        rows.append({"gene": g, "source": "the build's alignment, not realigned", "aligner": "as the build (MAFFT v7.526, FFT-NS-2)",
                     "sequences": len(d), "columns": len(next(iter(d.values()))), "informative_columns": informative_columns(d),
                     "sha256": hashlib.sha256(raw.read_bytes()).hexdigest()})
        feats += [{"gene": g, "slug": s, **gene_features(seq)} for s, seq in d.items()]
    pd.DataFrame(check).to_csv(out / "rrna_check.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(out / "alignments.tsv", sep="\t", index=False)
    pd.DataFrame(feats).sort_values(["gene", "slug"], kind="mergesort").to_csv(out / "rna_features.tsv", sep="\t", index=False, float_format="%.6g")
    print(pd.DataFrame(rows)[["gene", "sequences", "columns", "informative_columns"]].to_string(index=False))
    print(pd.DataFrame(check).to_string(index=False))
    print(f"MAFFT {version}")
