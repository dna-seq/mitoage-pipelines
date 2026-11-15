"""Alignments of the mitochondrial genes across the live species of the build, as the site's
`bun run align` makes them (site repository: pipeline/align.ts, pipeline/steps/alignments.ts,
pipeline/lib/align.ts; decision D5 in docs/DECISIONS.md).

- Each of the 13 proteins as RefSeq translates it (/translation), the first copy of a gene in the record,
  aligned with MAFFT v7.526 FFT-NS-2 (`--retree 2 --maxiterate 0 --quiet --thread N`), species sorted by
  slug; the codon alignment the protein alignment implies (each residue -> its codon, a gap -> three, the
  stop codon left out); the two rRNAs as nucleotides, aligned the same way.
- Records with only a CONTIG line take their sequence from efetch of the same version.
- Output: <GENE>.aa / .codon / .nt .fasta.gz (one line per sequence) and alignments.json with the SHA-256
  of each uncompressed FASTA, the absent and duplicated genes and the codon rows that do not translate.
"""
import gzip
import hashlib
import json
import re
import subprocess
import tempfile
from pathlib import Path

import polars as pl

from .. import external
from . import genome as gm
from . import pins

MAFFT_ARGS = ["--retree", "2", "--maxiterate", "0", "--quiet"]
RRNA_GENES = ("rRNA-12S", "rRNA-16S")
ALIGNED = gm.PROTEIN_GENES + RRNA_GENES


def fasta(seqs: list[tuple[str, str, str]]) -> str:
    return "".join(f">{i}{' ' + d if d else ''}\n{s}\n" for i, d, s in seqs)


def parse_fasta(text: str) -> list[tuple[str, str, str]]:
    out, cur, parts = [], None, []
    for raw in re.split(r"\r?\n", text):
        line = raw.strip()
        if not line: continue
        if line.startswith(">"):
            if cur: out.append((*cur, "".join(parts)))
            head = line[1:].strip(); sp = re.search(r"\s", head)
            cur = (head, "") if not sp else (head[:sp.start()], head[sp.start() + 1:].strip()); parts = []
        else:
            if cur is None: raise ValueError("FASTA: sequence before the first header")
            parts.append(line)
    if cur: out.append((*cur, "".join(parts)))
    return out


def mafft_version(binary: str) -> str:
    p = subprocess.run([binary, "--version"], capture_output=True, text=True)
    m = re.search(r"v\d+\.\d+", p.stderr + p.stdout)
    if not m: raise RuntimeError(f"could not read MAFFT's version from: {(p.stderr + p.stdout).strip()}")
    return m.group(0)


def align(binary: str, seqs: list[tuple[str, str, str]], threads: int) -> list[tuple[str, str, str]]:
    """MAFFT on the sequences; the result must hold exactly them, in their order, apart from gaps."""
    with tempfile.TemporaryDirectory(prefix="mitoage-align-") as d:
        f = Path(d) / "in.fasta"; f.write_text(fasta(seqs))
        p = subprocess.run([binary, *MAFFT_ARGS, "--thread", str(threads), str(f)], capture_output=True, text=True)
        if p.returncode != 0: raise RuntimeError(f"MAFFT exited with {p.returncode}: {p.stderr.strip()[-400:]}")
    aligned = [(i, d, s.upper()) for i, d, s in parse_fasta(p.stdout)]
    if len(aligned) != len(seqs): raise RuntimeError(f"alignment has {len(aligned)} rows for {len(seqs)} sequences")
    cols = len(aligned[0][2])
    for (ai, _, a), (si, _, s) in zip(aligned, seqs):
        if ai != si or len(a) != cols or a.replace("-", "") != s.upper():
            raise RuntimeError(f"alignment row {ai} is not its input sequence {si} with gaps")
    return aligned


def collect(species: list[dict], records: dict[str, str]) -> dict:
    """Per gene: the protein (RefSeq's translation) and coding sequence, or the rRNA, of every species."""
    out = {"protein": {g: [] for g in gm.PROTEIN_GENES}, "cds": {g: {} for g in gm.PROTEIN_GENES},
           "rrna": {g: [] for g in RRNA_GENES}, "absent": [], "duplicates": []}
    for s in sorted(species, key=lambda x: x["slug"]):
        rec = gm.parse_genbank(records[s["refseq"]])
        if not rec.sequence: raise ValueError(f"record {s['refseq']} of {s['slug']} carries no sequence")
        gm.apply_gene_patches(rec, None)
        desc = f"{s['scientific_name']} {s['refseq']}"
        copies: dict[str, int] = {}
        for f in rec.features:
            if f.key not in ("CDS", "rRNA"): continue
            n = gm.feature_name(f)
            if not n: continue
            name = n[0]; seen = copies.get(name, 0); copies[name] = seen + 1
            if seen > 0: continue
            if f.key == "CDS" and name in out["protein"]:
                tr = f.qualifiers.get("translation", [None])[0]
                if not tr:
                    copies[name] = seen; continue
                out["protein"][name].append((s["slug"], desc, tr.upper()))
                start = int(f.qualifiers.get("codon_start", ["1"])[0])
                out["cds"][name][s["slug"]] = gm.extract(rec.sequence, f.segments)[start - 1:].upper()
            elif f.key == "rRNA" and name in out["rrna"]:
                out["rrna"][name].append((s["slug"], desc, gm.extract(rec.sequence, f.segments).upper()))
        for g in ALIGNED:
            c = copies.get(g, 0)
            if c == 0: out["absent"].append({"gene": g, "slug": s["slug"], "reason": "not annotated in the record"})
            elif c > 1: out["duplicates"].append({"gene": g, "slug": s["slug"], "copies": c})
    return out


def thread_codons(aligned_protein: str, cds: str) -> str:
    """The codon row a protein row implies: a residue -> its codon, a gap -> three gaps; a stop codon is left out."""
    residues = len(aligned_protein.replace("-", ""))
    if len(cds) < residues * 3: raise ValueError(f"coding sequence of {len(cds)} nt is too short for {residues} residues")
    if len(cds) > residues * 3 + 3: raise ValueError(f"coding sequence of {len(cds)} nt is too long for {residues} residues")
    out, k = [], 0
    for ch in aligned_protein:
        if ch == "-": out.append("---")
        else: out.append(cds[k:k + 3]); k += 3
    return "".join(out)


def codon_mismatches(aligned_protein: str, aligned_codons: str, table: int) -> int:
    """Codons that do not translate to the residue above them (first codon and ambiguity codes not counted)."""
    code = gm.genetic_code(table)
    bad, first = 0, True
    for i, aa in enumerate(aligned_protein):
        if aa == "-": continue
        if first: first = False; continue
        expect = code.get(aligned_codons[i * 3:i * 3 + 3].replace("T", "U"))
        if expect and aa != "X" and expect != aa: bad += 1
    return bad


def summarise(aligned: list[tuple[str, str, str]]) -> dict:
    n, cols = len(aligned), len(aligned[0][2])
    occupied = invariant = conserved = 0
    for i in range(cols):
        counts: dict[str, int] = {}; gaps = 0
        for _, _, s in aligned:
            ch = s[i].upper()
            if ch == "-": gaps += 1
            else: counts[ch] = counts.get(ch, 0) + 1
        if gaps * 2 > n: continue
        occupied += 1
        top = max(counts.values()) / (n - gaps)
        if top == 1: invariant += 1
        if top >= 0.9: conserved += 1
    return {"columns": cols, "occupiedColumns": occupied, "invariantColumns": invariant, "conserved90Columns": conserved}


def main(species_parquet: Path, release: Path, release_label: str, contig: Path, out: Path, threads: int) -> dict:
    from .download import md5_of
    from .refseq import scan
    binary = external.find("mafft", env="MAFFT")
    v = mafft_version(binary)
    if v != pins.MAFFT_VERSION:
        raise external.MissingProgram(f"MAFFT {v} found at {binary}; the alignments are pinned to {pins.MAFFT_VERSION}")
    if not species_parquet.is_file():
        raise FileNotFoundError(f"no {species_parquet}: run `mitoage build2 assemble` first")
    sp = pl.read_parquet(species_parquet).filter(pl.col("record_status") != "retired")
    species = sp.select("slug", "scientific_name", "refseq", "transl_table").to_dicts()
    want = {s["refseq"] for s in species}
    records = {h.version: t for h, t in scan(release, keep=lambda h: h.version in want) if t is not None}
    for ver, t in list(records.items()):
        if "\nORIGIN" not in t:
            p = contig / f"{ver}.gb"
            if not p.is_file(): raise FileNotFoundError(f"{ver} has only a CONTIG line: `mitoage download contig`")
            records[ver] = p.read_text()
    missing = want - set(records)
    if missing: raise ValueError(f"not in the release: {', '.join(sorted(missing)[:5])}")
    seqs = collect(species, records)
    table = {s["slug"]: s["transl_table"] for s in species}
    out.mkdir(parents=True, exist_ok=True)
    files, notes = [], []

    def write(gene, kind, suffix, aligned):
        text = fasta(aligned)
        name = f"{gene}.{suffix}.fasta.gz"
        gz = gzip.compress(text.encode(), compresslevel=9, mtime=0)
        (out / name).write_bytes(gz)
        files.append({"file": name, "gene": gene, "kind": kind, "sequences": len(aligned), **summarise(aligned),
                      "sha256": hashlib.sha256(text.encode()).hexdigest(), "bytes": len(gz)})
        print(f"  {name}: {len(aligned)} sequences, {files[-1]['columns']} columns", flush=True)

    for gene in gm.PROTEIN_GENES:
        aligned = align(binary, seqs["protein"][gene], threads)
        write(gene, "protein", "aa", aligned)
        codons = []
        for i, d, s in aligned:
            try:
                c = thread_codons(s, seqs["cds"][gene][i])
                bad = codon_mismatches(s, c, table[i])
                if bad > 0: notes.append({"gene": gene, "slug": i, "reason": f"{bad} codons do not translate to the aligned residue (RNA editing, a frameshift or an exception in the record)"})
                codons.append((i, d, c))
            except ValueError as e:
                notes.append({"gene": gene, "slug": i, "reason": f"left out of the codon alignment: {e}"})
        write(gene, "codon", "codon", codons)
    for gene in RRNA_GENES:
        write(gene, "rRNA", "nt", align(binary, seqs["rrna"][gene], threads))
    manifest = {"aligner": {"program": "MAFFT", "version": v, "arguments": " ".join(MAFFT_ARGS), "method": "FFT-NS-2 (progressive)"},
                "source": {"name": "NCBI RefSeq mitochondrion release", "version": release_label, "md5": md5_of(release)},
                "species": len(species), "files": files, "absent": seqs["absent"], "duplicates": seqs["duplicates"], "codonNotes": notes}
    (out / "alignments.json").write_text(json.dumps(manifest, indent=1) + "\n")
    return manifest
