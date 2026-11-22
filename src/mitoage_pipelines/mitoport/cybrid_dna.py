"""Nucleotide divergence between a cybrid's host and donor, beside H16's amino acid counts (exploratory).

Added on 2026-11-22, after H28's run had been read, on the maintainer's question whether DNA or amino acid divergence
matters more and whether the peptides read from the rRNA genes could add to the proteins. Nothing here can pass or
fail a rule: the scores go through H16's readouts T1 to T4 beside `d_all` and `f_all`, and the table says so.

Per (host, donor, complex): `nt_cds_d` and `nt_cds_f`, the differing columns (a gap is a difference) of the global
nucleotide alignments of the complex's mtDNA-encoded genes, as a count and as a share of the columns. Per (host,
donor), repeated on each of its complexes: the same for the two rRNA genes (`nt_rrna_d`, `nt_rrna_f`), for the 22
tRNA genes (`nt_trna_d`, `nt_trna_f`), and the count within the reading frames of humanin (in the 16S rRNA gene)
and of MOTS-c (in the 12S rRNA gene), whose places are those of the human record carried to each species by
alignment of the rRNA gene (`nt_humanin_d`, `nt_motsc_d`).

Added the same day, on the maintainer's question which regions matter most. Per complex, the differences of its
genes split by what they do: within codons that code the same amino acid in host and donor (`nt_syn_d`, silent)
and within codons that code another (`nt_nonsyn_d`), read on the alignment of the two proteins, vertebrate
mitochondrial code; `_f` is the count over three times the codons compared. Per pair: the control region as the
records annotate it (`nt_dloop_d`, `nt_dloop_f`; not assessed where a record annotates none), and, where the host's
record annotates them (the mouse's does, the human's does not), the differences within the host's conserved sequence
blocks 1 to 3 carried to the donor by alignment of the control regions (`nt_csb_d`). Promoters and binding sites
are annotated in no record here and are not scored. `dna_correlations.tsv` holds the rank correlations between the
scores across the donors.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

from . import cybrid

# the reading frames on the human record NC_012920.1, 1-based and inclusive, stop codon included
HUMANIN = (2633, 2707)      # in the 16S rRNA gene; 24 residues, MAPRGFSCLLLLTSEIDLPVKRRA
MOTSC = (1343, 1393)        # in the 12S rRNA gene; 16 residues, MRWQEMGYIFYPRKLR
RRNA = {"12s ribosomal rna": "12S", "s-rrna": "12S", "16s ribosomal rna": "16S", "l-rrna": "16S"}
PAIR_SCORES = ("nt_rrna_d", "nt_rrna_f", "nt_trna_d", "nt_trna_f", "nt_humanin_d", "nt_motsc_d",
               "nt_dloop_d", "nt_dloop_f", "nt_csb_d")
COMPLEX_SCORES = ("nt_cds_d", "nt_cds_f", "nt_syn_d", "nt_syn_f", "nt_nonsyn_d", "nt_nonsyn_f")
SCORES = (*COMPLEX_SCORES, *PAIR_SCORES)
CSB = ("CSB1", "CSB2", "CSB3")
BASELINE = ("d_all", "f_all")


def aligner() -> PairwiseAligner:
    al = PairwiseAligner()
    al.mode = "global"
    al.match_score, al.mismatch_score = 1, -1
    al.open_gap_score, al.extend_gap_score = -5, -1
    return al


def read_features(gb: Path) -> dict[str, dict]:
    """{accession.version: {"cds": {gene: nt}, "rrna": {"12S"|"16S": (nt, start)}, "trna": {(product, i): nt},
    "dloop": (nt, start) or None, "blocks": {"CSB1"..: (start, end)}}}. A coding sequence starts at its first whole
    codon; `start` is a 0-based place on the record's plus strand; tRNA genes of one product are numbered in the
    order of the record; the control region is its plus strand read in the genome's order (across the record's end
    where it is annotated in two parts), and `blocks` are places on the record. Refuses two rRNA genes of one kind
    or two control regions."""
    out = {}
    for rec in SeqIO.parse(str(gb), "genbank"):
        cds, rrna, trna, seen, dloop, blocks = {}, {}, {}, {}, None, {}
        for f in rec.features:
            product = f.qualifiers.get("product", [None])[0]
            if f.type == "CDS":
                g = cybrid.gene_name(f.qualifiers.get("gene", [None])[0], product)
                if g is not None:
                    cds[g] = str(f.extract(rec.seq)).upper()[int(f.qualifiers.get("codon_start", ["1"])[0]) - 1:]
            elif f.type == "rRNA":
                k = RRNA.get((product or "").strip().lower())
                if k is None:
                    raise ValueError(f"{rec.id}: an rRNA gene with the product {product!r}")
                if k in rrna:
                    raise ValueError(f"{rec.id}: two {k} rRNA genes")
                if f.location.strand != 1:
                    raise ValueError(f"{rec.id}: the {k} rRNA gene is not on the plus strand")
                rrna[k] = (str(f.extract(rec.seq)).upper(), int(f.location.start))
            elif f.type == "tRNA":
                i = seen.get(product, 0); seen[product] = i + 1
                trna[(product, i)] = str(f.extract(rec.seq)).upper()
            elif f.type == "D-loop" or (f.type == "misc_feature"
                                        and f.qualifiers.get("note", [""])[0].strip().lower() == "control region"):
                if dloop is not None:
                    raise ValueError(f"{rec.id}: two control regions")
                parts = sorted(f.location.parts, key=lambda q: int(q.start) == 0)      # the part at the record's
                dloop = ("".join(str(rec.seq[int(q.start):int(q.end)]) for q in parts).upper(),     # start comes last
                         int(parts[0].start))
            elif f.type == "misc_feature":
                name = f.qualifiers.get("note", [""])[0].split(";")[0].strip()
                if name in CSB:
                    blocks[name] = (int(f.location.start), int(f.location.end))
        out[rec.id] = {"cds": cds, "rrna": rrna, "trna": trna, "dloop": dloop, "blocks": blocks}
    return out


def columns(a: str, b: str, al: PairwiseAligner) -> list[tuple[str, str]]:
    """The columns of the best global alignment of two sequences, a gap as '-'."""
    if not a or not b:
        return [(x, "-") for x in a] + [("-", y) for y in b]
    aln = al.align(a, b)[0]
    return list(zip(str(aln[0]), str(aln[1])))


def differences(a: str, b: str, al: PairwiseAligner) -> tuple[int, int]:
    """(differing columns, columns) of the alignment; a gap is a difference."""
    cols = columns(a, b, al)
    return sum(x != y for x, y in cols), len(cols)


def carried(human: str, other: str, start: int, end: int, al: PairwiseAligner) -> str:
    """The part of `other` aligned to the places start..end (0-based, end excluded) of `human`, with what `other`
    has inserted between two of those places and without what it has inserted before the first or after the last."""
    out, i = [], 0
    for x, y in columns(human, other, al):
        if x != "-":
            i += 1
        inside = start < i < end if x == "-" else start < i <= end      # an insertion counts only between two places
        if inside and y != "-":
            out.append(y)
    return "".join(out)


def codons(cds: str) -> list[str]:
    """The whole codons of a coding sequence, a final stop codon left out."""
    out = [cds[i:i + 3] for i in range(0, len(cds) - len(cds) % 3, 3)]
    if out and str(Seq(out[-1]).translate(table=2)) == "*":
        out.pop()
    return out


def silent_changing(a: str, b: str, al: PairwiseAligner | None = None) -> tuple[int, int, int]:
    """(differing bases within codons that code the same amino acid, within codons that code another, codons
    compared) for two coding sequences, the codons paired by the alignment of the two proteins; a codon opposite a
    gap, or holding a base that is not A, C, G or T, is not compared."""
    ca, cb = codons(a), codons(b)
    pa, pb = (str(Seq("".join(c)).translate(table=2)) for c in (ca, cb))
    syn = non = n = 0
    for (a0, a1), (b0, b1) in zip(*(al or cybrid.aligner()).align(pa, pb)[0].aligned):
        for k in range(a1 - a0):
            x, y = ca[a0 + k], cb[b0 + k]
            if set(x + y) - set("ACGT"):
                continue
            n += 1
            d = sum(p != q for p, q in zip(x, y))
            if pa[a0 + k] == pb[b0 + k]:
                syn += d
            else:
                non += d
    return syn, non, n


def control_region(h: dict, d: dict, al: PairwiseAligner) -> dict:
    """The pair's control-region scores; NaN where a record annotates no control region, and for the blocks where
    the host's record annotates none."""
    out = {"nt_dloop_d": np.nan, "nt_dloop_f": np.nan, "nt_csb_d": np.nan}
    if h["dloop"] is None or d["dloop"] is None:
        return out
    (hs, start), (ds, _) = h["dloop"], d["dloop"]
    x, n = differences(hs, ds, al)
    out.update(nt_dloop_d=x, nt_dloop_f=x / n)
    if set(h["blocks"]) == set(CSB):
        total = 0
        for a, b in h["blocks"].values():
            if not 0 <= a - start < b - start <= len(hs):
                raise ValueError("a conserved sequence block lies outside the control region")
            total += differences(hs[a - start:b - start], carried(hs, ds, a - start, b - start, al), al)[0]
        out["nt_csb_d"] = total
    return out


def frames(features: dict, slug: str, human: dict, al: PairwiseAligner) -> dict[str, str]:
    """{"humanin": nt, "motsc": nt} of one species: the human record's reading frames carried to its rRNA genes."""
    mine = features[cybrid.GENOMES[slug]]["rrna"]
    out = {}
    for name, (first, last), kind in (("humanin", HUMANIN, "16S"), ("motsc", MOTSC, "12S")):
        seq, start = human["rrna"][kind]
        a, b = first - 1 - start, last - start
        if not 0 <= a < b <= len(seq):
            raise ValueError(f"the frame of {name} is not inside the human {kind} rRNA gene")
        out[name] = carried(seq, mine[kind][0], a, b, al)
    return out


def cds_of(features: dict, slug: str, gene: str) -> str:
    """The gene's nucleotides from the record H16 takes its protein from (the cell line's own where one is listed)."""
    acc = cybrid.OWN_LINE.get((slug, gene), cybrid.GENOMES[slug])
    return features[acc]["cds"][gene]


def pair_rows(features: dict, host: str, donor: str, al: PairwiseAligner | None = None) -> list[dict]:
    al = al or aligner()
    human = features[cybrid.GENOMES[cybrid.HUMAN]]
    h, d = features[cybrid.GENOMES[host]], features[cybrid.GENOMES[donor]]
    if set(h["rrna"]) != {"12S", "16S"} or set(d["rrna"]) != {"12S", "16S"}:
        raise ValueError(f"{host} / {donor}: an rRNA gene is missing")
    if set(h["trna"]) != set(d["trna"]) or len(h["trna"]) != 22:
        raise ValueError(f"{host} / {donor}: the tRNA genes do not pair ({len(h['trna'])} and {len(d['trna'])})")
    rr = [differences(h["rrna"][k][0], d["rrna"][k][0], al) for k in ("12S", "16S")]
    tr = [differences(h["trna"][k], d["trna"][k], al) for k in sorted(h["trna"])]
    fh, fd = frames(features, host, human, al), frames(features, donor, human, al)
    pair = {"nt_rrna_d": sum(x for x, _ in rr), "nt_rrna_f": sum(x for x, _ in rr) / sum(n for _, n in rr),
            "nt_trna_d": sum(x for x, _ in tr), "nt_trna_f": sum(x for x, _ in tr) / sum(n for _, n in tr),
            "nt_humanin_d": differences(fh["humanin"], fd["humanin"], al)[0],
            "nt_motsc_d": differences(fh["motsc"], fd["motsc"], al)[0], **control_region(h, d, al)}
    rows, pal = [], cybrid.aligner()
    for cx, genes in cybrid.COMPLEX_GENES.items():
        cd = [differences(cds_of(features, host, g), cds_of(features, donor, g), al) for g in genes]
        sc = [silent_changing(cds_of(features, host, g), cds_of(features, donor, g), pal) for g in genes]
        syn, non, n = (sum(x[i] for x in sc) for i in range(3))
        rows.append({"host": host, "donor": donor, "complex": cx, "nt_cds_d": sum(x for x, _ in cd),
                     "nt_cds_f": sum(x for x, _ in cd) / sum(n for _, n in cd), "nt_syn_d": syn,
                     "nt_syn_f": syn / (3 * n), "nt_nonsyn_d": non, "nt_nonsyn_f": non / (3 * n), **pair})
    return rows


def main(gb: Path, out: Path, measurements: Path = cybrid.MEASUREMENTS) -> pd.DataFrame:
    """Writes dna_scores.tsv, dna_calibration.tsv, dna_summary.tsv and dna_correlations.tsv into `out`, which must
    hold H16's scores.tsv (`mitoport cybrid`) for the amino acid baselines."""
    key = ["host", "donor", "complex"]
    if not (out / "scores.tsv").is_file():
        raise SystemExit(f"{out / 'scores.tsv'} is missing: run `mitoage mitoport cybrid` first")
    features = read_features(gb)
    meas = cybrid.read_measurements(measurements)
    al = aligner()
    rows = [r for h, d in sorted(meas[["host", "donor"]].drop_duplicates().itertuples(index=False))
            for r in pair_rows(features, h, d, al)]
    dna = pd.DataFrame(rows)
    dna.to_csv(out / "dna_scores.tsv", sep="\t", index=False, float_format="%.4f")
    scores = pd.read_csv(out / "scores.tsv", sep="\t")[key + list(BASELINE)].merge(dna, on=key, how="left")
    names = SCORES + BASELINE
    # T2b and T3 compare complexes or points within a complex: only the scores that differ between complexes enter;
    # a pair's RNA scores are the same on each of its complexes
    per_complex = (*COMPLEX_SCORES, *BASELINE)
    cal, summary = cybrid.summarise(scores, meas, names, t2b=per_complex, t3=per_complex,
                                    t4=("nt_cds_d", "nt_syn_d", "nt_nonsyn_d", "nt_rrna_d", "nt_trna_d",
                                        "nt_humanin_d", "nt_motsc_d", "nt_dloop_d", "d_all"))
    # how far the scores can be told apart at all: rank correlations across the donors, a complex score as its sum
    # over the complexes (hosts together; the two hosts' donors are few)
    per_pair = scores.groupby(["host", "donor"], sort=False).agg(
        {**{s: "sum" for s in ("nt_syn_d", "nt_nonsyn_d", "d_all")}, **{s: "first" for s in PAIR_SCORES if s.endswith("_d")}})
    corr = per_pair.corr(method="spearman", min_periods=5)
    corr.to_csv(out / "dna_correlations.tsv", sep="\t", float_format="%.3f")
    summary["reads"] = summary.reads + "; exploratory, added after H28's run"
    cal.to_csv(out / "dna_calibration.tsv", sep="\t", index=False, float_format="%.4f")
    summary.to_csv(out / "dna_summary.tsv", sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.width", 220, "display.max_colwidth", 130, "display.max_rows", 400):
        print(dna.drop_duplicates(["host", "donor"])[["host", "donor", *PAIR_SCORES]].to_string(index=False))
        print(summary[summary.test.isin(["T1", "T2"])].to_string(index=False))
        print(corr.round(2).to_string())
    print("wrote", out / "dna_scores.tsv", out / "dna_calibration.tsv", out / "dna_summary.tsv", out / "dna_correlations.tsv")
    return summary
