"""Possible duplicates and possible mislabels, as candidates for a person to check, never verdicts.

Near-identical pairs. Nucleotide identity of two records over the 13 protein genes, read from the build's codon
alignments concatenated (stop codons are not in them): matching bases over the columns where both have A, C, G or
T. A pair is reported when the two records are of different species (every record of the build is), share at least
MIN_NT such columns and have identity >= NEAR (0.99); at >= IDENTICAL (0.999) it is also marked near-identical.
Why 0.99, chosen before any identity was computed: mitochondrial divergence within an animal species is usually
under 1 %, while sister species usually differ by more than 2 % at COX1 (the "barcode gap", Hebert et al. 2003, Proc
R Soc B 270:313-321). Two species' records within 1 % over all 13 genes are where one species would be: the same
animal under two names, a misidentified specimen, a hybrid's or an introgressed mitogenome, or a very young pair of
species. Pairs of records whose whole sequences are equal are reported whatever the alignment says.

Out-of-place records. Amino-acid identity over the 13 proteins, from the build's protein alignments concatenated,
over columns where both have one of the 20 amino acids, at least MIN_AA of them. A record's nearest neighbour is the
record of highest identity (the first slug on a tie). It is out of its order when some record outside its order is
strictly closer than every other record of its order; likewise for its class. A record whose order (class) has no
other record compared with it is unknown, not "in place". Proteins, not COX1 alone: more sites, and less saturated
than nucleotides at the depths of orders and classes. The build's class scheme has paraphyletic classes (Reptilia
without birds, Actinopterygii beside Sarcopterygii), which make some out-of-class candidates expected.

Identities are computed as integer counts by matrix products of 0/1 matrices in float32: every partial sum is an
integer below 2^24, so the counts are exact whatever the order of summation, and a rerun gives the same bytes.
"""
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from ..build2.genome import PROTEIN_GENES
from ..seqio import read_fasta
from . import records as rec
from .common import NO, UNKNOWN, YES, Warnings, write_tsv

NEAR = 0.99
IDENTICAL = 0.999
MIN_NT = 5000
MIN_AA = 1500
NT = "ACGT"
AA = "ACDEFGHIKLMNPQRSTVWY"
LEVELS = ("genus", "family", "order", "class")


def concatenate(alignments: dict[str, dict[str, str]], slugs: list[str]) -> list[str]:
    """Each species' rows of every gene, in gene order, gaps where a gene has no row for it."""
    out = {s: [] for s in slugs}
    for g in sorted(alignments):
        aln = alignments[g]
        width = len(next(iter(aln.values())))
        for s in slugs:
            out[s].append(aln.get(s, "-" * width).upper())
    return ["".join(out[s]) for s in slugs]


def identity_counts(seqs: list[str], alphabet: str) -> tuple[np.ndarray, np.ndarray]:
    """(matches, compared): for every pair, the columns where both have the same letter of the alphabet, and the
    columns where both have a letter of the alphabet."""
    if len({len(s) for s in seqs}) > 1:
        raise ValueError("rows of unequal length")
    n = len(seqs)
    x = np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8).reshape(n, -1)
    match = np.zeros((n, n), dtype=np.float32)
    valid = np.zeros(x.shape, dtype=np.float32)
    for ch in alphabet:
        o = (x == ord(ch)).astype(np.float32)
        match += o @ o.T
        valid += o
    compared = valid @ valid.T
    if max(x.shape[1], 1) >= 2 ** 24:
        raise ValueError("too many columns for exact float32 counts")
    return np.rint(match).astype(np.int64), np.rint(compared).astype(np.int64)


def level(a: dict, b: dict) -> str:
    """The rank at which two records' taxonomy first agrees, read upward from genus."""
    for lv in LEVELS:
        if a[lv] is not None and not pd.isna(a[lv]) and a[lv] == b[lv]:
            return "same genus" if lv == "genus" else f"same {lv}, different {LEVELS[LEVELS.index(lv) - 1]}"
    return "different class"


def pairs(tax: pd.DataFrame, m: np.ndarray, k: np.ndarray, seq_hash: dict[str, str]) -> pd.DataFrame:
    slugs = list(tax.index)
    with np.errstate(divide="ignore", invalid="ignore"):
        ident = np.where(k >= MIN_NT, m / np.maximum(k, 1), np.nan)
    iu = np.triu_indices(len(slugs), 1)
    hit = set(zip(*[a[(ident[iu] >= NEAR)] for a in iu]))
    pos = {s: i for i, s in enumerate(slugs)}
    by_hash: dict[str, list[str]] = {}
    for s, h in sorted(seq_hash.items()):
        by_hash.setdefault(h, []).append(s)
    same_seq = set()
    for group in by_hash.values():
        for i, a in enumerate(group):
            for b in group[i + 1:]:
                same_seq.add(tuple(sorted((a, b))))
    keys = {tuple(sorted((slugs[i], slugs[j]))) for i, j in hit} | same_seq
    rows = []
    for a, b in keys:
        ta, tb = tax.loc[a].to_dict(), tax.loc[b].to_dict()
        i, j = pos.get(a), pos.get(b)
        known = i is not None and j is not None and k[i, j] >= MIN_NT
        rows.append({"slug_a": a, "slug_b": b, "refseq_a": ta["refseq"], "refseq_b": tb["refseq"],
                     "species_a": ta["scientific_name"], "species_b": tb["scientific_name"], "relation": level(ta, tb),
                     "identity": round(float(ident[i, j]), 6) if known else None,
                     "compared_nt": int(k[i, j]) if i is not None and j is not None else None,
                     "differences": int(k[i, j] - m[i, j]) if known else None,
                     "near_identical": (YES if ident[i, j] >= IDENTICAL else NO) if known else UNKNOWN,
                     "same_sequence": YES if (a, b) in same_seq else NO,
                     **{f"{lv}_a": ta[lv] for lv in ("class", "order", "family", "genus")},
                     **{f"{lv}_b": tb[lv] for lv in ("class", "order", "family", "genus")}})
    cols = ["slug_a", "slug_b", "refseq_a", "refseq_b", "species_a", "species_b", "relation", "identity", "compared_nt",
            "differences", "near_identical", "same_sequence", "class_a", "order_a", "family_a", "genus_a", "class_b", "order_b",
            "family_b", "genus_b"]
    return pd.DataFrame(rows, columns=cols)


def nearest(tax: pd.DataFrame, m: np.ndarray, k: np.ndarray) -> pd.DataFrame:
    slugs = list(tax.index)
    n = len(slugs)
    with np.errstate(divide="ignore", invalid="ignore"):
        ident = np.where(k >= MIN_AA, m / np.maximum(k, 1), -1.0)
    np.fill_diagonal(ident, -1.0)
    rows = []
    for i, s in enumerate(slugs):
        t = tax.iloc[i]
        row = {"slug": s, "refseq": t.refseq, "class": t["class"], "order": t["order"], "family": t["family"], "genus": t["genus"],
               "compared_records": int((ident[i] >= 0).sum())}
        if row["compared_records"] == 0:
            rows.append({**row, "out_of_order": UNKNOWN, "out_of_class": UNKNOWN, "reason": f"no record shares {MIN_AA} residues with it"})
            continue
        j = int(np.argmax(ident[i]))
        nb = tax.iloc[j]
        row.update({"nn_slug": slugs[j], "nn_identity": round(float(ident[i, j]), 6), "nn_class": nb["class"], "nn_order": nb["order"],
                    "nn_family": nb["family"]})
        reasons = []
        for lv in ("order", "class"):
            own = t[lv]
            if own is None or pd.isna(own):
                row[f"out_of_{lv}"] = UNKNOWN; reasons.append(f"no {lv} in the build"); continue
            same = (tax[lv].to_numpy() == own) & (ident[i] >= 0)
            other = (tax[lv].to_numpy() != own) & tax[lv].notna().to_numpy() & (ident[i] >= 0)
            if not same.any():
                row[f"out_of_{lv}"] = UNKNOWN; reasons.append(f"no other record of its {lv} compared"); continue
            js = int(np.flatnonzero(same)[np.argmax(ident[i][same])])
            row[f"best_in_{lv}_slug"] = slugs[js]; row[f"best_in_{lv}_identity"] = round(float(ident[i, js]), 6)
            if other.any():
                jo = int(np.flatnonzero(other)[np.argmax(ident[i][other])])
                row[f"margin_{lv}"] = round(float(ident[i, jo] - ident[i, js]), 6)
                row[f"out_of_{lv}"] = YES if ident[i, jo] > ident[i, js] else NO
            else:
                row[f"out_of_{lv}"] = NO
        rows.append({**row, "reason": "; ".join(reasons)})
    cols = ["slug", "refseq", "class", "order", "family", "genus", "compared_records", "nn_slug", "nn_identity", "nn_class", "nn_order",
            "nn_family", "best_in_order_slug", "best_in_order_identity", "margin_order", "out_of_order", "best_in_class_slug",
            "best_in_class_identity", "margin_class", "out_of_class", "reason"]
    return pd.DataFrame(rows).reindex(columns=cols)


def load_alignments(build: Path, kind: str) -> dict[str, dict[str, str]]:
    out = {}
    for g in PROTEIN_GENES:
        p = build / "alignments" / f"{g}.{kind}.fasta.gz"
        if not p.is_file():
            raise FileNotFoundError(f"no {p}: run `mitoage download build` first")
        out[g] = read_fasta(p)
    return out


def main(interim: Path, build: Path, out: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = rec.load(interim)
    warn = Warnings("audit duplicates")
    nt = load_alignments(build, "codon")
    aa = load_alignments(build, "aa")
    in_aln = sorted(set().union(*[set(a) for a in nt.values()]) | set().union(*[set(a) for a in aa.values()]))
    tax = df.set_index("slug")[["refseq", "scientific_name", "class", "order", "family", "genus"]]
    unknown = sorted(set(in_aln) - set(tax.index))
    if unknown:
        raise ValueError(f"alignment rows of species not in the records: {', '.join(unknown[:5])}")
    for s in sorted(set(tax.index) - set(in_aln)):
        warn.add("in no alignment (not compared)", s)
    t = tax.loc[in_aln]
    seq_hash = {r.slug: hashlib.sha256(rec.genome_of(r).record.sequence.encode()).hexdigest()
                for r in df.itertuples(index=False) if r.in_release}
    m, k = identity_counts(concatenate(nt, in_aln), NT)
    universe = sorted(set(in_aln) | set(seq_hash))
    pr = pairs(tax.loc[universe], *_pad(m, k, in_aln, universe), seq_hash)
    ma, ka = identity_counts(concatenate(aa, in_aln), AA)
    nn = nearest(t, ma, ka)
    for r in nn.itertuples(index=False):
        if r.out_of_order == UNKNOWN or r.out_of_class == UNKNOWN:
            warn.add(r.reason or "unknown", r.slug)
    nn = pd.concat([nn, pd.DataFrame([{"slug": s, "refseq": tax.loc[s, "refseq"], "class": tax.loc[s, "class"], "order": tax.loc[s, "order"],
                                        "family": tax.loc[s, "family"], "genus": tax.loc[s, "genus"], "compared_records": 0,
                                        "out_of_order": UNKNOWN, "out_of_class": UNKNOWN, "reason": "in no alignment of the build"}
                                       for s in sorted(set(tax.index) - set(in_aln))], columns=nn.columns)], ignore_index=True)
    warn.emit()
    write_tsv(pr, out / "pairs.tsv", ["identity", "slug_a", "slug_b"], [False, True, True])
    write_tsv(nn, out / "nearest.tsv", ["slug"])
    print(f"{len(in_aln)} records compared; pairs at identity >= {NEAR}: {int((pr.identity >= NEAR).sum())} "
          f"({pr.relation.value_counts().sort_index().to_dict()}); same sequence: {int((pr.same_sequence == YES).sum())}; "
          f"out of order: {int((nn.out_of_order == YES).sum())}, out of class: {int((nn.out_of_class == YES).sum())}")
    return pr, nn


def _pad(m: np.ndarray, k: np.ndarray, rows: list[str], want: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """The count matrices over `want`: species in no alignment get zero compared columns."""
    idx = {s: i for i, s in enumerate(rows)}
    sel = np.array([idx.get(s, -1) for s in want])
    mm = np.zeros((len(want), len(want)), dtype=np.int64); kk = np.zeros_like(mm)
    ok = sel >= 0
    mm[np.ix_(ok, ok)] = m[np.ix_(sel[ok], sel[ok])]; kk[np.ix_(ok, ok)] = k[np.ix_(sel[ok], sel[ok])]
    return mm, kk
