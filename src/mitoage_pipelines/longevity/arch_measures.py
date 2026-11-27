"""H38, step 1: the architecture of every live species' mitochondrial genome, as numbers.

Per species (one row of `architecture.tsv`), read from the RefSeq record exactly as the audit of H13 reads it
(`audit.records.genome_of`, the rules of build 2's `genome.annotate`, the pinned release 237):

(a) Rearrangement. The record's signed circular arrangement (`audit.order.arrangement`, canonical under rotation and
    under reversal, so a record deposited from the other strand is not rearranged) against the most common arrangement
    of its class, at the two tiers of the audit: `full` (13 protein genes, 2 rRNAs and the tRNAs; unknown for a record
    with no tRNA annotated, never 0) and `genes` (the 13 and the 2 rRNAs). The distance is the breakpoint distance
    (`breakpoint_distance`): the adjacencies of the circular order, each read the same from either strand, as a multiset,
    and half the size of the difference of the two multisets. For two records with the same genes it is the number of
    adjacencies of one that the other lacks, the standard breakpoint distance of signed circular genomes; for records with
    different genes (a gene absent, a gene twice) the adjacencies of the absent or extra gene count on their side, so the
    distance is still symmetric and still zero only for equal arrangements. A record with no arrangement at the tier is
    unknown, as is a record whose class is unknown; the class's most common arrangement counts every live record of
    the build that has a class and an arrangement, the ties going to the arrangement first in string order.
(b) The control region: the number of control-region features the record carries (`cr_copies`) and their total length
    (`cr_length`). A record with no control-region feature is unknown for both, since an unannotated control region
    cannot be told from an absent one; `cr_annotated` keeps the plain count, zero included.
(c) Tandem repeats: Tandem Repeats Finder at version 4.09.1 (the bioconda package, build h031d066_4; the binary's
    SHA-256 is pinned in `TRF_SHA256` and any other binary is refused) with Benson's standard parameters `2 7 7 80 10 50
    500` (its redundancy elimination left on); the share of the genome, and of the control region, covered by
    a repeat (the union of the repeats' intervals). The raw output of every record is kept in `trf_repeats.tsv`, and its
    SHA-256 is written to `architecture_inputs.tsv`. TRF is the standard finder; none is written here. The circular
    genome is read as a line, so a repeat across the origin is seen as two.
(d) The non-coding content the build stores (`non_coding_bp`, and its share of `mtdna_size_bp`), and the total length of
    gene overlaps: the positions covered by two or more features among the CDS, rRNA and tRNA ones, a feature annotated
    twice at the same coordinates on the same strand counting once.

Nothing here reads a lifespan.
"""
import hashlib
import os
import subprocess
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl

from ..audit import order as aorder
from ..audit import records as rec
from .. import external

TRF_VERSION = "4.09.1"
TRF_BUILD = "bioconda trf 4.09.1 build h031d066_4 (the program reports Version 4.09)"
TRF_SHA256 = "dc8e610afe6674f770639fa0bcf026e09ba1b51588de68b8710c3d3466859fa8"
TRF_PARAMS = (2, 7, 7, 80, 10, 50, 500)          # match, mismatch, indel, match probability, indel probability, minimum score, maximum period
TRF_CHUNK = 40                                    # records per TRF process
TIERS = ("full", "genes")
CODING = ("CDS", "rRNA", "tRNA")
REPEAT_COLUMNS = ["slug", "start", "end", "period", "copies", "consensus_size", "percent_match", "percent_indel", "score",
                  "consensus"]


# --- (a) the distance ------------------------------------------------------------------------------
def _flip(tok: str) -> str:
    return ("-" if tok[0] == "+" else "+") + tok[1:]


def adjacencies(order: str) -> Counter:
    """The adjacencies of a circular signed order ('+COX1 -tRNA-Ser(UCN) ...'), as a multiset: each neighbouring pair
    (a, b), with the last next to the first, read the same from the other strand, where it is (-b, -a)."""
    toks = order.split(" ") if order else []
    n = len(toks)
    out: Counter = Counter()
    for i in range(n):
        a, b = toks[i], toks[(i + 1) % n]
        out[min((a, b), (_flip(b), _flip(a)))] += 1
    return out


def breakpoint_distance(a: str | None, b: str | None) -> float | None:
    """Half the size of the difference of the two orders' adjacency multisets; None when either order is unknown."""
    if a is None or b is None or (isinstance(a, float) and np.isnan(a)) or (isinstance(b, float) and np.isnan(b)):
        return None
    ca, cb = adjacencies(a), adjacencies(b)
    return sum(abs(ca[k] - cb[k]) for k in set(ca) | set(cb)) / 2


def modal_arrangement(orders: pd.Series) -> str | None:
    """The most common value of the arrangements; the first in string order on a tie; None when there are none."""
    known = orders.dropna()
    if known.empty:
        return None
    c = known.value_counts()
    return min(c.index[c.values == c.max()])


def class_distances(tab: pd.DataFrame, tier: str, group: str = "class") -> pd.DataFrame:
    """For each record of `tab` (columns `group` and order_<tier>_unoriented), the group's most common arrangement and the
    breakpoint distance to it. Unknown (missing) when the record has no arrangement or no group."""
    col = f"order_{tier}_unoriented"
    modal = {}
    for g, d in tab[tab[group].notna()].groupby(group, sort=True):
        modal[g] = modal_arrangement(d[col])
    ref = tab[group].map(modal)
    dist = [breakpoint_distance(a, b) if isinstance(b, str) else None for a, b in zip(tab[col], ref)]
    return pd.DataFrame({f"modal_{tier}": ref, f"rearr_{tier}": pd.array(dist, dtype="Float64")}, index=tab.index)


# --- (b) and (d), from one annotated genome ---------------------------------------------------------
def control_region(g) -> tuple[int, int | None]:
    """(control-region features, their total length in bases); the length is None when there is no such feature."""
    cr = [f for f in g.features if f.type == "D-loop"]
    if not cr:
        return 0, None
    return len(cr), sum(s.end - s.start + 1 for f in cr for s in f.segments)


def gene_overlap(g) -> int:
    """Positions covered by two or more CDS, rRNA or tRNA features; a feature repeated with the same segments and strand
    counts once."""
    n = len(g.record.sequence)
    depth = np.zeros(n + 2, dtype=np.int32)
    seen = set()
    for f in g.features:
        if f.type not in CODING:
            continue
        key = (f.type, f.name, f.strand, tuple((s.start, s.end) for s in f.segments))
        if key in seen:
            continue
        seen.add(key)
        for s in f.segments:
            depth[max(s.start, 1):min(s.end, n) + 1] += 1
    return int((depth[1:n + 1] >= 2).sum())


def cr_mask(g) -> np.ndarray:
    n = len(g.record.sequence)
    mask = np.zeros(n + 1, dtype=bool)
    for f in g.features:
        if f.type == "D-loop":
            for s in f.segments:
                mask[max(s.start, 1):min(s.end, n) + 1] = True
    return mask[1:]


# --- (c) tandem repeats ----------------------------------------------------------------------------
def find_trf() -> str:
    """The pinned TRF binary: on PATH, with the SHA-256 of build h031d066_4 of the bioconda package 4.09.1."""
    binary = external.find("trf")
    h = hashlib.sha256(Path(binary).resolve().read_bytes()).hexdigest()
    if h != TRF_SHA256:
        raise SystemExit(f"{binary} is not the pinned Tandem Repeats Finder ({TRF_BUILD}): its SHA-256 is {h}, not {TRF_SHA256}")
    return binary


def run_trf(binary: str, seqs: dict[str, str], workdir: Path | None = None) -> dict[str, list[dict]]:
    """Tandem Repeats Finder over `seqs` (name -> bases) in one process; name -> its repeats, in the order TRF reports them.
    TRF exits with status 1 after a normal run, so the output file and its parameter line are checked instead; a sequence
    TRF did not report on is an error, never an empty answer."""
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        tmp = Path(tmp)
        fa = tmp / "in.fa"
        fa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        p = subprocess.run([binary, "in.fa", *map(str, TRF_PARAMS), "-d", "-h"], cwd=tmp, capture_output=True, text=True)
        dat = tmp / ("in.fa." + ".".join(map(str, TRF_PARAMS)) + ".dat")
        if not dat.is_file():
            raise RuntimeError(f"TRF wrote no output (exit {p.returncode}): {p.stdout[-300:]} {p.stderr[-300:]}")
        return parse_dat(dat.read_text(), list(seqs))


def parse_dat(text: str, names: list[str]) -> dict[str, list[dict]]:
    """The .dat file of TRF -d: `Sequence: <name>`, a `Parameters:` line, then one line per repeat (start end period copies
    consensus-size percent-matches percent-indels score A C G T entropy consensus repeat). Every name must have its block
    and the parameters must be those asked for."""
    out: dict[str, list[dict]] = {}
    cur = None
    want = " ".join(map(str, TRF_PARAMS))
    for line in text.splitlines():
        if line.startswith("Sequence: "):
            cur = line[len("Sequence: "):].split()[0]
            out[cur] = []
        elif line.startswith("Parameters: "):
            if line[len("Parameters: "):].strip() != want:
                raise RuntimeError(f"TRF ran with parameters {line!r}, not {want}")
        elif cur is not None and line.strip() and line[0].isdigit():
            f = line.split()
            out[cur].append({"start": int(f[0]), "end": int(f[1]), "period": int(f[2]), "copies": float(f[3]),
                             "consensus_size": int(f[4]), "percent_match": int(f[5]), "percent_indel": int(f[6]),
                             "score": int(f[7]), "consensus": f[13]})
    missing = [n for n in names if n not in out]
    if missing:
        raise RuntimeError(f"TRF gave no block for {len(missing)} sequences, first {missing[0]}")
    return out


def covered(repeats: list[dict], n: int) -> np.ndarray:
    """A boolean mask over positions 1..n: inside at least one repeat."""
    mask = np.zeros(n + 1, dtype=bool)
    for r in repeats:
        mask[max(r["start"], 1):min(r["end"], n) + 1] = True
    return mask[1:]


def repeat_shares(repeats: list[dict], n: int, cr: np.ndarray | None) -> tuple[float, float | None]:
    """(share of the genome covered by a repeat, share of the control region covered); the second is None without a
    control region."""
    m = covered(repeats, n)
    whole = float(m.sum()) / n if n else float("nan")
    if cr is None or not cr.any():
        return whole, None
    return whole, float((m & cr).sum()) / float(cr.sum())


# --- one record -------------------------------------------------------------------------------------
def measure_record(row: dict) -> dict:
    """Everything one record gives, apart from the class distances (they need the whole table) and the repeats (run in
    chunks): the audit's arrangement row, the control region, the gene overlaps, the bases."""
    g = rec.genome_of(pd.Series(row))
    a = aorder.row_of({"slug": row["slug"]}, g)
    copies, length = control_region(g)
    mask = cr_mask(g)
    seq = g.record.sequence.upper()
    return {"slug": row["slug"], "status": a["status"], "genes_present": a["genes_present"], "trnas_annotated": a["trnas_annotated"],
            "order_genes_unoriented": a["order_genes_unoriented"], "order_full_unoriented": a["order_full_unoriented"],
            "cr_annotated": copies, "cr_copies": copies if copies else None, "cr_length": length,
            "cr_bases": int(mask.sum()), "gene_overlap_bp": gene_overlap(g), "genome_length": len(seq), "_seq": seq,
            "_cr": mask}


def _measure_chunk(args) -> list[dict]:
    rows, binary, workdir = args
    out = [measure_record(r) for r in rows]
    reps = run_trf(binary, {m["slug"]: m["_seq"] for m in out}, workdir)
    for m in out:
        r = reps[m["slug"]]
        m["_repeats"] = r
        m["repeat_genome"], m["repeat_cr"] = repeat_shares(r, m["genome_length"], m["_cr"])
        del m["_seq"], m["_cr"]
    return out


def measure(records: pd.DataFrame, workers: int, workdir: Path | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The per-record table and the repeats table, from `audit.records.load`'s table; records not in the release are left
    out here (they are not live species with a record, and H37 refuses them as well)."""
    binary = find_trf()
    live = records[records.in_release & records.genbank.notna()].sort_values("slug", kind="mergesort")
    rows = live.to_dict("records")
    chunks = [(rows[i:i + TRF_CHUNK], binary, workdir) for i in range(0, len(rows), TRF_CHUNK)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        done = [m for part in pool.map(_measure_chunk, chunks) for m in part]
    rep = pd.DataFrame([{"slug": m["slug"], **r} for m in done for r in m["_repeats"]], columns=REPEAT_COLUMNS)
    tab = pd.DataFrame([{k: v for k, v in m.items() if not k.startswith("_")} for m in done])
    return tab, rep


def table(records: pd.DataFrame, species: pd.DataFrame, workers: int, workdir: Path | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The architecture table of the live species and the repeats table. `species`: the build's species table."""
    sp = species[["slug", "class", "order", "family", "non_coding_bp", "mtdna_size_bp"]][species.record_status.ne("retired")].copy()
    tab, rep = measure(records[records.slug.isin(set(sp.slug))], workers, workdir)
    tab = sp.merge(tab, on="slug", how="inner", validate="one_to_one").reset_index(drop=True)
    for tier in TIERS:
        tab = tab.join(class_distances(tab, tier))
    tab["non_coding_frac"] = tab.non_coding_bp / tab.mtdna_size_bp
    tab["non_coding_bp"] = tab.non_coding_bp.astype(float)
    return tab.sort_values("slug", kind="mergesort").reset_index(drop=True), rep.sort_values(["slug", "start", "end"], kind="mergesort").reset_index(drop=True)


def sha256_of(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(data: Path, interim: Path, out: Path, workers: int | None = None) -> pd.DataFrame:
    """Reads the audit's records.parquet (`audit records`) and the build's species table; writes architecture.tsv,
    trf_repeats.tsv and architecture_inputs.tsv into `out`."""
    workers = workers or len(os.sched_getaffinity(0))
    records = rec.load(interim)
    species = pl.read_parquet(data / "species.parquet").to_pandas()
    out.mkdir(parents=True, exist_ok=True)
    tab, reps = table(records, species, workers)
    tab.to_csv(out / "architecture.tsv", sep="\t", index=False, float_format="%.6g", na_rep="")
    reps.to_csv(out / "trf_repeats.tsv", sep="\t", index=False, na_rep="")
    pd.DataFrame([
        {"item": "species.parquet", "sha256": sha256_of(data / "species.parquet")},
        {"item": "audit records.parquet", "sha256": sha256_of(interim / "records.parquet")},
        {"item": "trf binary", "sha256": TRF_SHA256},
        {"item": f"trf {TRF_VERSION}", "sha256": "parameters " + " ".join(map(str, TRF_PARAMS)) + "; " + TRF_BUILD},
        {"item": "trf_repeats.tsv", "sha256": sha256_of(out / "trf_repeats.tsv")},
        {"item": "architecture.tsv", "sha256": sha256_of(out / "architecture.tsv")},
    ]).to_csv(out / "architecture_inputs.tsv", sep="\t", index=False)
    n = len(tab)
    print(f"{n} live species; repeats {len(reps)}; unknown at the full tier {int(tab.rearr_full.isna().sum())}, at the genes tier "
          f"{int(tab.rearr_genes.isna().sum())}; no control region {int(tab.cr_length.isna().sum())}")
    return tab
