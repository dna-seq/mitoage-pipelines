"""RM16's first question for the catalogued peptides: is each published frame there in each species?

For a peptide of n codons placed in the human row (catalogue.place), its window is the human start codon's
alignment column and the 3n columns of the human frame. In a species the frame is read from the species' first
base at or after the start column, on the peptide's strand, codon by codon, until a stop. The answer is one of
three, never folded into another:
  kept           a start codon of the reading there, and no stop in the first ceil(0.9 n) codons
  lost           a stop in the first ceil(0.9 n) codons ("early stop"), or open without a start ("no start")
  not assessed   under half the window's columns hold a base ("gap"), a base other than ACGT before the decision
                 ("ambiguous base"), or the gene ends before ceil(0.9 n) codons with no stop ("runs off the gene")
"open" is the frame without the start condition (no early stop). Species: the build's rows that are not retired.
Each frame is read twice (codes.py): with the species' own code and with the standard code.
Writes states.tsv (one row per peptide, reading and species) and kept.tsv (per peptide, reading and class).
"""
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from ..seqio import read_fasta
from . import codes
from .catalogue import PEPTIDES, Peptide, place

OPEN_SHARE = 0.9        # a frame is open when no stop falls in its first ceil(0.9 n) codons
COVERAGE_MIN = 0.5      # under half the window's columns with a base: not assessed
SPECIES_MIN = 20        # a class is summarised (and tested) with at least this many assessed species
READINGS = ("own", "standard")
KEPT, LOST, NOT_ASSESSED = 1, 0, -1
BIG = 1 << 30


def needed(n: int) -> int:
    """Codons that must be open: ceil(0.9 n)."""
    return math.ceil(OPEN_SHARE * n - 1e-9)


@dataclass
class Host:
    """One alignment as arrays, on the peptide's strand (for "-" every row is reverse-complemented, gaps kept)."""
    slugs: list[str]
    nongap: np.ndarray          # (species, columns) bool
    cum: np.ndarray             # (species, columns + 1): bases before each column
    nt: np.ndarray              # (species, longest + 3) bases 0..3, 4 unknown or padding
    length: np.ndarray          # (species,) ungapped length

    def row(self, slug: str) -> int:
        return self.slugs.index(slug)

    def sequence(self, i: int) -> str:
        return "".join("ACGTN"[b] for b in self.nt[i, :self.length[i]])


def load_host(path: Path, strand: str, live: set[str]) -> Host:
    aln = read_fasta(path)
    slugs = sorted(s for s in aln if s in live)
    W = len(aln[slugs[0]])
    raw = np.frombuffer("".join(aln[s].upper() for s in slugs).encode(), dtype=np.uint8).reshape(len(slugs), W)
    nongap = raw != ord("-")
    base = codes.INDEX[raw]
    if strand == "-":
        nongap = nongap[:, ::-1]
        base = base[:, ::-1]
        base = np.where(base < 4, 3 - base, 4).astype(np.int8)
    length = nongap.sum(1)
    nt = np.full((len(slugs), int(length.max()) + 3), 4, dtype=np.int8)
    for i in range(len(slugs)):
        nt[i, :length[i]] = base[i][nongap[i]]
    cum = np.zeros((len(slugs), W + 1), dtype=np.int32)
    np.cumsum(nongap, axis=1, out=cum[:, 1:])
    return Host(slugs, np.ascontiguousarray(nongap), cum, nt, length)


@dataclass
class Reader:
    """Per species and position of the ungapped sequence: does a start codon begin here, and where is the next stop
    and the next unknown codon in the same frame (BIG when none)."""
    start: np.ndarray
    next_stop: np.ndarray
    next_bad: np.ndarray


def _next(flag: np.ndarray) -> np.ndarray:
    """For each position, the nearest position at or after it in the same frame where `flag` holds (BIG if none)."""
    S, L = flag.shape
    L3 = -(-L // 3) * 3
    f = np.zeros((S, L3), dtype=bool); f[:, :L] = flag
    pos = np.arange(L3).reshape(1, -1, 3).repeat(S, 0)
    arr = np.where(f.reshape(S, -1, 3), pos, BIG)
    nxt = np.minimum.accumulate(arr[:, ::-1, :], axis=1)[:, ::-1, :]
    return nxt.reshape(S, L3)[:, :L]


def reader(host: Host, tables: list[codes.Table]) -> Reader:
    """`tables[i]` is species i's reading."""
    S, L = host.nt.shape
    cidx = np.full((S, L), -2, dtype=np.int16)            # -2: no complete codon starts here (past the end)
    for i in range(S):
        cidx[i, :host.length[i]] = codes.codon_index(host.nt[i, :host.length[i]])
    stop = np.zeros((S, L), dtype=bool); start = np.zeros((S, L), dtype=bool)
    ok = cidx >= 0
    for t in {id(t): t for t in tables}.values():
        rows = np.array([tt is t for tt in tables])
        sub = cidx[rows]
        stop[rows] = np.where(sub >= 0, t.stop[np.clip(sub, 0, 63)], False)
        start[rows] = np.where(sub >= 0, t.start[np.clip(sub, 0, 63)], False)
    return Reader(start, _next(stop), _next(cidx == -1))


def state_at(host: Host, rd: Reader, c0: int, cols: np.ndarray, n: int) -> dict[str, np.ndarray]:
    """The frame of n codons whose human start is at column c0 and whose human bases are at `cols`, in every species.
    -> arrays over species: state, reason code, p (start position in the species' ungapped sequence), orf (codons
    before the first stop, capped by the gene's end), start (bool), open (bool), coverage."""
    S = len(host.slugs)
    m = needed(n)
    r = np.arange(S)
    cov = host.nongap[:, cols].mean(1)
    p = host.cum[:, c0]
    inside = p < host.length
    pp = np.where(inside, p, 0)
    k_stop = np.where(rd.next_stop[r, pp] >= BIG, BIG, (rd.next_stop[r, pp] - pp) // 3)
    k_bad = np.where(rd.next_bad[r, pp] >= BIG, BIG, (rd.next_bad[r, pp] - pp) // 3)
    avail = (host.length - pp) // 3
    start = rd.start[r, pp] & inside
    early = (k_stop < m) & (k_stop < k_bad)
    reason = np.full(S, "", dtype=object)
    state = np.full(S, NOT_ASSESSED, dtype=np.int8)
    gap = (cov < COVERAGE_MIN) | ~inside
    amb = ~gap & ~early & (k_bad < m)
    off = ~gap & ~early & ~amb & (avail < m)
    decided = ~gap & ~amb & ~off
    reason[gap] = "gap"; reason[amb] = "ambiguous base"; reason[off] = "runs off the gene"
    state[decided & ~early & start] = KEPT
    state[decided & (early | ~start)] = LOST
    reason[decided & ~early & start] = "kept"
    reason[decided & early & start] = "early stop"
    reason[decided & early & ~start & (k_bad > 0)] = "no start, early stop"
    reason[decided & early & ~start & (k_bad == 0)] = "early stop"        # start codon unknown, stop decides
    reason[decided & ~early & ~start] = "no start"
    orf = np.minimum(k_stop, avail)
    return {"state": state, "reason": reason, "p": p, "orf": orf, "start": start, "open": decided & ~early,
            "coverage": cov}


@dataclass
class Placed:
    """A catalogued peptide on its host alignment: the human window in the host's (strand-oriented) columns."""
    pep: Peptide
    host: Host
    q0: int                 # the human start's index in the strand-oriented human sequence
    hcols: np.ndarray       # the human row's columns, strand-oriented
    n: int

    @property
    def c0(self) -> int:
        return int(self.hcols[self.q0])

    def cols(self, q: int | None = None) -> np.ndarray:
        q = self.q0 if q is None else q
        return self.hcols[q:q + 3 * self.n]


def species_table(data: Path) -> pl.DataFrame:
    return pl.read_parquet(data / "species.parquet").filter(pl.col("record_status") != "retired").select(
        "slug", "class", "order", "family", "genus", "transl_table")


def placed(data: Path, pep: Peptide, live: set[str], hosts: dict | None = None) -> Placed:
    key = (pep.alignment, pep.strand)
    if hosts is not None and key in hosts:
        host = hosts[key]
    else:
        host = load_host(data / "alignments" / f"{pep.alignment}.fasta.gz", pep.strand, live)
        if hosts is not None:
            hosts[key] = host
    pl_ = place(pep, data)
    h = host.row("homo-sapiens")
    L = int(host.length[h])
    q0 = pl_["human_index"] if pep.strand == "+" else L - 1 - pl_["human_index"]
    return Placed(pep, host, q0, np.flatnonzero(host.nongap[h]), len(pep.sequence))


def tables_for(pc: Placed, reading: str, code_of: dict[str, int]) -> list[codes.Table]:
    name = codes.species_reading(reading, pc.pep.reading)
    return [codes.tables(code_of[s], name) for s in pc.host.slugs]


def observe(pc: Placed, reading: str, code_of: dict[str, int]) -> tuple[dict, list[codes.Table]]:
    tabs = tables_for(pc, reading, code_of)
    rd = reader(pc.host, tabs)
    return state_at(pc.host, rd, pc.c0, pc.cols(), pc.n), tabs


def identity(pc: Placed, i: int, p: int, table: codes.Table) -> tuple[float, str]:
    """Residues 2..n of the species' first n codons against the published peptide's, position by position."""
    nt = "".join("ACGTN"[b] for b in pc.host.nt[i, p:p + 3 * pc.n])
    aa = codes.translate(nt, table)
    if len(aa) < 2:
        return float("nan"), aa
    same = sum(a == b for a, b in zip(aa[1:], pc.pep.sequence[1:pc.n]))
    return same / (pc.n - 1), aa


def main(data: Path, out: Path) -> tuple[pl.DataFrame, pl.DataFrame]:
    out.mkdir(parents=True, exist_ok=True)
    sp = species_table(data)
    meta = {r["slug"]: r for r in sp.iter_rows(named=True)}
    code_of = {s: int(r["transl_table"]) for s, r in meta.items()}
    live = set(meta)
    hosts: dict = {}
    rows = []
    for pep in PEPTIDES:
        pc = placed(data, pep, live, hosts)
        for reading in READINGS:
            st, tabs = observe(pc, reading, code_of)
            for i, slug in enumerate(pc.host.slugs):
                ident, aa = (identity(pc, i, int(st["p"][i]), tabs[i]) if st["state"][i] != NOT_ASSESSED
                             else (float("nan"), ""))
                start = "".join("ACGTN"[b] for b in pc.host.nt[i, st["p"][i]:st["p"][i] + 3]) if st["p"][i] < pc.host.length[i] else ""
                m = meta[slug]
                rows.append({"peptide": pep.name, "reading": reading, "slug": slug, "class": m["class"],
                             "order": m["order"], "family": m["family"], "transl_table": code_of[slug],
                             "state": {KEPT: "kept", LOST: "lost", NOT_ASSESSED: "not assessed"}[int(st["state"][i])],
                             "reason": st["reason"][i], "open": bool(st["open"][i]),
                             "start": bool(st["start"][i]) and st["state"][i] != NOT_ASSESSED, "start_codon": start,
                             "orf_codons": int(st["orf"][i]) if st["state"][i] != NOT_ASSESSED else None,
                             "coverage": round(float(st["coverage"][i]), 4),
                             "identity": None if ident != ident else round(ident, 4), "translation": aa})
        # species with no row in the host alignment: not assessed, said so
        for slug in sorted(live - set(pc.host.slugs)):
            m = meta[slug]
            for reading in READINGS:
                rows.append({"peptide": pep.name, "reading": reading, "slug": slug, "class": m["class"],
                             "order": m["order"], "family": m["family"], "transl_table": code_of[slug],
                             "state": "not assessed", "reason": "no host row", "open": False, "start": False, "start_codon": "",
                             "orf_codons": None, "coverage": 0.0, "identity": None, "translation": ""})
    states = pl.DataFrame(rows, infer_schema_length=None).sort(["peptide", "reading", "slug"])
    states.write_csv(out / "states.tsv", separator="\t")
    kept = summarise(states)
    kept.write_csv(out / "kept.tsv", separator="\t")
    return states, kept


def summarise(states: pl.DataFrame) -> pl.DataFrame:
    """Per peptide, reading and class (and "all"): counts of the three answers and the shares among the assessed."""
    allc = states.with_columns(pl.lit("all").alias("class"))
    both = pl.concat([states, allc])
    g = both.group_by(["peptide", "reading", "class"]).agg(
        pl.len().alias("species"),
        (pl.col("state") == "kept").sum().alias("kept"),
        (pl.col("state") == "lost").sum().alias("lost"),
        (pl.col("state") == "not assessed").sum().alias("not_assessed"),
        (pl.col("open") & (pl.col("state") != "not assessed")).sum().alias("open"),
        (pl.col("start") & (pl.col("state") != "not assessed")).sum().alias("start"),
        (pl.col("reason") != "no host row").sum().alias("with_region"),
        ((pl.col("reason") == "early stop") | (pl.col("reason") == "no start, early stop")).sum().alias("early_stop"),
        (pl.col("reason") == "no start").sum().alias("no_start_only"),
        pl.col("identity").filter(pl.col("state") == "kept").median().alias("median_identity_kept"),
    )
    g = g.with_columns((pl.col("kept") + pl.col("lost")).alias("assessed"))
    g = g.with_columns(
        (pl.col("kept") / pl.col("assessed")).alias("kept_share"),
        (pl.col("open") / pl.col("assessed")).alias("open_share"),
        pl.when(pl.col("assessed") >= SPECIES_MIN).then(pl.lit("summarised"))
          .otherwise(pl.lit(f"not assessed: under {SPECIES_MIN} assessed species")).alias("status"))
    return g.select("peptide", "reading", "class", "species", "with_region", "assessed", "open", "start", "kept",
                    "lost", "not_assessed", "early_stop", "no_start_only", "kept_share", "open_share", "median_identity_kept", "status"
                    ).sort(["peptide", "reading", "class"])
