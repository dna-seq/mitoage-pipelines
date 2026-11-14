"""H7: substitution rates of the mitochondrial proteins against lifespan, by the relative-rate test.

Why not branch rates. The published disagreement in this literature is about normalisation: whether
lifespan predicts the substitution rate depends on how the rate is divided by divergence, generation
time or the tree (Hernandez-Pacheco et al. 2026, GBE, doi:10.1093/gbe/evag067). A cherry of the dated
tree - two tips sharing a parent - removes the question: both species have had exactly the same time
since they split. Measuring each of the two against one common outgroup turns their rate difference
into a contrast that needs no divergence-time correction at all.

For every cherry (A, B) with the nearest outgroup O:

    d_dS  = dS(A,O) - dS(B,O)        and likewise for dN and for the dN/dS ratio
    d_ll  = log10 lifespan(A) - log10 lifespan(B)
    d_lm  = log10 body mass(A) - log10 body mass(B)

and d_dS is regressed on d_ll through the origin, as contrasts are. Four normalisations are reported,
because the choice is what moves results in the published record: the raw difference, the difference
per million years of the cherry's depth, the dN/dS ratio, and the raw difference with d_lm in the
model. A sign test on the same pairs asks the question directly: does the longer-lived sister carry
the lower rate?

dN and dS are Nei-Gojobori counts with the Jukes-Cantor correction, in each species' own
mitochondrial genetic code (table 2 where the three species disagree). Codon pairs differing at more
than one position are averaged over every mutational pathway that avoids a stop. A pair whose dS
exceeds 1 is flagged `saturated`, and every test is reported with and without those pairs.

Species: the default dataset of H1, so this test and the amino-acid tests see the same animals.
"""
import gzip
import itertools
import math
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Data import CodonTable
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ..genes import MT_GENES
from ..seqio import read_fasta
from .aa_pgls import default_dataset

CONCAT = "CONCAT"
NORMALISATIONS = ("dS_raw", "dS_per_myr", "dN_raw", "dNdS", "dS_mass_adj")
SATURATED_DS = 1.0
MIN_CODONS = 30
MIN_CONTRASTS = 12


# --- Nei-Gojobori -------------------------------------------------------------------------------
_CODE: dict[int, tuple[dict[str, str], dict[str, float]]] = {}


def code_tables(table: int) -> tuple[dict[str, str], dict[str, float]]:
    """Codon -> amino acid, and codon -> number of synonymous sites, for one genetic code."""
    if table in _CODE:
        return _CODE[table]
    tab = CodonTable.unambiguous_dna_by_id[table]
    aa = dict(tab.forward_table)
    for stop in tab.stop_codons:
        aa[stop] = "*"
    syn = {}
    for codon, a in aa.items():
        s = 0.0
        for i in range(3):
            for b in "ACGT":
                if b == codon[i]:
                    continue
                if aa.get(codon[:i] + b + codon[i + 1:]) == a and a != "*":
                    s += 1 / 3
        syn[codon] = s
    _CODE[table] = (aa, syn)
    return _CODE[table]


def counted_codons(s1: str, s2: str, table: int):
    """The aligned codon pairs both sequences resolve: no gap, no ambiguity, no stop."""
    aa, _ = code_tables(table)
    for i in range(0, min(len(s1), len(s2)) - 2, 3):
        c1, c2 = s1[i:i + 3].upper(), s2[i:i + 3].upper()
        if c1 not in aa or c2 not in aa or aa[c1] == "*" or aa[c2] == "*":
            continue
        yield c1, c2


def sites(s1: str, s2: str, table: int) -> tuple[float, float, int]:
    """Synonymous and nonsynonymous sites, averaged over the two sequences, and the codons counted."""
    _, syn = code_tables(table)
    S = 0.0
    n = 0
    for c1, c2 in counted_codons(s1, s2, table):
        S += (syn[c1] + syn[c2]) / 2
        n += 1
    return S, 3 * n - S, n


def differences(s1: str, s2: str, table: int) -> tuple[float, float]:
    """Synonymous and nonsynonymous differences, averaging codon pairs over their mutational paths."""
    aa, _ = code_tables(table)
    Sd = Nd = 0.0
    for c1, c2 in counted_codons(s1, s2, table):
        if c1 == c2:
            continue
        at = [k for k in range(3) if c1[k] != c2[k]]
        if len(at) == 1:
            if aa[c1] == aa[c2]:
                Sd += 1
            else:
                Nd += 1
            continue
        paths = 0
        sd = nd = 0.0
        for order in itertools.permutations(at):
            cur, s_, n_, ok = c1, 0.0, 0.0, True
            for k in order:
                nxt = cur[:k] + c2[k] + cur[k + 1:]
                if aa.get(nxt) == "*":
                    ok = False
                    break
                if aa[nxt] == aa[cur]:
                    s_ += 1
                else:
                    n_ += 1
                cur = nxt
            if ok:
                paths += 1
                sd += s_
                nd += n_
        if paths:
            Sd += sd / paths
            Nd += nd / paths
    return Sd, Nd


def jukes_cantor(p: float) -> float:
    """The Jukes-Cantor correction; nan when the proportion is at or past saturation."""
    if not np.isfinite(p) or p < 0:
        return np.nan
    if p == 0:
        return 0.0
    if p >= 0.75:
        return np.nan
    return -0.75 * math.log(1 - 4 * p / 3)


def distance(s1: str, s2: str, table: int) -> dict:
    """dN and dS between two aligned coding sequences; nan when too few codons resolve."""
    S, N, codons = sites(s1, s2, table)
    if codons < MIN_CODONS or S <= 0 or N <= 0:
        return {"dN": np.nan, "dS": np.nan, "codons": codons}
    Sd, Nd = differences(s1, s2, table)
    return {"dN": jukes_cantor(Nd / N), "dS": jukes_cantor(Sd / S), "codons": codons}


# --- cherries -----------------------------------------------------------------------------------
def cherries(tree, usable: set[str]):
    """Every pair of tips that share a parent and are both usable, with that parent."""
    for clade in tree.get_nonterminals():
        kids = clade.clades
        if len(kids) == 2 and all(k.is_terminal() for k in kids):
            a, b = kids[0].name, kids[1].name
            if a in usable and b in usable:
                yield clade, a, b


def parent_map(tree) -> dict:
    return {child: parent for parent in tree.find_clades(order="level") for child in parent.clades}


def outgroup_for(tree, parents: dict, node, exclude: set[str], usable: set[str]) -> str | None:
    """The nearest usable tip outside the cherry, walking up one ancestor at a time."""
    current = node
    while current in parents:
        parent = parents[current]
        candidates = [t for t in parent.get_terminals() if t.name not in exclude and t.name in usable]
        if candidates:
            return min(candidates, key=lambda t: tree.distance(node, t)).name
        current = parent
    return None


# --- the contrasts and the tests ------------------------------------------------------------------
def build_contrasts(data: Path, species: pd.DataFrame, tree_file: Path) -> pd.DataFrame:
    """One row per cherry and gene: dN and dS of each sister against the outgroup, and the contrasts."""
    alignments = {g: read_fasta(data / "alignments" / f"{g}.codon.fasta.gz") for g in MT_GENES}
    aligned = set.intersection(*[set(a) for a in alignments.values()])
    usable = aligned & set(species.slug)
    tree = Phylo.read(tree_file, "newick")
    parents = parent_map(tree)
    info = species.set_index("slug")
    rows = []
    for node, a, b in cherries(tree, usable):
        o = outgroup_for(tree, parents, node, {a, b}, usable)
        if o is None:
            continue
        tables = {int(info.loc[s, "transl_table"]) for s in (a, b, o) if not pd.isna(info.loc[s, "transl_table"])}
        table = tables.pop() if len(tables) == 1 else 2
        if table not in CodonTable.unambiguous_dna_by_id:
            table = 2
        time_my = float(np.mean([k.branch_length or 0.0 for k in node.clades]))
        for gene in [*MT_GENES, CONCAT]:
            if gene == CONCAT:
                sa, sb, so = ("".join(alignments[g][s] for g in MT_GENES) for s in (a, b, o))
            else:
                sa, sb, so = alignments[gene][a], alignments[gene][b], alignments[gene][o]
            ra, rb = distance(sa, so, table), distance(sb, so, table)
            rows.append({"class": info.loc[a, "class"], "gene": gene, "sp_a": a, "sp_b": b, "outgroup": o,
                         "transl_table": table, "time_my": time_my, "codons": ra["codons"],
                         "dS_a": ra["dS"], "dS_b": rb["dS"], "dN_a": ra["dN"], "dN_b": rb["dN"],
                         "d_ll": info.loc[a, "ll"] - info.loc[b, "ll"],
                         "d_lm": info.loc[a, "lm"] - info.loc[b, "lm"]})
    c = pd.DataFrame(rows)
    if c.empty:
        return c
    c["d_dS"] = c.dS_a - c.dS_b
    c["d_dN"] = c.dN_a - c.dN_b
    c["w_a"] = c.dN_a / c.dS_a.replace(0, np.nan)
    c["w_b"] = c.dN_b / c.dS_b.replace(0, np.nan)
    c["d_w"] = c.w_a - c.w_b
    c["saturated"] = c[["dS_a", "dS_b"]].max(axis=1) > SATURATED_DS
    return c


def test_contrasts(df: pd.DataFrame, y: str, per_time: bool = False, covariates: list[str] | None = None) -> dict | None:
    """Regression of a rate contrast on the lifespan contrast through the origin, plus a sign test."""
    cols = [y, "d_ll", *(covariates or []), *(["time_my"] if per_time else [])]
    d = df.dropna(subset=cols).copy()
    if per_time:
        d = d[d.time_my > 0]
        d[y] = d[y] / d.time_my
    d = d[np.isfinite(d[y])]
    if len(d) < MIN_CONTRASTS:
        return None
    X = d[["d_ll", *(covariates or [])]].to_numpy()
    yv = d[y].to_numpy()
    beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
    resid = yv - X @ beta
    dof = len(yv) - X.shape[1]
    se = np.sqrt(np.diag(np.linalg.pinv(X.T @ X) * (resid @ resid / dof)))
    # A fit with no residual leaves the slope with no standard error: that is a degenerate test, not a
    # certain one, so t and p are reported as unknown rather than as infinity and zero.
    t = beta[0] / se[0] if se[0] > 0 else np.nan
    signed = d[d.d_ll != 0]
    lower = int(((signed.d_ll > 0) & (signed[y] < 0)).sum() + ((signed.d_ll < 0) & (signed[y] > 0)).sum())
    return {"n": len(d), "slope": beta[0], "se": se[0], "t": t,
            "p": 2 * stats.t.sf(abs(t), dof) if np.isfinite(t) else np.nan,
            "lower_rate_in_longer_lived": lower, "sign_n": len(signed),
            "sign_p": stats.binomtest(lower, len(signed), 0.5).pvalue if len(signed) else np.nan}


def select(data: Path, composition: Path, fourfold: Path, tree_file: Path) -> pd.DataFrame:
    """The default dataset of H1, across classes, with the genetic code of each species."""
    sp = pd.read_csv(composition, sep="\t").merge(
        pd.read_parquet(data / "species.parquet",
                        columns=["slug", "record_status", "lifespan_quality", "lifespan_record_type", "transl_table"]),
        on="slug", how="left")
    sp = sp[default_dataset(sp)].dropna(subset=["lm"])
    sp = sp.merge(pd.read_csv(fourfold, sep="\t"), on="slug")
    present = {t.name for t in Phylo.read(tree_file, "newick").get_terminals()}
    return sp[sp.slug.isin(present)].reset_index(drop=True)


def main(data: Path, composition: Path, fourfold: Path, tree_file: Path, out: Path) -> pd.DataFrame:
    """Writes rate_contrasts.tsv (one row per cherry and gene) and rates.tsv (the tests)."""
    out.mkdir(parents=True, exist_ok=True)
    sp = select(data, composition, fourfold, tree_file)
    print(f"default dataset with a tip in the dated tree: {len(sp)} species")
    contrasts = build_contrasts(data, sp, tree_file)
    if contrasts.empty:
        raise SystemExit("no cherry of the dated tree has two species of the default dataset")
    contrasts.to_csv(out / "rate_contrasts.tsv", sep="\t", index=False, float_format="%.6g")
    concat = contrasts[contrasts.gene == CONCAT]
    print(f"cherries: {concat.shape[0]} pairs; rows {len(contrasts)}; "
          f"dS > {SATURATED_DS} in {concat.saturated.mean():.1%} of concatenated pairs")

    rows = []
    for cls in [*sorted(contrasts["class"].dropna().unique()), "all"]:
        by_class = contrasts if cls == "all" else contrasts[contrasts["class"] == cls]
        for with_saturated in (False, True):
            sub = by_class if with_saturated else by_class[~by_class.saturated]
            for gene in [*MT_GENES, CONCAT]:
                g = sub[sub.gene == gene]
                for label in NORMALISATIONS:
                    kw = {"dS_raw": {"y": "d_dS"},
                          "dS_per_myr": {"y": "d_dS", "per_time": True},
                          "dN_raw": {"y": "d_dN"},
                          "dNdS": {"y": "d_w"},
                          "dS_mass_adj": {"y": "d_dS", "covariates": ["d_lm"]}}[label]
                    r = test_contrasts(g, **kw)
                    if r:
                        rows.append({"class": cls, "gene": gene, "normalisation": label,
                                     "includes_saturated": with_saturated, **r})
    res = pd.DataFrame(rows)
    res["q"] = res.groupby(["class", "normalisation", "includes_saturated"]).p.transform(
        lambda p: pd.Series(multipletests(p, method="fdr_bh")[1], index=p.index))
    res = res.sort_values("p").reset_index(drop=True)
    res.to_csv(out / "rates.tsv", sep="\t", index=False, float_format="%.5g")
    head = res[(~res.includes_saturated) & (res.gene == CONCAT)].sort_values("p")
    with pd.option_context("display.width", 220):
        print("\nconcatenated genes, saturated pairs excluded:")
        print(head.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(f"\n{len(res)} tests; q < 0.05 in {(res.q < 0.05).sum()}")
    return res
