"""H37 part A: the per-branch rate of the RNA genes, by H7's relative-rate test, and the positive control that
must pass before any of them is read.

H7 measures the rate difference of two sister tips of the dated tree (a cherry) against one common outgroup, so
that no divergence time is needed, and regresses it on the lifespan contrast through the origin (`rates.py`). Its
distances are Nei-Gojobori dN and dS, which a tRNA or an rRNA does not have. The rate here is the same contrast
with the distance H7 gives a stretch with no codons: the Jukes-Cantor distance over the alignment columns that
all three species resolve (A, C, G or T), the same shared-columns rule H7 uses for its codons, with a floor of
MIN_SITES columns where H7's is MIN_CODONS codons. Cherries, outgroup, contrast, regression, saturation flag
(distance above 1), the three normalisations that apply (the raw difference, per million years of the cherry's
depth, with the body-mass contrast in the model; dN and dN/dS do not exist for an RNA) and the sign test are
H7's own functions.

Positive control (`control`). The same driver, in codon mode, on the 13 protein genes and their concatenation, with
the cherries of H7 (species with all 13 proteins), writes a table that must equal H7's `rates.tsv` for the same
build row by row and column by column at the precision the table is written. If it does not, the pipeline is wrong
and `require_control` refuses every RNA rate and composition test.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Data import CodonTable
from statsmodels.stats.multitest import multipletests

from ..genes import MT_GENES
from ..seqio import read_fasta
from . import rates, species_sets

CONCAT = rates.CONCAT
MIN_SITES = 30                                   # the floor of shared columns for an RNA distance (H7's floor is 30 codons)
NT_NORMALISATIONS = {"d_raw": {"y": "d_dist"},
                     "d_per_myr": {"y": "d_dist", "per_time": True},
                     "d_mass_adj": {"y": "d_dist", "covariates": ["d_lm"]}}
CODON_NORMALISATIONS = {"dS_raw": {"y": "d_dS"}, "dS_per_myr": {"y": "d_dS", "per_time": True}, "dN_raw": {"y": "d_dN"},
                        "dNdS": {"y": "d_w"}, "dS_mass_adj": {"y": "d_dS", "covariates": ["d_lm"]}}
CONTROL_FILE = "h7_control.json"
BASES = "ACGT"


# --- the distance ----------------------------------------------------------------------------------
def nt_distance(s1: str, s2: str) -> dict:
    """Jukes-Cantor distance of two aligned rows over the columns both resolve (A, C, G, T in both); nan below the floor."""
    pairs = [(a, b) for a, b in zip(s1.upper(), s2.upper()) if a in BASES and b in BASES]
    if len(pairs) < MIN_SITES:
        return {"dist": np.nan, "sites": len(pairs)}
    p = sum(a != b for a, b in pairs) / len(pairs)
    return {"dist": rates.jukes_cantor(p), "sites": len(pairs)}


def nt_outgroup_distances(s1: str, s2: str, outgroup: str) -> tuple[dict, dict]:
    """Both sisters' distances to one outgroup on the same columns, resolved in all three rows (H7's rule for codons)."""
    if len({len(s1), len(s2), len(outgroup)}) != 1:
        raise ValueError("the three aligned rows must have equal lengths")
    rows = (s1.upper(), s2.upper(), outgroup.upper())
    keep = [i for i in range(len(s1)) if all(r[i] in BASES for r in rows)]
    a, b, o = ("".join(r[i] for i in keep) for r in rows)
    return nt_distance(a, o), nt_distance(b, o)


# --- the contrasts ---------------------------------------------------------------------------------
def _table_of(info: pd.DataFrame, names: tuple[str, ...]) -> int:
    tables = {int(info.loc[s, "transl_table"]) for s in names if not pd.isna(info.loc[s, "transl_table"])}
    table = tables.pop() if len(tables) == 1 else 2
    return table if table in CodonTable.unambiguous_dna_by_id else 2


def codon_contrasts(species: pd.DataFrame, tree_file: Path, alignments: dict[str, dict[str, str]]) -> pd.DataFrame:
    """H7's contrasts, for the genes of `alignments` (codon alignments) and their concatenation, in the order of the
    genes given: one row per cherry and gene. Cherries are those of the species all genes have."""
    genes = list(alignments)
    usable = set.intersection(*[set(a) for a in alignments.values()]) & set(species.slug)
    tree = Phylo.read(tree_file, "newick")
    parents = rates.parent_map(tree)
    info = species.set_index("slug")
    rows = []
    for node, a, b in rates.cherries(tree, usable):
        o = rates.outgroup_for(tree, parents, node, {a, b}, usable)
        if o is None:
            continue
        table = _table_of(info, (a, b, o))
        time_my = float(np.mean([k.branch_length or 0.0 for k in node.clades]))
        for gene in [*genes, CONCAT]:
            if gene == CONCAT:
                sa, sb, so = ("".join(alignments[g][s] for g in genes) for s in (a, b, o))
            else:
                sa, sb, so = alignments[gene][a], alignments[gene][b], alignments[gene][o]
            ra, rb = rates.outgroup_distances(sa, sb, so, table)
            rows.append({"class": info.loc[a, "class"], "gene": gene, "sp_a": a, "sp_b": b, "outgroup": o,
                         "transl_table": table, "time_my": time_my, "codons": ra["codons"],
                         "dS_a": ra["dS"], "dS_b": rb["dS"], "dN_a": ra["dN"], "dN_b": rb["dN"],
                         "d_ll": info.loc[a, "ll"] - info.loc[b, "ll"], "d_lm": info.loc[a, "lm"] - info.loc[b, "lm"]})
    c = pd.DataFrame(rows)
    if c.empty:
        return c
    c["d_dS"] = c.dS_a - c.dS_b
    c["d_dN"] = c.dN_a - c.dN_b
    c["w_a"] = c.dN_a / c.dS_a.replace(0, np.nan)
    c["w_b"] = c.dN_b / c.dS_b.replace(0, np.nan)
    c["d_w"] = c.w_a - c.w_b
    c["saturated"] = c[["dS_a", "dS_b"]].max(axis=1) > rates.SATURATED_DS
    return c


def nt_contrasts(species: pd.DataFrame, tree_file: Path, alignments: dict[str, dict[str, str]]) -> pd.DataFrame:
    """The same contrasts for nucleotide alignments: one row per cherry and gene. A cherry needs the gene in both
    sisters and in an outgroup, so each gene has the cherries of its own species."""
    tree = Phylo.read(tree_file, "newick")
    parents = rates.parent_map(tree)
    info = species.set_index("slug")
    rows = []
    for gene, ali in alignments.items():
        usable = set(ali) & set(species.slug)
        for node, a, b in rates.cherries(tree, usable):
            o = rates.outgroup_for(tree, parents, node, {a, b}, usable)
            if o is None:
                continue
            ra, rb = nt_outgroup_distances(ali[a], ali[b], ali[o])
            rows.append({"class": info.loc[a, "class"], "gene": gene, "sp_a": a, "sp_b": b, "outgroup": o,
                         "time_my": float(np.mean([k.branch_length or 0.0 for k in node.clades])), "sites": ra["sites"],
                         "dist_a": ra["dist"], "dist_b": rb["dist"],
                         "d_ll": info.loc[a, "ll"] - info.loc[b, "ll"], "d_lm": info.loc[a, "lm"] - info.loc[b, "lm"]})
    c = pd.DataFrame(rows)
    if c.empty:
        return c
    c["d_dist"] = c.dist_a - c.dist_b
    c["saturated"] = c[["dist_a", "dist_b"]].max(axis=1) > rates.SATURATED_DS
    return c


# --- the tests -------------------------------------------------------------------------------------
def partial_r(t: float, dof: float) -> float:
    """The partial correlation a t statistic with `dof` degrees of freedom stands for: the effect size the bar compares."""
    if not np.isfinite(t) or dof <= 0:
        return np.nan
    return t / math.sqrt(t * t + dof)


def test_table(contrasts: pd.DataFrame, genes: list[str], norms: dict[str, dict], sort: bool) -> pd.DataFrame:
    """H7's loop over classes, saturation, genes and normalisations, with `rates.test_contrasts`. `sort`: H7 orders
    its table by p; the RNA tables are ordered by their keys. Adds q (Benjamini-Hochberg over the genes, as H7), the
    Holm-adjusted p over the same genes (the family-wise correction the bar names) and the partial correlation r."""
    rows = []
    for cls in [*sorted(contrasts["class"].dropna().unique()), "all"]:
        by_class = contrasts if cls == "all" else contrasts[contrasts["class"] == cls]
        for with_saturated in (False, True):
            sub = by_class if with_saturated else by_class[~by_class.saturated]
            for gene in genes:
                g = sub[sub.gene == gene]
                for label, kw in norms.items():
                    r = rates.test_contrasts(g, **kw)
                    if r:
                        rows.append({"class": cls, "gene": gene, "normalisation": label, "includes_saturated": with_saturated, **r})
    res = pd.DataFrame(rows)
    if res.empty:
        return res
    grp = ["class", "normalisation", "includes_saturated"]
    res["q"] = res.groupby(grp).p.transform(lambda p: pd.Series(multipletests(p, method="fdr_bh")[1], index=p.index))
    if not sort:
        res["p_holm"] = res.groupby(grp).p.transform(lambda p: pd.Series(multipletests(p, method="holm")[1], index=p.index))
        cov = res.normalisation.map(lambda n: len(norms[n].get("covariates", [])))
        res["r"] = [partial_r(t, n - 1 - c) for t, n, c in zip(res.t, res.n, cov)]
        return res.sort_values(["class", "normalisation", "includes_saturated", "gene"], kind="mergesort").reset_index(drop=True)
    return res.sort_values("p").reset_index(drop=True)


# --- the species of a run --------------------------------------------------------------------------
def species_for(data: Path, composition: Path, fourfold: Path, tree_file: Path, species_set: str | None) -> pd.DataFrame:
    """H7's species (`rates.select`: the default dataset of H1 with body mass, a tip in the dated tree, all 13 proteins),
    restricted to one of H27's sets when `species_set` is named. The environment variable that restricts every analysis
    at once must be unset: the restriction is made here, once, and visibly."""
    if species_sets.active() is not None:
        raise SystemExit(f"{species_sets.ENV} is set: H37 names its species set itself, unset it")
    sp = rates.select(data, composition, fourfold, tree_file)
    if species_set is not None:
        sp = sp[sp.slug.isin(species_sets.slugs(data, species_set))].reset_index(drop=True)
    return sp


def codon_alignments(data: Path) -> dict[str, dict[str, str]]:
    return {g: read_fasta(data / "alignments" / f"{g}.codon.fasta.gz") for g in MT_GENES}


# --- the control -----------------------------------------------------------------------------------
KEY = ["class", "gene", "normalisation", "includes_saturated"]


def compare(mine: pd.DataFrame, reference: pd.DataFrame) -> dict:
    """Row-by-row comparison of two H7 tables at the precision they are written (five significant digits). -> counts."""
    f = lambda v: f"{v:.5g}" if isinstance(v, (float, np.floating)) else str(v)
    cols = [c for c in reference.columns]
    m = {c: [f(v) for v in mine[c]] for c in cols if c in mine.columns}
    missing_cols = [c for c in cols if c not in mine.columns]
    ref = {c: [f(v) for v in reference[c]] for c in cols}
    key = lambda d, i: tuple(d[k][i] for k in KEY)
    mk = {key(m, i): i for i in range(len(mine))}
    rk = {key(ref, i): i for i in range(len(reference))}
    only_mine, only_ref = sorted(set(mk) - set(rk)), sorted(set(rk) - set(mk))
    differing, cells = 0, 0
    worst = {}
    for k in set(mk) & set(rk):
        bad = [c for c in cols if c in m and m[c][mk[k]] != ref[c][rk[k]]]
        differing += bool(bad)
        cells += len(bad)
        for c in bad:
            worst[c] = worst.get(c, 0) + 1
    return {"rows_mine": len(mine), "rows_reference": len(reference), "only_in_mine": len(only_mine),
            "only_in_reference": len(only_ref), "rows_that_differ": differing, "cells_that_differ": cells,
            "columns_missing_in_mine": missing_cols, "differing_columns": worst,
            "reproduced": not (only_mine or only_ref or differing or missing_cols)}


def control(data: Path, composition: Path, fourfold: Path, tree_file: Path, reference: Path, out: Path) -> dict:
    """The H7 pipeline of this driver on the 13 protein genes, against `reference` (H7's rates.tsv for this build).
    Writes `h7_control_rates.tsv` and `h7_control.json` into `out`; the json is what `require_control` reads."""
    out.mkdir(parents=True, exist_ok=True)
    sp = species_for(data, composition, fourfold, tree_file, None)
    al = codon_alignments(data)
    c = codon_contrasts(sp, tree_file, al)
    if c.empty:
        raise SystemExit("no cherry of the dated tree has two species of the default dataset")
    res = test_table(c, [*MT_GENES, CONCAT], CODON_NORMALISATIONS, sort=True)
    res.to_csv(out / "h7_control_rates.tsv", sep="\t", index=False, float_format="%.5g")
    ref = pd.read_csv(reference, sep="\t")
    cmp = compare(pd.read_csv(out / "h7_control_rates.tsv", sep="\t"), ref)
    verdict = {"species_sha256": _sha(data / "species.parquet"),
               "reference": reference.name, "reference_sha256": _sha(reference), "species_in_default_dataset": len(sp),
               "cherries": int((c.gene == CONCAT).sum()), **cmp}
    (out / CONTROL_FILE).write_text(json.dumps(verdict, indent=2, sort_keys=True) + "\n")
    head = res[(~res.includes_saturated) & (res.gene == CONCAT)].sort_values("p")
    print("concatenated genes, saturated pairs excluded:")
    print(head.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print(json.dumps(verdict, indent=2, sort_keys=True))
    return verdict


def _sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require_control(data: Path, out: Path) -> dict:
    """Refuse unless the H7 control of this build, in `out`, was reproduced. No RNA rate or test is read before."""
    p = out / CONTROL_FILE
    if not p.is_file():
        raise SystemExit(f"no {p}: run `mitoage longevity rna-control` first; the RNA rates are not used before H7 is reproduced")
    v = json.loads(p.read_text())
    if v.get("species_sha256") != _sha(data / "species.parquet"):
        raise SystemExit(f"{p} was made on another build than {data}: run `mitoage longevity rna-control` again")
    if v.get("reproduced") is not True:
        raise SystemExit(f"the pipeline did not reproduce H7 ({p}): {v.get('rows_that_differ')} rows differ, "
                         f"{v.get('only_in_mine')} only here, {v.get('only_in_reference')} only in the reference. "
                         "A rate pipeline that cannot reproduce H7 is not used on the RNA genes.")
    return v


# --- the RNA rates -----------------------------------------------------------------------------------
def rna_alignments(data: Path, genes_dir: Path, genes: list[str]) -> dict[str, dict[str, str]]:
    """The alignment of each gene: the tRNAs from `genes_dir` (made by `rna-genes`), the rRNAs from the build."""
    from .rna_genes import RRNA_GENES, file_name, read_alignment
    return {g: read_alignment(data / "alignments" / f"{g}.nt.fasta.gz") if g in RRNA_GENES
            else read_alignment(genes_dir / f"{file_name(g)}.nt.fasta") for g in genes}


def rna_rates(data: Path, composition: Path, fourfold: Path, tree_file: Path, genes_dir: Path, out: Path,
              genes: list[str] | None = None) -> pd.DataFrame:
    """Rates of the RNA genes, and of the 13 protein genes for the bar's comparison (dS with the body-mass contrast,
    H7's own measure), on H7's species and on H27's `high` set. Writes rna_rate_contrasts.<set>.tsv, rna_rates.tsv
    and protein_rates.tsv (each one table with a column `set`)."""
    require_control(data, out)
    from .rna_genes import RNA_GENES
    genes = genes or RNA_GENES
    al = rna_alignments(data, genes_dir, genes)
    codon = codon_alignments(data)
    tables, ptables = [], []
    for name in ("default", "high"):
        sp = species_for(data, composition, fourfold, tree_file, None if name == "default" else "high")
        c = nt_contrasts(sp, tree_file, al)
        c.insert(0, "set", name)
        c.to_csv(out / f"rna_rate_contrasts.{name}.tsv", sep="\t", index=False, float_format="%.6g")
        t = test_table(c, genes, NT_NORMALISATIONS, sort=False)
        t.insert(0, "set", name)
        tables.append(t)
        pc = codon_contrasts(sp, tree_file, codon)
        pt = test_table(pc[pc.gene != CONCAT], list(MT_GENES), CODON_NORMALISATIONS, sort=False)
        pt.insert(0, "set", name)
        ptables.append(pt)
        print(f"{name}: {len(sp)} species, {len(c)} RNA contrast rows, {len(t)} RNA tests; {len(pt)} protein tests")
    res = pd.concat(tables, ignore_index=True)
    res.to_csv(out / "rna_rates.tsv", sep="\t", index=False, float_format="%.5g")
    pd.concat(ptables, ignore_index=True).to_csv(out / "protein_rates.tsv", sep="\t", index=False, float_format="%.5g")
    return res
