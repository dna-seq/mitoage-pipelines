"""H19: mito-nuclear coupling class by class, its fall with distance from the mtDNA-encoded core, and what the
mitochondrial side of it is (docs/HYPOTHESES.md, H19).

The statistic is H14's R0 (`seq_control_other` of `mitoport iface`), computed within one class of the build: over
pairs of the class's species, the Spearman correlation of a nuclear protein's distance with the distance of a set of
mitochondrial genes, the patristic distance of the dated tree and its logarithm partialled out. Nuclear side: the
share of differing residues over the columns the human protein fills, the columns that touch a mitochondrial subunit
left out. Mitochondrial side: the burial-matched non-contact control residues of `mitoport burial`, the genes of a
set pooled (mean of the per-gene distances, each scaled to mean 1 over the class's species). No contact residue is
on either side.

Per class, nuclear protein and mitochondrial set (the genes of complex I, III, IV or V, or all 13):
- `coupling`: the statistic above.
- `given_other`: the same with the distance of the other complexes' mitochondrial genes partialled out too, and
  `other_given`: the coupling with the other complexes' genes, this complex's partialled out (sets I to V).
- with all 13 genes: `tv4`, the coupling with the transversion distance at fourfold-degenerate third positions of the
  12 heavy-strand genes (a neutral rate that saturates late); `aa_given_tv4`, the coupling with the mitochondrial
  proteins given that distance; `tv4_given_aa`, the coupling with that distance given the proteins.
An extra covariate enters as its rank over the pairs used.

The distance of a nuclear subunit from the core: over its heavy atoms in the solved human complex, the median
distance to the nearest heavy atom of any mtDNA-encoded chain; a gene in several chains takes the chain with the
smallest median.

Species pairs are not independent, so every test compares groups of proteins (Mann-Whitney over proteins) or ranks
the subunits within their complex and permutes them there. A value that cannot be computed is written empty with
its reason in `status`; it is never read as zero.

Writes core_distance.tsv, class_coupling.tsv, class_matrix.tsv and class_readouts.tsv.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.PDB import MMCIFParser
from scipy import stats
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests

from ..genes import COMPLEX_II, CYTOCHROME_C, CYTOSOLIC_CONTROLS, INNER_MEMBRANE, MATRIX, MT_GENES, MT_HEAVY_STRAND
from ..seqio import read_fasta
from . import complexes
from .iface import MIN_IFACE, hamd, interface_columns, resid

CLASSES = ("Mammalia", "Aves", "Actinopterygii", "Reptilia")
COMPLEXES = ("I", "III", "IV", "V")
ALL = "all"                      # the mitochondrial set of all 13 genes
MIN_SPECIES = 40                 # species a protein needs in a class
MIN_PAIRS = 200                  # usable species pairs a correlation needs (as `iface.pcorr`)
MIN_PROTEINS = 5                 # proteins a group needs on each side of a comparison
MIN_STRATUM = 3                  # subunits a complex needs to enter a stratified statistic
MIN_FOURFOLD = 100               # fourfold sites two species must share
PERMUTATIONS = 10_000
SEED = 20261116
# first two bases of the codon families that are fourfold degenerate in the vertebrate mitochondrial code (table 2)
FOURFOLD_PREFIXES = ("CT", "GT", "TC", "CC", "AC", "GC", "CG", "GG")
CONTROL = "cytosolic control"
OTHER_GROUPS = {"cytochrome c": CYTOCHROME_C, "complex II (no mt subunit)": COMPLEX_II,
                "matrix, not respiratory": MATRIX, "inner membrane, not respiratory": INNER_MEMBRANE,
                CONTROL: CYTOSOLIC_CONTROLS}
PRIMARY = ("C1", "G1", "S1", "T2", "T3")          # the readouts Holm's correction covers, within class
VALUES = ["coupling", "given_other", "other_given", "tv4", "aa_given_tv4", "tv4_given_aa"]
COUPLING_COLUMNS = ["class", "group", "complex", "gene", "mt", "status", "species", "pairs", "columns", "iface_cols",
                    "divergence", *VALUES]
READOUT_COLUMNS = ["class", "readout", "set", "measure", "status", "n_a", "n_b", "median_a", "median_b", "estimate",
                   "p", "p_holm", "predicted"]


def subunit_group(name: str) -> str:
    return f"complex {name} subunit"


# --- the structure ---------------------------------------------------------------------------------
def core_distance(cx: "complexes.Complex") -> pd.DataFrame:
    """Per nuclear gene of the complex: the chain read, the median and the smallest distance (A) of its heavy atoms
    to the nearest heavy atom of any chain of a mitochondrial gene (both copies where the entry holds two). A gene in
    several chains takes the chain with the smallest median (the first in sorted order on a tie)."""
    gene_of = complexes.chain_genes(cx.cif, cx)
    model = MMCIFParser(QUIET=True).get_structure("x", str(cx.cif))[0]
    xyz = {ch.id: np.array([a.coord for a in ch.get_atoms() if a.element != "H"]) for ch in model if ch.id in gene_of}
    core = cKDTree(np.vstack([xyz[c] for c in sorted(xyz) if gene_of[c] in cx.mt]))
    rows = []
    for g in cx.nuclear:
        best = None
        for c in sorted(c for c in xyz if gene_of[c] == g):
            d, _ = core.query(xyz[c])
            cand = (float(np.median(d)), c, float(d.min()))
            if best is None or cand[0] < best[0]:
                best = cand
        rows.append({"complex": cx.name, "gene": g, "chain": best[1], "median_dist": best[0], "min_dist": best[2]})
    return pd.DataFrame(rows)


# --- distances -------------------------------------------------------------------------------------
def fourfold_tv(codon_alns: dict[str, dict[str, str]], sp: list[str]) -> np.ndarray:
    """Share of transversions over the third positions two species share in fourfold-degenerate codons of one family
    (both codons start with the same two bases of FOURFOLD_PREFIXES and end in A, C, G or T), the genes of
    `codon_alns` pooled; NaN where they share fewer than MIN_FOURFOLD such sites."""
    n = len(sp); tv = np.zeros((n, n)); sites = np.zeros((n, n))
    for g, aln in codon_alns.items():
        A = np.array([list(aln[s].upper()) for s in sp])
        if A.shape[1] % 3:
            raise SystemExit(f"{g}: the codon alignment has {A.shape[1]} columns, not a multiple of three")
        first = np.char.add(A[:, 0::3], A[:, 1::3]); third = A[:, 2::3]
        purine, pyrimidine = np.isin(third, ["A", "G"]), np.isin(third, ["C", "T"])
        for prefix in FOURFOLD_PREFIXES:
            fam = first == prefix
            R, Y = (fam & purine).astype(float), (fam & pyrimidine).astype(float)
            both = R + Y
            sites += both @ both.T
            tv += R @ Y.T + Y @ R.T
    return np.where(sites >= MIN_FOURFOLD, tv / np.maximum(sites, 1), np.nan)


def mt_distance(mt_aln: dict[str, dict[str, str]], sp: list[str], cols: dict[str, np.ndarray], genes: list[str]) -> np.ndarray:
    """Distance over species pairs at `cols` of `genes`, pooled as `mitoport iface` pools a complex's genes: the mean
    of the per-gene shares of differing residues, each scaled to mean 1 over `sp`; a gene with under 5 columns is
    left out; all NaN when no gene is left."""
    Ds = []
    for g in genes:
        if len(cols.get(g, ())) < 5:
            continue
        D = hamd(mt_aln[g], sp, cols[g])
        Ds.append(D / np.nanmean(D))
    if not Ds:
        return np.full((len(sp), len(sp)), np.nan)
    with np.errstate(all="ignore"):
        return np.nanmean(np.stack(Ds), axis=0)


def pcorr(x: np.ndarray, y: np.ndarray, T: np.ndarray, extra: tuple[np.ndarray, ...] = ()) -> tuple[float, int]:
    """(Spearman correlation of x and y with T, log(1 + T) and the ranks of each array of `extra` partialled out, the
    pairs used). NaN under MIN_PAIRS usable pairs, or when a side has nothing left after the covariates. Without
    `extra` it is `iface.pcorr`."""
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(T)
    for e in extra:
        ok &= np.isfinite(e)
    n = int(ok.sum())
    if n < MIN_PAIRS:
        return np.nan, n
    cov = [T[ok], np.log1p(T[ok])] + [stats.rankdata(e[ok]) for e in extra]
    rx, ry = resid(stats.rankdata(x[ok]), cov), resid(stats.rankdata(y[ok]), cov)
    if rx.std() < 1e-9 * n or ry.std() < 1e-9 * n:
        return np.nan, n
    return float(np.corrcoef(rx, ry)[0, 1]), n


def shared_paths(tree_file: Path) -> tuple[dict[str, int], np.ndarray]:
    """({tip: index}, shared path length from the root for every pair of tips), as `mitoport iface` reads the tree."""
    tree = Phylo.read(tree_file, "newick")
    tips = {t.name: i for i, t in enumerate(tree.get_terminals())}
    H = np.zeros((len(tips), len(tips)))

    def walk(cl, depth):
        for ch in cl.clades:
            d = depth + (ch.branch_length or 0.0)
            ids = [tips[t.name] for t in ch.get_terminals()]
            H[np.ix_(ids, ids)] = d
            if not ch.is_terminal():
                walk(ch, d)
    walk(tree.root, 0.0)
    return tips, H


def patristic(H: np.ndarray, idx: np.ndarray) -> np.ndarray:
    S = H[np.ix_(idx, idx)]; h = np.diag(S)
    return h[:, None] + h[None, :] - 2 * S


# --- one class -------------------------------------------------------------------------------------
def class_species(species: pd.DataFrame, cls: str, mt_aln: dict, codon_aln: dict, tips: dict, nuc_alns: dict) -> list[str]:
    """The class's species that are in every mitochondrial protein alignment, in every codon alignment read, a tip
    of the tree, and in at least one nuclear alignment; sorted."""
    have = set.intersection(*(set(a) for a in mt_aln.values()), *(set(a) for a in codon_aln.values())) & set(tips)
    nuclear = set().union(*(set(a) for a in nuc_alns.values())) if nuc_alns else set()
    return sorted(set(species.slug[species["class"].eq(cls)]) & have & nuclear)


def couple_class(cls: str, sp: list[str], genes: list[tuple[str, str, str]], nuc_alns: dict, iface_cols: dict,
                 mt_aln: dict, codon_aln: dict, control_cols: dict, tips: dict, H: np.ndarray) -> list[dict]:
    """The rows of class_coupling.tsv for one class: per (group, complex, gene) of `genes` and mitochondrial set."""
    def blank(group, cx_name, g, mt, status, **kw):
        return {"class": cls, "group": group, "complex": cx_name, "gene": g, "mt": mt, "status": status, **kw}

    sets = [*[c for c in COMPLEXES if c in control_cols], ALL]
    if len(sp) < MIN_SPECIES:
        return [blank(group, cx_name, g, mt, f"not assessed: {len(sp)} species of the class with the mitochondrial "
                      "genes, a tip and a nuclear sequence") for group, cx_name, g in genes for mt in sets]
    all_cols = {g: c for cols in control_cols.values() for g, c in cols.items()}
    T_all = patristic(H, np.array([tips[s] for s in sp]))
    D = {c: mt_distance(mt_aln, sp, control_cols[c], list(control_cols[c])) for c in sets if c != ALL}
    D_not = {c: mt_distance(mt_aln, sp, all_cols, [g for g in all_cols if g not in control_cols[c]]) for c in D}
    D[ALL] = mt_distance(mt_aln, sp, all_cols, list(all_cols))
    S = fourfold_tv(codon_aln, sp)
    where = {s: i for i, s in enumerate(sp)}
    rows = []
    for group, cx_name, g in genes:
        aln = nuc_alns.get(g)
        if aln is None or "homo-sapiens" not in aln:
            why = "no alignment" if aln is None else "the alignment holds no human protein"
            rows += [blank(group, cx_name, g, mt, f"not assessed: {why}") for mt in sets]
            continue
        mine = [s for s in sp if s in aln]
        if len(mine) < MIN_SPECIES:
            rows += [blank(group, cx_name, g, mt, f"not assessed: {len(mine)} species", species=len(mine)) for mt in sets]
            continue
        hrow = aln["homo-sapiens"]
        touching = set(int(c) for c in iface_cols.get(g, ()))
        cols = np.array([c for c in range(len(hrow)) if hrow[c] != "-" and c not in touching], dtype=int)
        i = np.array([where[s] for s in mine]); iu = np.triu_indices(len(mine), 1)
        dn = hamd(aln, mine, cols)[iu]
        T = T_all[np.ix_(i, i)][iu]
        base = {"species": len(mine), "columns": len(cols), "iface_cols": len(touching),
                "divergence": float(np.nanmean(dn)) if np.isfinite(dn).any() else np.nan}
        for mt in sets:
            dm = D[mt][np.ix_(i, i)][iu]
            r, n = pcorr(dn, dm, T)
            row = blank(group, cx_name, g, mt, "computed" if np.isfinite(r) else f"not assessed: {n} usable pairs",
                        pairs=n, coupling=r, **base)
            if mt == ALL:
                s4 = S[np.ix_(i, i)][iu]
                row["tv4"] = pcorr(dn, s4, T)[0]
                row["aa_given_tv4"] = pcorr(dn, dm, T, (s4,))[0]
                row["tv4_given_aa"] = pcorr(dn, s4, T, (dm,))[0]
            else:
                do = D_not[mt][np.ix_(i, i)][iu]
                row["given_other"] = pcorr(dn, dm, T, (do,))[0]
                row["other_given"] = pcorr(dn, do, T, (dm,))[0]
            rows.append(row)
        print(f"  {cls}: {group} {g}: {len(mine)} species", flush=True)
    return rows


# --- statistics over proteins ------------------------------------------------------------------------
def within_ranks(v: np.ndarray, strata: np.ndarray) -> np.ndarray:
    """Ranks within each stratum, scaled to (rank - 0.5) / n."""
    out = np.full(len(v), np.nan)
    for s in np.unique(strata):
        m = strata == s
        out[m] = (stats.rankdata(v[m]) - 0.5) / m.sum()
    return out


def _usable(strata: np.ndarray, ok: np.ndarray, need: int) -> np.ndarray:
    keep = ok.copy()
    for s in np.unique(strata):
        if (ok & (strata == s)).sum() < need:
            keep &= strata != s
    return keep


def stratified_corr(x, y, strata, covariate=None, seed: int = SEED) -> dict:
    """Correlation of x and y over items ranked within their stratum (a stratified Spearman): the Pearson correlation
    of the within-stratum scaled ranks, pooled. With `covariate`, both rank vectors are first regressed, stratum by
    stratum, on the covariate's within-stratum rank. Two-sided p by permuting y within strata, PERMUTATIONS times:
    (1 + permutations at least as far from zero) / (1 + permutations). Strata with under MIN_STRATUM usable items
    are left out; {"n": 0} when none is left or a side does not vary."""
    x, y, strata = np.asarray(x, float), np.asarray(y, float), np.asarray(strata)
    ok = np.isfinite(x) & np.isfinite(y)
    if covariate is not None:
        covariate = np.asarray(covariate, float); ok &= np.isfinite(covariate)
    keep = _usable(strata, ok, MIN_STRATUM)
    if not keep.any():
        return {"n": 0}
    x, y, strata = x[keep], y[keep], strata[keep]
    rx, ry = within_ranks(x, strata), within_ranks(y, strata)
    if covariate is not None:
        rc = within_ranks(covariate[keep], strata)
        for s in np.unique(strata):
            m = strata == s
            rx[m], ry[m] = resid(rx[m], [rc[m]]), resid(ry[m], [rc[m]])
    if rx.std() < 1e-12 or ry.std() < 1e-12:
        return {"n": 0}
    r = float(np.corrcoef(rx, ry)[0, 1])
    rng = np.random.default_rng(seed)
    groups = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    hits = 0
    for _ in range(PERMUTATIONS):
        p = ry.copy()
        for gidx in groups:
            p[gidx] = ry[rng.permutation(gidx)]
        hits += abs(np.corrcoef(rx, p)[0, 1]) >= abs(r) - 1e-12
    return {"n": int(len(x)), "strata": len(groups), "r": r, "p": (1 + hits) / (1 + PERMUTATIONS)}


def stratified_diff(v, label, strata, seed: int = SEED) -> dict:
    """Mean within-stratum scaled rank of the items with `label` minus that of the items without, pooled over the
    strata that hold at least two of each; two-sided p by permuting the labels within strata."""
    v, label, strata = np.asarray(v, float), np.asarray(label, bool), np.asarray(strata)
    keep = np.isfinite(v)
    for s in np.unique(strata):
        m = keep & (strata == s)
        if (m & label).sum() < 2 or (m & ~label).sum() < 2:
            keep &= strata != s
    if not keep.any():
        return {"n": 0}
    v, label, strata = v[keep], label[keep], strata[keep]
    r = within_ranks(v, strata)

    def stat(lab):
        return r[lab].mean() - r[~lab].mean()
    obs = float(stat(label))
    rng = np.random.default_rng(seed)
    groups = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    hits = 0
    for _ in range(PERMUTATIONS):
        lab = label.copy()
        for gidx in groups:
            lab[gidx] = label[rng.permutation(gidx)]
        hits += abs(stat(lab)) >= abs(obs) - 1e-12
    return {"n_a": int(label.sum()), "n_b": int((~label).sum()), "strata": len(groups), "diff": obs,
            "median_a": float(np.median(v[label])), "median_b": float(np.median(v[~label])),
            "p": (1 + hits) / (1 + PERMUTATIONS)}


def group_test(base: dict, a: pd.Series, b: pd.Series, predicted: str) -> dict:
    """Mann-Whitney (two-sided) of `a` against `b`, NaN left out; not assessed under MIN_PROTEINS on a side."""
    a, b = a.dropna(), b.dropna()
    if len(a) < MIN_PROTEINS or len(b) < MIN_PROTEINS:
        return {**base, "status": f"not assessed: {len(a)} and {len(b)} proteins with a value"}
    return {**base, "status": "tested", "n_a": len(a), "n_b": len(b), "median_a": a.median(), "median_b": b.median(),
            "estimate": a.median() - b.median(), "p": float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue),
            "predicted": predicted}


def readouts(coupling: pd.DataFrame, dist: pd.DataFrame) -> pd.DataFrame:
    """The readouts of H19 from class_coupling.tsv and core_distance.tsv, class by class; Holm over PRIMARY within class."""
    rows = []
    names = [c for c in COMPLEXES if c in set(coupling.mt)]
    for k, cls in enumerate(c for c in CLASSES if c in set(coupling["class"])):
        t = coupling[coupling["class"].eq(cls)]
        if not t.status.eq("computed").any():
            rows.append({"class": cls, "readout": "*", "set": "*", "measure": "every readout",
                         "status": "not assessed: " + "; ".join(sorted({s.removeprefix('not assessed: ') for s in t.status}))[:300]})
            continue

        def val(mt, col, group=None, where=None):
            d = t[t.mt.eq(mt)]
            if group is not None:
                d = d[d.group.eq(group)]
            if where is not None:
                d = d[where(d)]
            return d.set_index("gene")[col].astype(float)

        subunit = t.complex.isin(names)
        for c in names:                        # C1, S1, S2: one complex at a time, against the cytosolic controls
            for readout, col, what in (("C1", "coupling", "coupling with the complex's mitochondrial genes"),
                                       ("S1", "given_other", "the same, the other complexes' mitochondrial genes partialled out"),
                                       ("S2", "other_given", "coupling with the other complexes' mitochondrial genes, this complex's partialled out")):
                rows.append(group_test({"class": cls, "readout": readout, "set": f"complex {c}",
                                        "measure": f"{what}: the complex's subunits vs cytosolic controls"},
                                       val(c, col, subunit_group(c)), val(c, col, CONTROL), "positive"))
        # the subunits' own coupling, beside their distance from the core
        own = t[subunit & t.mt.eq(t.complex)].merge(dist, on=["complex", "gene"], how="left")
        st = own.complex.to_numpy()
        for readout, cov, what in (("G1", None, "coupling against distance from the mtDNA-encoded core, ranked within complex"),
                                   ("G2", own.divergence.to_numpy(float), "the same, the subunit's own divergence partialled out")):
            base = {"class": cls, "readout": readout, "set": "subunits of " + ", ".join(names), "measure": what}
            got = stratified_corr(own.median_dist, own.coupling, st, cov, seed=SEED + 10 * k + (cov is not None))
            rows.append({**base, "status": "not assessed: no complex with enough subunits"} if got["n"] == 0 else
                        {**base, "status": "tested", "n_a": got["n"], "n_b": got["strata"], "estimate": got["r"],
                         "p": got["p"], "predicted": "negative"})
        base = {"class": cls, "readout": "G3", "set": "subunits of " + ", ".join(names),
                "measure": f"coupling of subunits with {MIN_IFACE} or more interface columns vs subunits with none, ranked within complex"}
        two = own[(own.iface_cols >= MIN_IFACE) | (own.iface_cols == 0)]
        got = stratified_diff(two.coupling, two.iface_cols >= MIN_IFACE, two.complex.to_numpy(), seed=SEED + 10 * k + 2)
        rows.append({**base, "status": "not assessed: no complex with two subunits of each kind"} if "p" not in got else
                    {**base, "status": "tested", "n_a": got["n_a"], "n_b": got["n_b"], "median_a": got["median_a"],
                     "median_b": got["median_b"], "estimate": got["diff"], "p": got["p"], "predicted": "positive"})
        # groups against each other, all 13 mitochondrial genes
        touch = lambda d: d.complex.isin(names) & (d.iface_cols >= MIN_IFACE)      # noqa: E731
        apart = lambda d: d.complex.isin(names) & (d.iface_cols == 0)              # noqa: E731
        im = "inner membrane, not respiratory"
        for readout, a, b, what, predicted in (
                ("M1", val(ALL, "coupling", im), val(ALL, "coupling", CONTROL), "inner-membrane proteins without a mitochondrial partner vs cytosolic controls", "none"),
                ("M2", val(ALL, "coupling", where=touch), val(ALL, "coupling", im), "subunits that touch a mitochondrial one vs inner-membrane proteins without a partner", "positive"),
                ("O1", val(ALL, "coupling", "complex II (no mt subunit)"), val(ALL, "coupling", CONTROL), "complex II vs cytosolic controls", "none"),
                ("O2", val(ALL, "coupling", "matrix, not respiratory"), val(ALL, "coupling", CONTROL), "matrix enzymes vs cytosolic controls", "none"),
                ("O3", val(ALL, "coupling", where=apart), val(ALL, "coupling", CONTROL), "subunits that touch no mitochondrial one vs cytosolic controls", "none"),
                ("T1", val(ALL, "tv4", where=lambda d: d.complex.isin(names)), val(ALL, "tv4", CONTROL), "coupling with the fourfold transversion distance: all subunits vs cytosolic controls", "positive"),
                ("T2", val(ALL, "aa_given_tv4", where=lambda d: d.complex.isin(names)), val(ALL, "aa_given_tv4", CONTROL), "coupling with the mitochondrial proteins given the fourfold transversion distance: all subunits vs cytosolic controls", "positive"),
                ("T3", val(ALL, "tv4_given_aa", where=lambda d: d.complex.isin(names)), val(ALL, "tv4_given_aa", CONTROL), "coupling with the fourfold transversion distance given the mitochondrial proteins: all subunits vs cytosolic controls", "positive")):
            rows.append(group_test({"class": cls, "readout": readout, "set": "all 13 mitochondrial genes",
                                    "measure": f"{what}"}, a, b, predicted))
    res = pd.DataFrame(rows).reindex(columns=READOUT_COLUMNS)
    for cls, idx in res[res.status.eq("tested") & res.readout.isin(PRIMARY)].groupby("class").groups.items():
        res.loc[idx, "p_holm"] = multipletests(res.loc[idx, "p"], method="holm")[1]
    return res


def matrix(coupling: pd.DataFrame) -> pd.DataFrame:
    """Median coupling of each group of proteins with each mitochondrial set, per class (descriptive)."""
    d = coupling[coupling.status.eq("computed")]
    rows = [{"class": cls, "group": grp, "mt": mt, "proteins": len(x), "median_coupling": x.coupling.median(),
             "median_given_other": x.given_other.median(), "median_other_given": x.other_given.median()}
            for (cls, grp, mt), x in d.groupby(["class", "group", "mt"], sort=True)]
    return pd.DataFrame(rows, columns=["class", "group", "mt", "proteins", "median_coupling", "median_given_other",
                                       "median_other_given"])


# --- the run ---------------------------------------------------------------------------------------
def main(data: Path, nucaln: Path, tree_file: Path, out: Path, classes: tuple[str, ...] = CLASSES,
         names: tuple[str, ...] = COMPLEXES, work: dict[str, Path] | None = None) -> pd.DataFrame:
    """`names`: the complexes read; `work`: {complex: folder of its contact map} when not the default places."""
    cxs = [complexes.get(n) for n in names]
    work = work or {cx.name: cx.work for cx in cxs}
    out.mkdir(parents=True, exist_ok=True)
    genes, iface_cols, control_cols, dists = [], {}, {}, []
    for cx in cxs:
        w = work[cx.name]
        for f in ("sites.tsv", "chains.tsv", "contacts.tsv"):
            if not (w / f).is_file():
                raise SystemExit(f"needs complex {cx.name}'s contact map ({w / f}): run `mitoport contacts --complex "
                                 f"{cx.name}` and `mitoport burial --complex {cx.name}` first")
        sites = pd.read_csv(w / "sites.tsv", sep="\t", dtype={"chain": str})
        unplaced = sites[sites.group.eq("control") & sites.aln_col.le(0)]
        if len(unplaced):
            raise SystemExit(f"complex {cx.name}: {len(unplaced)} control residues have no alignment column in "
                             f"{w / 'sites.tsv'}; refusing to read them as the alignment's last column")
        control_cols[cx.name] = {g: sites[sites.gene.eq(g) & sites.group.eq("control")].aln_col.to_numpy() - 1 for g in cx.mt}
        ct = pd.read_csv(w / "contacts.tsv", sep="\t", dtype={"chain": str, "partner_chain": str})
        model = MMCIFParser(QUIET=True).get_structure("x", str(cx.cif))[0]
        chains = pd.read_csv(w / "chains.tsv", sep="\t", dtype={"chain": str})
        headless = [g for g in cx.nuclear if (nucaln / f"{g}.afa").is_file() and "homo-sapiens" not in read_fasta(nucaln / f"{g}.afa")]
        if headless:      # its columns cannot be placed on the structure: written as not assessed below
            print(f"NOTE: complex {cx.name}: no human protein in the alignment of {' '.join(headless)}")
        _, cols, _ = interface_columns(model, chains[~chains.gene.isin(headless)], ct[ct.partner_is_mt == 0], nucaln)
        iface_cols.update(cols)
        genes += [(subunit_group(cx.name), cx.name, g) for g in cx.nuclear]
        dists.append(core_distance(cx))
    seen = {g for _, _, g in genes}
    twice = [g for g in seen if sum(g == x for _, _, x in genes) > 1]
    if twice:
        raise SystemExit(f"a gene is a subunit of two complexes read: {' '.join(sorted(twice))}")
    for group, members in OTHER_GROUPS.items():
        genes += [(group, "", g) for g in members if g not in seen]
    dist = pd.concat(dists, ignore_index=True)
    dist.to_csv(out / "core_distance.tsv", sep="\t", index=False, float_format="%.3f")

    nuc_alns = {g: read_fasta(nucaln / f"{g}.afa") for _, _, g in genes if (nucaln / f"{g}.afa").is_file()}
    missing = sorted(g for _, _, g in genes if g not in nuc_alns)
    if not nuc_alns:
        raise SystemExit(f"no alignment under {nucaln}: `download orthologs --ortholog-taxon vertebrata`, `mitoport "
                         "select --class ...` and `mitoport iface-align` make them")
    if missing:
        print(f"NOTE: no alignment for {len(missing)} genes, written as not assessed: {' '.join(missing)}")
    mt_aln = {g: read_fasta(data / "alignments" / f"{g}.aa.fasta.gz") for g in MT_GENES}
    codon_aln = {g: read_fasta(data / "alignments" / f"{g}.codon.fasta.gz") for g in MT_HEAVY_STRAND}
    tips, H = shared_paths(tree_file)
    species = pd.read_parquet(data / "species.parquet")
    rows = []
    for cls in classes:
        sp = class_species(species, cls, mt_aln, codon_aln, tips, nuc_alns)
        print(f"{cls}: {len(sp)} species with the 13 mitochondrial genes, a tip and a nuclear sequence", flush=True)
        rows += couple_class(cls, sp, genes, nuc_alns, iface_cols, mt_aln, codon_aln, control_cols, tips, H)
    coupling = pd.DataFrame(rows).reindex(columns=COUPLING_COLUMNS)
    coupling.to_csv(out / "class_coupling.tsv", sep="\t", index=False, float_format="%.4f")
    matrix(coupling).to_csv(out / "class_matrix.tsv", sep="\t", index=False, float_format="%.4f")
    res = readouts(coupling, dist)
    res.to_csv(out / "class_readouts.tsv", sep="\t", index=False, float_format="%.5g")
    with pd.option_context("display.width", 250, "display.max_colwidth", 90, "display.max_rows", 500):
        print(res.drop(columns=["predicted"]).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("wrote", out / "core_distance.tsv", out / "class_coupling.tsv", out / "class_matrix.tsv", out / "class_readouts.tsv")
    return res
