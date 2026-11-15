"""H14: the readouts fixed in advance, per complex, from the tables `burial`, `iface` and `partners` wrote.

Does H5's reading for complex IV, that the coupling of nuclear and mitochondrial subunits belongs to whole proteins
and not to the residues that touch, hold in complexes I, III and V? Per complex (docs/HYPOTHESES.md, H14):

- R0, whole-protein coupling (is there coupling to localise?): `seq_control_other` of `mitoport iface`, the coupling
  of a subunit's non-interface columns with the mitochondrial burial-matched controls, no contact residue on either
  side, for the complex's nuclear subunits against the cytosolic controls; also H5's group contrast,
  `seq_contact_iface`, for the subunits with an interface. Mann-Whitney, two-sided.
- R1a, mitochondrial side, conservation: identity to human at the contact against the control residues (sites.tsv,
  the complex's genes pooled), Mann-Whitney (as H5), and Wilcoxon's signed-rank test over the matched pairs
  (control_pairs.tsv). Descriptive.
- R1b, mitochondrial side, coupling: `seq_excess_mt` of the subunits with at least MIN_IFACE interface columns
  against the cytosolic controls. Mann-Whitney, two-sided; the interface hypothesis predicts a positive shift.
- R2, nuclear side: `seq_excess_nuc`, the same sets and test; predicts a positive shift.
- R3, pair concordance: per species, discordance at the contact pairs minus at the control pairs
  (partners_species.tsv of `mitoport partners`); their mean as the intercept of a phylogenetic regression with
  Pagel's lambda by maximum likelihood on the class's part of the dated tree, with its 95% interval; mammals, and
  birds, ray-finned fishes and reptiles where 25 species or more. Predicts a negative intercept (contacts more
  concordant).
- Holm's correction over R1b, R2 and R3 in mammals, all complexes read together (`p_holm`).

A readout that cannot be computed is written with status "not assessed: <reason>" and no p value; it is never
counted as a negative answer and is left out of the correction.

Writes one table, a row per complex, readout and set, sorted by complex, readout and set.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats
from statsmodels.stats.multitest import multipletests

from ..longevity.aa_pgls import class_covariance
from ..longevity.pgls import gls, nll
from ..longevity.structure_pgls import CLASSES, MIN_SPECIES
from . import complexes
from .iface import MIN_IFACE

MIN_PROTEINS = 3          # a group comparison needs at least this many proteins on each side
CONTROL = "cytosolic control"
HOLM = {"R1b", "R2", "R3"}            # R3 in mammals only
PRIMARY_CLASS = "Mammalia"
COLUMNS = ["complex", "readout", "set", "measure", "status", "n_a", "n_b", "median_a", "median_b", "estimate",
           "ci95_low", "ci95_high", "lambda", "p", "p_holm", "predicted"]


def not_assessed(base: dict, why: str) -> dict:
    return {**base, "status": f"not assessed: {why}"}


def group_test(base: dict, a: pd.Series, b: pd.Series, predicted: str) -> dict:
    """Mann-Whitney (two-sided) of the values of `a` against those of `b`, NaN left out of both."""
    a, b = a.dropna(), b.dropna()
    if len(a) < MIN_PROTEINS or len(b) < MIN_PROTEINS:
        return not_assessed(base, f"{len(a)} and {len(b)} proteins with a value")
    return {**base, "status": "tested", "n_a": len(a), "n_b": len(b), "median_a": a.median(), "median_b": b.median(),
            "estimate": a.median() - b.median(), "p": float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue),
            "predicted": predicted}


def conservation(cx: complexes.Complex, sites: pd.DataFrame, matched: pd.DataFrame) -> list[dict]:
    """R1a: identity to human, contact against control residues, as groups (Mann-Whitney) and as matched pairs."""
    base = {"complex": cx.name, "readout": "R1a", "set": "all species"}
    con = sites.identity_to_human[sites.group.eq("contact")]
    ctl = sites.identity_to_human[sites.group.eq("control")]
    rows = [{**group_test({**base, "measure": "identity to human, contact vs control residues"}, con, ctl, "none"),
             "predicted": "descriptive"}]
    by = sites.set_index(["gene", "aln_col"]).identity_to_human
    d = np.array([by[(r.gene, r.contact_aln_col)] - by[(r.gene, r.control_aln_col)] for r in matched.itertuples()])
    pb = {**base, "measure": "identity to human, contact minus its matched control"}
    if len(d) < MIN_PROTEINS or np.all(d == 0):
        rows.append(not_assessed(pb, f"{len(d)} matched pairs"))
    else:
        rows.append({**pb, "status": "tested", "n_a": len(d), "estimate": float(np.median(d)),
                     "p": float(stats.wilcoxon(d).pvalue), "predicted": "descriptive"})
    return rows


def coupling(cx: complexes.Complex, cp: pd.DataFrame) -> list[dict]:
    """R0, R1b, R2 from coupling.tsv of `mitoport iface`."""
    part = cp[cp.group.eq(cx.partner_group)]
    touch = part[part.iface_cols.ge(MIN_IFACE)]
    ctl = cp[cp.group.eq(CONTROL)]
    base = {"complex": cx.name, "set": "mammals"}
    return [
        group_test({**base, "readout": "R0", "measure": "seq_control_other: every subunit vs cytosolic controls"},
                   part.seq_control_other, ctl.seq_control_other, "positive"),
        group_test({**base, "readout": "R0", "measure": "seq_contact_iface: subunits with an interface vs cytosolic controls"},
                   touch.seq_contact_iface, ctl.seq_contact_iface, "positive"),
        group_test({**base, "readout": "R1b", "measure": "seq_excess_mt: subunits with an interface vs cytosolic controls"},
                   touch.seq_excess_mt, ctl.seq_excess_mt, "positive"),
        group_test({**base, "readout": "R2", "measure": "seq_excess_nuc: subunits with an interface vs cytosolic controls"},
                   touch.seq_excess_nuc, ctl.seq_excess_nuc, "positive"),
    ]


def intercept(C: np.ndarray, y: np.ndarray) -> dict:
    """Phylogenetic regression of y on an intercept alone, Pagel's lambda by maximum likelihood in [0, 1]."""
    X = np.ones((len(y), 1))
    lam = float(optimize.minimize_scalar(nll, bounds=(0.0, 1.0), args=(X, y, C), method="bounded").x)
    b, se, df = gls(X, y, C, lam)
    h = stats.t.ppf(0.975, df) * se[0]
    t = b[0] / se[0]
    return {"lambda": lam, "estimate": float(b[0]), "ci95_low": float(b[0] - h), "ci95_high": float(b[0] + h),
            "p": float(2 * stats.t.sf(abs(t), df))}


def concordance(cx: complexes.Complex, per: pd.DataFrame | None, tree_file: Path) -> list[dict]:
    """R3: per class, the mean of contact minus control discordance over species, on the tree."""
    rows = []
    for cls in CLASSES:
        base = {"complex": cx.name, "readout": "R3", "set": cls,
                "measure": "discordance, contact minus control pairs (phylogenetic mean)"}
        d = per[per["class"].eq(cls)].dropna(subset=["contact minus control"]) if per is not None else pd.DataFrame()
        if len(d) < MIN_SPECIES:
            rows.append(not_assessed(base, f"{len(d)} species with both discordances")); continue
        C, order, _ = class_covariance(tree_file, sorted(d.slug))
        d = d.set_index("slug").loc[order]
        r = intercept(C, d["contact minus control"].to_numpy(dtype=float))
        rows.append({**base, "status": "tested", "n_a": len(d), "median_a": d.contact.mean(), "median_b": d.control.mean(),
                     **r, "predicted": "negative"})
    return rows


def holm(res: pd.DataFrame) -> pd.Series:
    """Holm-adjusted p over the tested rows of the readouts in HOLM (R3 in mammals only), all complexes together."""
    fam = res.status.eq("tested") & res.readout.isin(HOLM) & (res.readout.ne("R3") | res.set.eq(PRIMARY_CLASS))
    out = pd.Series(np.nan, index=res.index)
    if fam.any():
        out[fam] = multipletests(res.p[fam], method="holm")[1]
    return out


def main(names: tuple[str, ...], tree_file: Path, out: Path) -> pd.DataFrame:
    rows = []
    for name in names:
        cx = complexes.get(name)
        need = [cx.work / "sites.tsv", cx.work / "control_pairs.tsv", cx.out / "coupling.tsv"]
        absent = [str(p) for p in need if not p.is_file()]
        if absent:
            raise SystemExit(f"complex {name}: missing {', '.join(absent)}: run `mitoport contacts`, `burial` and `iface` "
                             f"with --complex {name} first")
        rows += conservation(cx, pd.read_csv(cx.work / "sites.tsv", sep="\t", dtype={"chain": str}),
                             pd.read_csv(cx.work / "control_pairs.tsv", sep="\t"))
        rows += coupling(cx, pd.read_csv(cx.out / "coupling.tsv", sep="\t"))
        per_file = cx.partners_out / "partners_species.tsv"
        if not per_file.is_file():
            print(f"NOTE: complex {name}: no {per_file} (`mitoport partners --complex {name}`): R3 NOT ASSESSED")
        rows += concordance(cx, pd.read_csv(per_file, sep="\t") if per_file.is_file() else None, tree_file)
    res = pd.DataFrame(rows).reindex(columns=COLUMNS)
    order = {n: i for i, n in enumerate(complexes.COMPLEXES)}
    res = res.sort_values(["complex", "readout", "set", "measure"], key=lambda s: s.map(order) if s.name == "complex" else s,
                          kind="stable").reset_index(drop=True)
    res["p_holm"] = holm(res)
    out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, sep="\t", index=False, float_format="%.5g")
    with pd.option_context("display.width", 250, "display.max_columns", 30, "display.max_colwidth", 70):
        print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("wrote", out)
    return res
