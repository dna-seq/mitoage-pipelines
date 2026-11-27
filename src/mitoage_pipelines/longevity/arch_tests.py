"""H38, steps 2 and 3: the positive control of the rearrangement measure, and the architecture measures against maximum
lifespan under the controls H11 and H27 settled, with the bar of the entry.

The positive control (`control`, written to `arch_control.json`; the tests are refused without it, as H37's `rna-control`
gates `rna-rates`). The rearrangement measure must agree with H13's audit of the same records where the two overlap. The
audit speaks of a record against *its family* (a record is flagged when it is the one record of a family of three or more
with an arrangement of its own), so the control reads the measure's distance against the family's most common
arrangement, with the same function and the same arrangement strings the lifespan tests use against the class's:
  G0  the arrangements of `architecture.tsv` are the audit's `order.tsv`, record by record, at both tiers;
  G1  every record the audit flags as singly arranged, other than a reversed deposit, is at a distance above zero;
  G2  every record the audit flags as a reversed deposit is at zero (the arrangement is canonical under reversal: the
      record is the other strand of its family's arrangement, which is not a rearrangement);
  G3  every record of every family the audit finds uniform (three or more records, one arrangement) is at zero;
  G4  where a family's most common arrangement is its class's, a record's distance to the class equals its distance
      to the family: the measure used against lifespan is the audit's reading there, to the digit.
Every one is checked with no exception, at both tiers. Departures that are not failures are counted, not hidden: a family
whose own arrangement is not its class's (a marsupial family against the mammals, say) is at a distance above zero from the
class and zero from itself; its records are listed in `arch_control_class.tsv`.

The tests (`run_cells`). Per taxon group (Mammalia, Aves, Actinopterygii, Reptilia: the within-group control), per
measure, on the species of H1's default dataset (`rna_tests.test_species`: a lifespan, a body mass, all 13 proteins, a
tip of the dated tree) in two sets: `not_unreliable`, H27's confidence-filtered records, which are the primary form, and
`high`, H27's high-confidence subset. A phylogenetic regression on the dated tree (`aa_pgls.fit`) of the standardised
measure on log10 lifespan and log10 adult mass, with Pagel's lambda estimated, with lambda 1 (the regression equal to
phylogenetically independent contrasts: the primary form here), and with lambda 0 (the raw regression, beside them for
comparison). Each of them again with the annotation flag in the model: 1 for a species whose record the audit flags in any
of four ways (a singly arranged record of its family at either tier; strand reversed or mixed; neither genetic code clean;
an annotation of the build that differs from the record's). A cell the flag cannot enter (no flagged species, or only flagged
ones, among the species of the cell) has no covariate row.

The family-wise correction (`correct_run`) is over every row of the whole run: every measure, taxon group, species set,
model and covariate variant, Holm and Benjamini-Hochberg, in the same function that returns the rows and before any
verdict reads them: `p_holm_run` and `q_run` are columns of the table the p-values are in. `judge` refuses a table without them.

The bar (`judge`), per measure and taxon group:
  (1) survives the controls: raw p below 0.05 with lambda estimated and with lambda 1, the same sign;
  (2) holds on the high-confidence subset: the same on the `high` set, with the sign of the primary result;
  (3) survives the correction over the whole run: `p_holm_run` below 0.05 in both lambda models, on the primary set and on
      the high set;
  (4) does not move when the annotation flag is added: in both lambda models `p_holm_run` of the row with the flag is below
      0.05 and the slope changes by no more than MOVE_TOLERANCE (a quarter, so less than the slope itself: a change of sign is a
      move) of itself.
Verdicts: `counts` (all four); `explained by annotation` (1 to 3 and not 4); `not counting` (1 and not 2 or not 3);
`null` (not 1); `not assessed` where a condition that matters could not be evaluated (a set too small, or no covariate row),
with the reason. A measure counts only where one of its taxon groups counts.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from . import aa_pgls, species_sets
from . import arch_measures as am
from .rna_rates import partial_r
from .rna_tests import test_species
from ..audit import order as aorder

CLASSES = ["Mammalia", "Aves", "Actinopterygii", "Reptilia"]
SETS = ("not_unreliable", "high")
MEASURES = ["rearr_full", "rearr_genes", "cr_length", "cr_copies", "repeat_genome", "repeat_cr", "non_coding_bp",
            "non_coding_frac", "gene_overlap_bp"]
LOGGED = ("cr_length", "non_coding_bp")
MODELS = (("lambda_ml", None), ("lambda_1", 1.0), ("raw", 0.0))
MIN_SPECIES = 25
MIN_MINORITY = 5                 # species not at the measure's most common value, or the measure is almost constant
MOVE_TOLERANCE = 0.25
ALPHA = 0.05
CONTROL_FILE = "arch_control.json"
FLAG_COLUMNS = ["flag_order", "flag_strand", "flag_code", "flag_build"]


# --- the annotation flag -------------------------------------------------------------------------------
def annotation_flags(audit: Path) -> pd.DataFrame:
    """The audit's per-record flags (`audit order|strand|code`, `audit records`), one row per species of the audit's
    tables, with the four reasons and `annot_flag`, their any. A species the strand, code and annotation tables do not
    all have is refused: a species with no audit row is never taken as unflagged."""
    order = pd.read_csv(audit / "order.tsv", sep="\t", usecols=["slug"])
    flags = pd.read_csv(audit / "order_flags.tsv", sep="\t")
    strand = pd.read_csv(audit / "strand.tsv", sep="\t", usecols=["slug", "orientation"])
    code = pd.read_csv(audit / "code.tsv", sep="\t", usecols=["slug", "neither_clean"])
    build = pd.read_csv(audit / "annotation_vs_build.tsv", sep="\t", usecols=["slug", "build_annotation"], keep_default_na=False)
    slugs = set(order.slug)
    for name, d in (("strand", strand), ("code", code), ("annotation_vs_build", build)):
        if set(d.slug) != slugs:
            raise ValueError(f"the audit's {name}.tsv and order.tsv do not hold the same species")
    out = order[["slug"]].copy()
    out["flag_order"] = out.slug.isin(set(flags.slug))
    out["flag_strand"] = out.slug.isin(set(strand.slug[strand.orientation.isin(["reversed", "mixed"])]))
    out["flag_code"] = out.slug.isin(set(code.slug[code.neither_clean == "yes"]))
    out["flag_build"] = out.slug.isin(set(build.slug[build.build_annotation == "differs"]))
    out["annot_flag"] = out[FLAG_COLUMNS].any(axis=1).astype(int)
    return out.sort_values("slug", kind="mergesort").reset_index(drop=True)


# --- the positive control ------------------------------------------------------------------------------
def control_tables(arch: pd.DataFrame, audit_order: pd.DataFrame, audit_flags: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """The checks G0 to G4 for both tiers. -> (the one-row-per-check table, the class-level departures, the counts)."""
    rows, depart, counts = [], [], {}
    merged = arch.merge(audit_order, on="slug", how="left", suffixes=("", "_audit"), indicator=True)
    exceptions = int((merged._merge != "both").sum())              # a live species the audit has no row for
    both = merged[merged._merge == "both"]
    for tier in am.TIERS:
        a, b = both[f"order_{tier}_unoriented"], both[f"order_{tier}_unoriented_audit"]
        exceptions += int(((a != b) & ~(a.isna() & b.isna())).sum())
    rows.append({"check": "G0 the arrangements are the audit's", "tier": "both", "records": int(len(merged)), "exceptions": exceptions})
    for tier in am.TIERS:
        col = f"order_{tier}"
        known = audit_order[audit_order[col].notna()].copy()
        fam = am.class_distances(known.reset_index(drop=True), tier, "family")
        known = known.reset_index(drop=True).join(fam)
        fam_d = known.set_index("slug")[f"rearr_{tier}"].astype(float)
        fam_modal = known.set_index("slug")[f"modal_{tier}"]
        cls = arch.set_index("slug")[[f"rearr_{tier}", f"modal_{tier}"]]
        fl = audit_flags[audit_flags.tier == tier]
        rev, odd = fl[fl.kind == "reversed deposit"], fl[fl.kind != "reversed deposit"]
        d1 = fam_d.reindex(odd.slug)
        d2 = fam_d.reindex(rev.slug)
        rows.append({"check": "G1 a record flagged as singly arranged is above zero from its family", "tier": tier,
                     "records": len(odd), "exceptions": int((~(d1 > 0)).sum())})
        rows.append({"check": "G2 a reversed deposit is at zero", "tier": tier, "records": len(rev), "exceptions": int((~(d2 == 0)).sum())})
        uni = []
        for f, d in known[known.family.notna()].groupby("family", sort=True):
            if len(d) >= aorder.MIN_FAMILY and d[col].nunique() == 1:
                uni.append((f, list(d.slug)))
        uni_slugs = [s for _, ss in uni for s in ss]
        rows.append({"check": "G3 a record of a uniform family is at zero from its family", "tier": tier,
                     "records": len(uni_slugs), "exceptions": int((fam_d.reindex(uni_slugs) != 0).sum()), "families": len(uni)})
        # G4: where the family's arrangement is the class's, the two distances agree
        involved = sorted(set(uni_slugs) | set(fl.slug))
        same, differ, g4_bad = 0, 0, 0
        for s in involved:
            if s not in cls.index or pd.isna(cls.at[s, f"rearr_{tier}"]) or s not in fam_d.index:
                continue
            famrow_modal = fam_modal.at[s]
            if famrow_modal == cls.at[s, f"modal_{tier}"]:
                same += 1
                g4_bad += int(float(fam_d.at[s]) != float(cls.at[s, f"rearr_{tier}"]))
            else:
                differ += 1
                depart.append({"tier": tier, "slug": s, "in": "uniform family" if s in set(uni_slugs) else "flagged",
                               "family_distance": float(fam_d.at[s]), "class_distance": float(cls.at[s, f"rearr_{tier}"])})
        rows.append({"check": "G4 where the family's arrangement is its class's, class distance equals family distance", "tier": tier,
                     "records": same, "exceptions": g4_bad})
        counts[tier] = {"flagged_records": len(fl), "reversed_deposits": len(rev), "uniform_families": len(uni),
                        "uniform_family_records": len(uni_slugs), "records_in_families_that_differ_from_their_class": differ,
                        "flagged_at_zero_from_class": int(sum(1 for r in depart if r["tier"] == tier and r["in"] == "flagged" and r["class_distance"] == 0)),
                        "uniform_family_records_above_zero_from_class": int(sum(1 for r in depart if r["tier"] == tier and r["in"] == "uniform family" and r["class_distance"] > 0))}
    return pd.DataFrame(rows), pd.DataFrame(depart, columns=["tier", "slug", "in", "family_distance", "class_distance"]), counts


def compare_published(flags: pd.DataFrame, published: pd.DataFrame) -> dict:
    """The flags of this run's audit against the flags H13 published (an earlier build): by tier, how many are in both, only in
    the published table, only here. Said, not gated: the two builds are not the same records."""
    out = {}
    for tier in am.TIERS:
        a = set(flags.slug[flags.tier == tier])
        b = set(published.slug[published.tier == tier])
        out[tier] = {"here": len(a), "published": len(b), "in_both": len(a & b), "only_published": len(b - a), "only_here": len(a - b)}
    return out


def control(arch_file: Path, audit: Path, species_file: Path, out: Path, published: Path | None = None) -> dict:
    """Runs the checks, writes arch_control.tsv, arch_control_class.tsv and arch_control.json; the verdict says whether every
    check holds with no exception. `published`: H13's order_flags.tsv of an earlier build, compared and not gated."""
    arch = pd.read_csv(arch_file, sep="\t")
    order = pd.read_csv(audit / "order.tsv", sep="\t")
    flags = pd.read_csv(audit / "order_flags.tsv", sep="\t")
    tab, depart, counts = control_tables(arch, order, flags)
    out.mkdir(parents=True, exist_ok=True)
    tab.to_csv(out / "arch_control.tsv", sep="\t", index=False)
    depart.sort_values(["tier", "slug"]).to_csv(out / "arch_control_class.tsv", sep="\t", index=False, float_format="%.6g")
    bad = tab[tab.exceptions > 0]
    verdict = {"reproduced": bool(len(bad) == 0), "checks": int(len(tab)), "checks_with_exceptions": int(len(bad)),
               "species_sha256": am.sha256_of(species_file), "architecture_sha256": am.sha256_of(arch_file),
               "class_level": counts}
    if published is not None:
        verdict["against_published_flags"] = compare_published(flags, pd.read_csv(published, sep="\t"))
    (out / CONTROL_FILE).write_text(json.dumps(verdict, indent=2, sort_keys=True) + "\n")
    print(tab.to_string(index=False))
    print("the positive control " + ("reproduced the audit's flags" if verdict["reproduced"] else "FAILED: the measure is not used against lifespan"))
    return verdict


def require_control(data: Path, out: Path) -> dict:
    """Refuse unless the control of this build, in `out`, reproduced the audit's flags. No lifespan number is read before."""
    p = out / CONTROL_FILE
    if not p.is_file():
        raise SystemExit(f"no {p}: run `mitoage longevity arch-control` first; the architecture is not tested before the audit's flags are reproduced")
    v = json.loads(p.read_text())
    if v.get("species_sha256") != am.sha256_of(data / "species.parquet"):
        raise SystemExit(f"{p} was made on another build than {data}: run `mitoage longevity arch-control` again")
    if v.get("reproduced") is not True:
        raise SystemExit(f"the rearrangement measure did not reproduce the audit's flags ({p}: {v.get('checks_with_exceptions')} checks "
                         "with exceptions): it is not used against lifespan")
    return v


# --- the species and the measures ----------------------------------------------------------------------
def prepare(arch: pd.DataFrame, flags: pd.DataFrame, sp: pd.DataFrame) -> pd.DataFrame:
    """The test species (`rna_tests.test_species`) with the measures as the regressions read them (lengths in log10) and the
    annotation flag; the set columns `not_unreliable` and `high`. A test species with no architecture row or no audit row is
    refused."""
    a = arch.copy()
    for c in LOGGED:
        a[c] = np.log10(a[c].where(a[c] > 0))
    missing = sorted(set(sp.slug) - set(a.slug))
    if missing:
        raise ValueError(f"{len(missing)} test species have no architecture row, first {missing[0]}")
    unaudited = sorted(set(sp.slug) - set(flags.slug))
    if unaudited:
        raise ValueError(f"{len(unaudited)} test species have no audit row, first {unaudited[0]}")
    d = sp.drop(columns=[c for c in sp.columns if c in MEASURES]).merge(a[["slug", *MEASURES]], on="slug", how="left", validate="one_to_one") \
          .merge(flags[["slug", "annot_flag"]], on="slug", how="left", validate="one_to_one")
    d["not_unreliable"] = species_sets.members(d, "not_unreliable")
    d["high"] = species_sets.members(d, "high")
    return d


# --- the regressions -----------------------------------------------------------------------------------
def fit_one(C: np.ndarray, y: np.ndarray, ll: np.ndarray, lm: np.ndarray, flag: np.ndarray | None, lam: float | None) -> dict:
    """The standardised measure on log10 lifespan and log10 mass (and the flag); the lifespan slope."""
    z = (y - y.mean()) / y.std(ddof=1)
    cols = [ll, lm] if flag is None else [ll, lm, flag.astype(float)]
    r = aa_pgls.fit(C, z, np.column_stack(cols), lam)
    r["r"] = partial_r(r["t"], r["df"])
    return r


def run_cells(sp: pd.DataFrame, tree_file: Path, measures: list[str] | None = None, classes: list[str] | None = None) -> pd.DataFrame:
    """Every cell of the run, corrected over the whole family before it is returned. `sp`: `prepare`'s table."""
    measures = MEASURES if measures is None else measures
    classes = CLASSES if classes is None else classes
    rows = []
    for name in SETS:
        for cls in classes:
            cell = sp[(sp["class"] == cls) & sp[name]]
            if len(cell) < MIN_SPECIES:
                for m in measures:
                    for cov in (False, True):
                        for model, _ in MODELS:
                            rows.append({"set": name, "class": cls, "measure": m, "covariate": cov, "model": model, "n": len(cell),
                                         "reason": f"only {len(cell)} species in the set"})
                continue
            C, order, extra = aa_pgls.class_covariance(tree_file, list(cell.slug))
            cell = cell.set_index("slug").loc[order]
            print(f"{cls} {name}: {len(cell)} species ({extra} other tips under their common ancestor)", flush=True)
            for m in measures:
                y_all = cell[m].to_numpy(dtype=float)
                ok = np.isfinite(y_all)
                n = int(ok.sum())
                y = y_all[ok]
                common = int(pd.Series(y).value_counts().max()) if n else 0
                base = {"set": name, "class": cls, "measure": m, "n": n, "n_off_mode": n - common}
                if n < MIN_SPECIES:
                    reason = f"only {n} species with the measure"
                elif n - common < MIN_MINORITY:
                    reason = f"only {n - common} species off the most common value"
                else:
                    reason = ""
                flag = cell.annot_flag.to_numpy(dtype=float)[ok]
                base["n_flagged"] = int(flag.sum())
                estimable = 0 < flag.sum() < n
                ll, lm = cell.ll.to_numpy()[ok], cell.lm.to_numpy()[ok]
                Cg = C[np.ix_(np.where(ok)[0], np.where(ok)[0])]
                for cov in (False, True):
                    for model, lam in MODELS:
                        row = {**base, "covariate": cov, "model": model, "reason": reason}
                        if not reason and cov and not estimable:
                            row["reason"] = "the annotation flag cannot enter: no flagged species, or all flagged"
                        elif not reason:
                            r = fit_one(Cg, y, ll, lm, flag if cov else None, lam)
                            row.update({"lambda": r["lambda"], "slope": r["slope"], "se": r["se"], "t": r["t"], "df": r["df"],
                                        "p": r["p"], "r": r["r"]})
                        rows.append(row)
    return correct_run(pd.DataFrame(rows))


def correct_run(cells: pd.DataFrame) -> pd.DataFrame:
    """Holm and Benjamini-Hochberg over every row of the run that has a p-value: the whole family at once, the columns of
    the table the raw p are in."""
    cells = cells.copy()
    for c in ("lambda", "slope", "se", "t", "df", "p", "r"):          # a run in which no cell could be fitted still has the columns
        if c not in cells.columns:
            cells[c] = np.nan
    cells["p_holm_run"] = np.nan
    cells["q_run"] = np.nan
    have = cells.p.notna()
    cells.loc[have, "p_holm_run"] = multipletests(cells.p[have], method="holm")[1]
    cells.loc[have, "q_run"] = multipletests(cells.p[have], method="fdr_bh")[1]
    cells["family_size"] = int(have.sum())
    return cells.sort_values(["measure", "class", "set", "covariate", "model"], kind="mergesort").reset_index(drop=True)


# --- the bar -------------------------------------------------------------------------------------------
def _pair(cells: pd.DataFrame, measure: str, cls: str, name: str, cov: bool):
    d = cells[(cells.measure == measure) & (cells["class"] == cls) & (cells.set == name) & (cells.covariate == cov)].set_index("model")
    if not all(m in d.index and pd.notna(d.at[m, "p"]) for m in ("lambda_ml", "lambda_1")):
        return None
    return d


def judge(cells: pd.DataFrame) -> pd.DataFrame:
    """The bar, one row per measure and taxon group; see the module's text. Needs the columns of `correct_run`."""
    for c in ("p_holm_run", "q_run"):
        if c not in cells.columns:
            raise ValueError(f"judge reads the whole-run correction in {c!r}: it is computed with the p-values, by run_cells")
    out = []
    for (m, cls), g in cells.groupby(["measure", "class"], sort=True):
        prim, high = _pair(cells, m, cls, "not_unreliable", False), _pair(cells, m, cls, "high", False)
        cov = _pair(cells, m, cls, "not_unreliable", True)
        row = {"measure": m, "class": cls, "n": np.nan, "n_flagged": np.nan, "verdict": "", "reason": ""}
        if prim is None:
            why = g[(g.set == "not_unreliable") & g.reason.notna() & g.reason.ne("")].reason
            row.update(verdict="not assessed", reason=why.iloc[0] if len(why) else "no regression")
            out.append(row)
            continue
        ml, l1 = prim.loc["lambda_ml"], prim.loc["lambda_1"]
        raw = prim.loc["raw"] if "raw" in prim.index else None
        row.update(n=int(ml.n), n_flagged=int(ml.n_flagged), slope_ml=ml.slope, slope_pic=l1.slope, r_ml=ml.r, r_pic=l1.r,
                   lambda_ml=ml["lambda"], p_ml=ml.p, p_pic=l1.p, p_holm_run_ml=ml.p_holm_run, p_holm_run_pic=l1.p_holm_run,
                   q_run_ml=ml.q_run, q_run_pic=l1.q_run,
                   raw_slope=np.nan if raw is None else raw.slope, raw_p=np.nan if raw is None else raw.p,
                   raw_p_holm_run=np.nan if raw is None else raw.p_holm_run)
        sgn = np.sign(ml.slope)
        c1 = bool(ml.p < ALPHA and l1.p < ALPHA and sgn == np.sign(l1.slope))
        c2 = None if high is None else bool(all(high.at[k, "p"] < ALPHA and np.sign(high.at[k, "slope"]) == sgn for k in ("lambda_ml", "lambda_1")))
        c3_prim = bool(ml.p_holm_run < ALPHA and l1.p_holm_run < ALPHA)
        c3 = None if high is None else bool(c3_prim and all(high.at[k, "p_holm_run"] < ALPHA for k in ("lambda_ml", "lambda_1")))
        if high is not None:
            row.update(high_n=int(high.at["lambda_ml", "n"]), high_p_ml=high.at["lambda_ml", "p"], high_p_pic=high.at["lambda_1", "p"],
                       high_p_holm_run_ml=high.at["lambda_ml", "p_holm_run"], high_p_holm_run_pic=high.at["lambda_1", "p_holm_run"],
                       high_slope_ml=high.at["lambda_ml", "slope"], high_slope_pic=high.at["lambda_1", "slope"])
        c4 = None
        if cov is not None:
            moved = []
            for k in ("lambda_ml", "lambda_1"):
                s0, s1 = prim.at[k, "slope"], cov.at[k, "slope"]
                moved.append(not (cov.at[k, "p_holm_run"] < ALPHA and abs(s1 - s0) <= MOVE_TOLERANCE * abs(s0)))
            c4 = not any(moved)
            row.update(cov_slope_ml=cov.at["lambda_ml", "slope"], cov_slope_pic=cov.at["lambda_1", "slope"],
                       cov_p_holm_run_ml=cov.at["lambda_ml", "p_holm_run"], cov_p_holm_run_pic=cov.at["lambda_1", "p_holm_run"])
        row.update(c1_survives_controls=c1, c2_holds_on_high=c2, c3_survives_run_correction=c3, c4_unmoved_by_flag=c4)
        if not c1:
            row["verdict"] = "null"
        elif c2 is False or c3 is False:
            row["verdict"] = "not counting"
        elif c2 is None or c3 is None:
            row.update(verdict="not assessed", reason="no regression on the high-confidence subset")
        elif c4 is None:
            row.update(verdict="not assessed", reason="the annotation flag cannot enter this cell")
        elif not c4:
            row["verdict"] = "explained by annotation"
        else:
            row["verdict"] = "counts"
        out.append(row)
    return pd.DataFrame(out)


def summarise(v: pd.DataFrame) -> pd.DataFrame:
    """One row per measure: the verdicts over its taxon groups; the measure counts where one of its groups counts."""
    rows = []
    for m in MEASURES:
        g = v[v.measure == m]
        row = {"measure": m, "groups": len(g)}
        for name in ("counts", "explained by annotation", "not counting", "null", "not assessed"):
            row["n_" + name.replace(" ", "_")] = int((g.verdict == name).sum())
        row["measure_counts"] = bool((g.verdict == "counts").any())
        row["counts_in"] = ", ".join(g["class"][g.verdict == "counts"])
        rows.append(row)
    return pd.DataFrame(rows)


# --- the step ------------------------------------------------------------------------------------------
def main(data: Path, composition: Path, tree_file: Path, arch_dir: Path, audit: Path, out: Path) -> pd.DataFrame:
    """Writes arch_cells.tsv (every cell with its raw and whole-run corrected p), arch_verdict.tsv and arch_measure_summary.tsv."""
    require_control(data, out)
    arch = pd.read_csv(arch_dir / "architecture.tsv", sep="\t")
    sp = prepare(arch, annotation_flags(audit), test_species(data, composition, tree_file))
    cells = run_cells(sp, tree_file)
    out.mkdir(parents=True, exist_ok=True)
    cells.to_csv(out / "arch_cells.tsv", sep="\t", index=False, float_format="%.6g")
    v = judge(cells)
    v.to_csv(out / "arch_verdict.tsv", sep="\t", index=False, float_format="%.6g")
    summarise(v).to_csv(out / "arch_measure_summary.tsv", sep="\t", index=False)
    print(f"{int(cells.p.notna().sum())} cells in the family; verdicts: " + ", ".join(f"{k} {n}" for k, n in v.verdict.value_counts().sort_index().items()))
    return v
