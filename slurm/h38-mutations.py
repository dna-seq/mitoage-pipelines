"""Break each guard of H38 (genome architecture against lifespan: `longevity/arch_measures.py`, `arch_tests.py`, the command that
exits non-zero) in turn and check that a test of `tests/test_arch.py` fails. Each mutation replaces one piece of one source file
exactly once, runs the tests, and the file is put back. Prints one line per mutation and exits non-zero if any survived, was not
applied, or if the unbroken tests do not all pass with none skipped (a skipped test breaks nothing).

The tests read the build (`MITOAGE_PIPELINES_DATA` with `audit records|strand|order|code`, `longevity arch-measures`, `controls` and
`dated-tree` made), the cybrid records and Tandem Repeats Finder on PATH; run from the root of a clone:

    uv run --frozen --group dev python slurm/h38-mutations.py [REGEX]

With a regular expression, only the mutations whose description matches it are run (the unbroken tests are still checked). With
`--count` the pieces are only looked for, each once, and nothing is run.
"""
import re
import subprocess
import sys
from pathlib import Path

SRC = Path("src/mitoage_pipelines/longevity")
M, T, C = SRC / "arch_measures.py", SRC / "arch_tests.py", Path("src/mitoage_pipelines/cli.py")
TESTS = "tests/test_arch.py"

H = '"p_holm_run"'
MUTATIONS = [
    # --- arch_measures.py: the distance
    ("an adjacency read from one strand only", M, "out[min((a, b), (_flip(b), _flip(a)))] += 1", "out[(a, b)] += 1"),
    ("the last gene not next to the first", M, "a, b = toks[i], toks[(i + 1) % n]", "a, b = toks[i], toks[min(i + 1, n - 1)]"),
    ("the distance not halved", M, "for k in set(ca) | set(cb)) / 2", "for k in set(ca) | set(cb))"),
    ("only the adjacencies one record has and the other lacks", M, "sum(abs(ca[k] - cb[k]) for k in set(ca) | set(cb)) / 2", "sum(max(ca[k] - cb[k], 0) for k in set(ca) | set(cb))"),
    ("an unknown arrangement at distance zero", M, "        return None\n    ca, cb = adjacencies(a), adjacencies(b)", "        return 0.0\n    ca, cb = adjacencies(a), adjacencies(b)"),
    ("a missing arrangement of the float kind not unknown", M, " or (isinstance(a, float) and np.isnan(a))", ""),
    ("the modal arrangement the last in string order on a tie", M, "return min(c.index[c.values == c.max()])", "return max(c.index[c.values == c.max()])"),
    ("the modal arrangement the rarest", M, "return min(c.index[c.values == c.max()])", "return min(c.index[c.values == c.min()])"),
    ("a record with no class at distance zero", M, "if isinstance(b, str) else None for a, b in zip", "if isinstance(b, str) else 0.0 for a, b in zip"),
    ("an unknown distance filled with zero", M, 'pd.array(dist, dtype="Float64")', 'pd.array([0.0 if x is None else x for x in dist], dtype="Float64")'),
    ("the class consensus from the arrangement as deposited", M, 'col = f"order_{tier}_unoriented"\n    modal = {}', 'col = f"order_{tier}"\n    modal = {}'),
    # --- the control region and the overlaps
    ("the control region's length from its first piece's feature only", M, "return len(cr), sum(s.end - s.start + 1 for f in cr for s in f.segments)", "return len(cr), sum(s.end - s.start + 1 for s in cr[0].segments)"),
    ("a control region that wraps the origin measured by its ends", M, "return len(cr), sum(s.end - s.start + 1 for f in cr for s in f.segments)", "return len(cr), sum(f.end - f.start + 1 for f in cr)"),
    ("a missing control region of length zero", M, "        return 0, None\n    return len(cr)", "        return 0, 0\n    return len(cr)"),
    ("a missing control region counted as zero copies", M, '"cr_copies": copies if copies else None', '"cr_copies": copies'),
    ("the control region's length not kept as unknown", M, '"cr_length": length,', '"cr_length": length or 0,'),
    ("a position in one gene an overlap", M, "(depth[1:n + 1] >= 2).sum()", "(depth[1:n + 1] >= 1).sum()"),
    ("a feature repeated in the same place an overlap", M, "        if key in seen:\n            continue", "        if False:\n            continue"),
    ("any feature a gene", M, "        if f.type not in CODING:\n            continue", "        if False:\n            continue"),
    # --- the repeats
    ("an unpinned binary accepted", M, "    if h != TRF_SHA256:", "    if False:"),
    ("the parameters of the run not checked", M, '            if line[len("Parameters: "):].strip() != want:', "            if False:"),
    ("a sequence TRF did not report on read as having none", M, "    if missing:\n        raise RuntimeError(f\"TRF gave no block", "    if False:\n        raise RuntimeError(f\"TRF gave no block"),
    ("a repeat's start and end swapped", M, '"start": int(f[0]), "end": int(f[1])', '"start": int(f[1]), "end": int(f[0])'),
    ("the parameters of the run changed", M, "TRF_PARAMS = (2, 7, 7, 80, 10, 50, 500)", "TRF_PARAMS = (2, 7, 7, 80, 10, 30, 500)"),
    ("overlapping repeats counted twice", M, "    whole = float(m.sum()) / n if n else float(\"nan\")", '    whole = float(sum(r["end"] - r["start"] + 1 for r in repeats)) / n if n else float("nan")'),
    ("the control region's repeat share over the whole genome", M, "float((m & cr).sum()) / float(cr.sum())", "float((m & cr).sum()) / float(n)"),
    ("the control region's repeat share over repeats outside it", M, "float((m & cr).sum()) / float(cr.sum())", "float(m.sum()) / float(cr.sum())"),
    ("no control region a zero repeat share", M, "if cr is None or not cr.any():\n        return whole, None", "if cr is None or not cr.any():\n        return whole, 0.0"),
    ("a repeat's last base left out", M, 'mask[max(r["start"], 1):min(r["end"], n) + 1] = True', 'mask[max(r["start"], 1):min(r["end"], n)] = True'),
    ("retired species measured and kept", M, '[species.record_status.ne("retired")].copy()', ".copy()"),
    ("the table not in slug order", M, 'return tab.sort_values("slug", kind="mergesort").reset_index(drop=True), rep.sort_values', "return tab.reset_index(drop=True), rep.sort_values"),
    ("the stored non-coding share over the wrong size", M, 'tab["non_coding_frac"] = tab.non_coding_bp / tab.mtdna_size_bp', 'tab["non_coding_frac"] = tab.non_coding_bp / tab.genome_length'),
    # --- arch_tests.py: the annotation flag
    ("any of the four reasons not enough, all four needed", T, ".any(axis=1).astype(int)", ".all(axis=1).astype(int)"),
    ("a mixed strand not flagged", T, 'strand.orientation.isin(["reversed", "mixed"])', 'strand.orientation.isin(["reversed"])'),
    ("a record with a clean code flagged", T, 'code.neither_clean == "yes"', 'code.neither_clean == "no"'),
    ("a record with no current features in the build flagged", T, 'build.build_annotation == "differs"', 'build.build_annotation != "same"'),
    ("a species the tables do not share tolerated", T, "        if set(d.slug) != slugs:", "        if False:"),
    ("the order flag of one tier only", T, 'out["flag_order"] = out.slug.isin(set(flags.slug))', 'out["flag_order"] = out.slug.isin(set(flags.slug[flags.tier == "full"]))'),
    # --- the positive control
    ("the arrangements not compared with the audit's", T, '        exceptions += int(((a != b) & ~(a.isna() & b.isna())).sum())', "        exceptions += 0"),
    ("a species the audit has no row for tolerated", T, 'exceptions = int((merged._merge != "both").sum())', "exceptions = 0"),
    ("a flagged record at zero from its family passes", T, '"exceptions": int((~(d1 > 0)).sum())', '"exceptions": 0'),
    ("a flagged record at zero from its family counted above zero", T, "(~(d1 > 0)).sum()", "(~(d1 >= 0)).sum()"),
    ("a reversed deposit above zero passes", T, '"exceptions": int((~(d2 == 0)).sum())', '"exceptions": 0'),
    ("a reversed deposit checked as a rearrangement", T, 'rev, odd = fl[fl.kind == "reversed deposit"], fl[fl.kind != "reversed deposit"]', 'rev, odd = fl[fl.kind == "reversed deposit"], fl'),
    ("a uniform family above zero passes", T, '"exceptions": int((fam_d.reindex(uni_slugs) != 0).sum())', '"exceptions": 0'),
    ("a family of one is uniform", T, "if len(d) >= aorder.MIN_FAMILY and d[col].nunique() == 1:", "if len(d) >= 1 and d[col].nunique() == 1:"),
    ("a family that differs only by the strand is uniform", T, "if len(d) >= aorder.MIN_FAMILY and d[col].nunique() == 1:", 'if len(d) >= aorder.MIN_FAMILY and d[f"{col}_unoriented"].nunique() == 1:'),
    ("a class distance unlike the family's passes", T, "g4_bad += int(float(fam_d.at[s]) != float(cls.at[s, f\"rearr_{tier}\"]))", "g4_bad += 0"),
    ("the family's and the class's distance compared where they differ in reference", T, 'if famrow_modal == cls.at[s, f"modal_{tier}"]:', "if True:"),
    ("a failed control reproduced", T, '"reproduced": bool(len(bad) == 0)', '"reproduced": True'),
    ("a control with one exception reproduced", T, "bad = tab[tab.exceptions > 0]", "bad = tab[tab.exceptions > 1]"),
    ("a control of another build accepted", T, '    if v.get("species_sha256") != am.sha256_of(data / "species.parquet"):', "    if False:"),
    ("a failed control accepted", T, '    if v.get("reproduced") is not True:', "    if False:"),
    ("no control file needed", T, "    if not p.is_file():\n        raise SystemExit(f\"no {p}", "    if False:\n        raise SystemExit(f\"no {p}"),
    ("the gate before the tests removed", T, "    require_control(data, out)\n    arch = pd.read_csv", "    arch = pd.read_csv"),
    ("a failed control exits as a passed one", C, "        raise typer.Exit(1)          # the tables and the json are written; a job that depends on this one does not start", "        pass"),
    # --- the species and the regressions
    ("a species with no architecture row tolerated", T, "    if missing:\n        raise ValueError(f\"{len(missing)} test species have no architecture row", "    if False:\n        raise ValueError(f\"{len(missing)} test species have no architecture row"),
    ("a species with no audit row taken as unflagged", T, "    if unaudited:", "    if False:"),
    ("lengths not in log10", T, "a[c] = np.log10(a[c].where(a[c] > 0))", "a[c] = a[c]"),
    ("a length of zero a number", T, "a[c] = np.log10(a[c].where(a[c] > 0))", "a[c] = np.log10(a[c].where(a[c] >= 0).fillna(1))"),
    ("the high set the confidence-filtered one", T, 'd["high"] = species_sets.members(d, "high")', 'd["high"] = species_sets.members(d, "not_unreliable")'),
    ("the primary set not confidence-filtered", T, 'd["not_unreliable"] = species_sets.members(d, "not_unreliable")', 'd["not_unreliable"] = True'),
    ("the measure not standardised", T, "    z = (y - y.mean()) / y.std(ddof=1)\n    cols", "    z = y\n    cols"),
    ("body mass left out of the model", T, "cols = [ll, lm] if flag is None else [ll, lm, flag.astype(float)]", "cols = [ll] if flag is None else [ll, flag.astype(float)]"),
    ("the annotation flag left out of its model", T, "else [ll, lm, flag.astype(float)]", "else [ll, lm]"),
    ("the flag tested instead of lifespan", T, "else [ll, lm, flag.astype(float)]", "else [flag.astype(float), lm, ll]"),
    ("lambda 1 not applied", T, '("lambda_1", 1.0)', '("lambda_1", None)'),
    ("the raw regression estimating lambda", T, '("raw", 0.0)', '("raw", None)'),
    ("the confidence filter not applied to the species of a cell", T, 'cell = sp[(sp["class"] == cls) & sp[name]]', 'cell = sp[(sp["class"] == cls)]'),
    ("a measure of few species tested", T, "                if n < MIN_SPECIES:\n                    reason", "                if n < 3:\n                    reason"),
    ("a nearly constant measure tested", T, "                elif n - common < MIN_MINORITY:", "                elif False:"),
    ("the minimum off the mode not enough", T, "                elif n - common < MIN_MINORITY:", "                elif n - common <= MIN_MINORITY:"),
    ("a small set tested", T, "            if len(cell) < MIN_SPECIES:", "            if len(cell) < 3:"),
    ("a flag no species or all species have entered", T, "estimable = 0 < flag.sum() < n", "estimable = True"),
    ("the model with the flag fitted without it", T, "fit_one(Cg, y, ll, lm, flag if cov else None, lam)", "fit_one(Cg, y, ll, lm, None, lam)"),
    # --- the correction over the whole run
    ("no correction in the table of the p-values", T, "    return correct_run(pd.DataFrame(rows))", "    return pd.DataFrame(rows)"),
    ("Holm replaced by Benjamini-Hochberg", T, 'cells.loc[have, "p_holm_run"] = multipletests(cells.p[have], method="holm")[1]', 'cells.loc[have, "p_holm_run"] = multipletests(cells.p[have], method="fdr_bh")[1]'),
    ("q as Holm", T, 'cells.loc[have, "q_run"] = multipletests(cells.p[have], method="fdr_bh")[1]', 'cells.loc[have, "q_run"] = multipletests(cells.p[have], method="holm")[1]'),
    ("Holm within each measure", T, 'cells.loc[have, "p_holm_run"] = multipletests(cells.p[have], method="holm")[1]', 'for _, g in cells[have].groupby("measure"):\n        cells.loc[g.index, "p_holm_run"] = multipletests(g.p, method="holm")[1]'),
    ("Holm within each taxon group", T, 'cells.loc[have, "p_holm_run"] = multipletests(cells.p[have], method="holm")[1]', 'for _, g in cells[have].groupby("class"):\n        cells.loc[g.index, "p_holm_run"] = multipletests(g.p, method="holm")[1]'),
    ("Holm within each species set", T, 'cells.loc[have, "p_holm_run"] = multipletests(cells.p[have], method="holm")[1]', 'for _, g in cells[have].groupby("set"):\n        cells.loc[g.index, "p_holm_run"] = multipletests(g.p, method="holm")[1]'),
    ("Holm within each model", T, 'cells.loc[have, "p_holm_run"] = multipletests(cells.p[have], method="holm")[1]', 'for _, g in cells[have].groupby("model"):\n        cells.loc[g.index, "p_holm_run"] = multipletests(g.p, method="holm")[1]'),
    ("Holm within each variant of the covariate", T, 'cells.loc[have, "p_holm_run"] = multipletests(cells.p[have], method="holm")[1]', 'for _, g in cells[have].groupby("covariate"):\n        cells.loc[g.index, "p_holm_run"] = multipletests(g.p, method="holm")[1]'),
    ("the raw regression left out of the family", T, "    have = cells.p.notna()\n    cells.loc[have, \"p_holm_run\"]", "    have = cells.p.notna() & (cells.model != \"raw\")\n    cells.loc[have, \"p_holm_run\"]"),
    ("the rows with the flag left out of the family", T, "    have = cells.p.notna()\n    cells.loc[have, \"p_holm_run\"]", "    have = cells.p.notna() & ~cells.covariate.astype(bool)\n    cells.loc[have, \"p_holm_run\"]"),
    ("the high set left out of the family", T, "    have = cells.p.notna()\n    cells.loc[have, \"p_holm_run\"]", "    have = cells.p.notna() & (cells.set != \"high\")\n    cells.loc[have, \"p_holm_run\"]"),
    # --- the bar
    ("a table whose correction was not computed with the p-values judged", T, '    for c in ("p_holm_run", "q_run"):\n        if c not in cells.columns:\n            raise', '    for c in ():\n        if c not in cells.columns:\n            raise'),
    ("lambda 1 not required for the controls", T, "c1 = bool(ml.p < ALPHA and l1.p < ALPHA and sgn == np.sign(l1.slope))", "c1 = bool(ml.p < ALPHA and sgn == np.sign(l1.slope))"),
    ("lambda estimated not required for the controls", T, "c1 = bool(ml.p < ALPHA and l1.p < ALPHA and sgn == np.sign(l1.slope))", "c1 = bool(l1.p < ALPHA and sgn == np.sign(l1.slope))"),
    ("the sign across lambda forms not required", T, "c1 = bool(ml.p < ALPHA and l1.p < ALPHA and sgn == np.sign(l1.slope))", "c1 = bool(ml.p < ALPHA and l1.p < ALPHA)"),
    ("the threshold of the controls inclusive", T, "c1 = bool(ml.p < ALPHA and l1.p < ALPHA", "c1 = bool(ml.p <= ALPHA and l1.p <= ALPHA"),
    ("the raw regression decides", T, "c1 = bool(ml.p < ALPHA and", "c1 = bool(raw is not None and raw.p < ALPHA and ml.p < ALPHA and"),
    ("the threshold of 5 per cent changed to 10", T, "ALPHA = 0.05", "ALPHA = 0.10"),
    ("only lambda estimated on the high subset", T, 'for k in ("lambda_ml", "lambda_1")))\n        c3_prim', 'for k in ("lambda_ml",)))\n        c3_prim'),
    ("only lambda 1 on the high subset", T, 'for k in ("lambda_ml", "lambda_1")))\n        c3_prim', 'for k in ("lambda_1",)))\n        c3_prim'),
    ("the sign on the high subset not required", T, ' and np.sign(high.at[k, "slope"]) == sgn for k in', " for k in"),
    ("the threshold on the high subset inclusive", T, 'all(high.at[k, "p"] < ALPHA and', 'all(high.at[k, "p"] <= ALPHA and'),
    ("lambda 1 not required to survive the correction", T, "c3_prim = bool(ml.p_holm_run < ALPHA and l1.p_holm_run < ALPHA)", "c3_prim = bool(ml.p_holm_run < ALPHA)"),
    ("lambda estimated not required to survive the correction", T, "c3_prim = bool(ml.p_holm_run < ALPHA and l1.p_holm_run < ALPHA)", "c3_prim = bool(l1.p_holm_run < ALPHA)"),
    ("the threshold of the correction inclusive", T, "c3_prim = bool(ml.p_holm_run < ALPHA and l1.p_holm_run < ALPHA)", "c3_prim = bool(ml.p_holm_run <= ALPHA and l1.p_holm_run <= ALPHA)"),
    ("the correction not required on the high subset", T, 'bool(c3_prim and all(high.at[k, "p_holm_run"] < ALPHA for k in ("lambda_ml", "lambda_1")))', "bool(c3_prim)"),
    ("only lambda estimated corrected on the high subset", T, 'all(high.at[k, "p_holm_run"] < ALPHA for k in ("lambda_ml", "lambda_1"))', 'all(high.at[k, "p_holm_run"] < ALPHA for k in ("lambda_ml",))'),
    ("only lambda 1 corrected on the high subset", T, 'all(high.at[k, "p_holm_run"] < ALPHA for k in ("lambda_ml", "lambda_1"))', 'all(high.at[k, "p_holm_run"] < ALPHA for k in ("lambda_1",))'),
    ("the correction over the run not required", T, 'elif c2 is False or c3 is False:', 'elif c2 is False:'),
    ("the high subset not required", T, 'elif c2 is False or c3 is False:', 'elif c3 is False:'),
    ("a missing high subset passed", T, "        elif c2 is None or c3 is None:", "        elif False:"),
    ("a missing flag model passed", T, "        elif c4 is None:", "        elif False:"),
    ("a result the flag moves counts", T, "        elif not c4:", "        elif False:"),
    ("a result the flag moves reported as a finding", T, 'row["verdict"] = "explained by annotation"', 'row["verdict"] = "counts"'),
    ("only lambda estimated checked against the flag", T, 'for k in ("lambda_ml", "lambda_1"):\n                s0, s1', 'for k in ("lambda_ml",):\n                s0, s1'),
    ("only lambda 1 checked against the flag", T, 'for k in ("lambda_ml", "lambda_1"):\n                s0, s1', 'for k in ("lambda_1",):\n                s0, s1'),
    ("a flag that moves one lambda form is not a move", T, "c4 = not any(moved)", "c4 = not all(moved)"),
    ("a change of sign by the flag allowed", T, "not (np.sign(s1) == np.sign(s0) and cov.at", "not (cov.at"),
    ("the corrected threshold not required with the flag", T, 'and cov.at[k, "p_holm_run"] < ALPHA and abs(s1 - s0)', "and abs(s1 - s0)"),
    ("any change of the slope allowed", T, " and abs(s1 - s0) <= MOVE_TOLERANCE * abs(s0)", ""),
    ("the tolerance exclusive", T, "abs(s1 - s0) <= MOVE_TOLERANCE * abs(s0)", "abs(s1 - s0) < MOVE_TOLERANCE * abs(s0)"),
    ("the tolerance a half", T, "MOVE_TOLERANCE = 0.25", "MOVE_TOLERANCE = 0.5"),
    ("a failed first condition not a null", T, '        if not c1:\n            row["verdict"] = "null"', '        if False:\n            row["verdict"] = "null"'),
    ("a measure counting where all its groups count", T, 'row["measure_counts"] = bool((g.verdict == "counts").any())', 'row["measure_counts"] = bool((g.verdict == "counts").all())'),
]


def run(tests: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", "-rs", tests], capture_output=True, text=True)


originals = {p: p.read_text() for p in {m[1] for m in MUTATIONS}}
if "--count" in sys.argv:
    bad = [(m[0], originals[m[1]].count(m[2])) for m in MUTATIONS if originals[m[1]].count(m[2]) != 1]
    print(f"{len(MUTATIONS)} mutations; not found exactly once: {bad}")
    sys.exit(1 if bad else 0)
try:
    base = run(TESTS)
    last = base.stdout.strip().splitlines()[-1] if base.stdout.strip() else ""
    if base.returncode != 0 or "skipped" in last:
        print(f"the unbroken tests of {TESTS} are not all passing with none skipped:\n{base.stdout[-1500:]}")
        sys.exit(2)
    print(f"baseline {TESTS}: {last}")
    caught = 0
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    chosen = [m for m in MUTATIONS if re.search(args[0], m[0])] if args else MUTATIONS
    for what, f, old, new in chosen:
        text = originals[f]
        if text.count(old) != 1:
            raise SystemExit(f"{what}: the piece to break is in {f} {text.count(old)} times, not once")
        f.write_text(text.replace(old, new))
        r = run(TESTS)
        failed = r.returncode != 0
        caught += failed
        print(("caught: " if failed else "SURVIVED: ") + what + f"  [{f.name}]", flush=True)
        f.write_text(text)
finally:
    for p, t in originals.items():
        p.write_text(t)
print(f"{caught} of {len(chosen)} mutations caught; sources restored: {all(p.read_text() == t for p, t in originals.items())}")
sys.exit(0 if caught == len(chosen) else 1)
