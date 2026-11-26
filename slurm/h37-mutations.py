"""Break each guard of H37 part A (the RNA genes against lifespan: `longevity/rna_genes.py`, `rna_rates.py`, `rna_tests.py`) in
turn and check that a test of its test file fails. Each mutation replaces one piece of one source file exactly once, runs
the tests of that file, and the file is put back. Prints one line per mutation and exits non-zero if any survived, was
not applied, or if the unbroken tests do not all pass with none skipped (a skipped test breaks nothing).

The tests read the build (`MITOAGE_PIPELINES_DATA` with `longevity controls`, `fourfold`, `dated-tree` made), the cybrid
records and MAFFT on PATH; run from the root of a clone:

    uv run --frozen --group dev python slurm/h37-mutations.py [REGEX]

With a regular expression, only the mutations whose description matches it are run (the unbroken tests are still checked).
"""
import re
import subprocess
import sys
from pathlib import Path

SRC = Path("src/mitoage_pipelines/longevity")
G, R, T = SRC / "rna_genes.py", SRC / "rna_rates.py", SRC / "rna_tests.py"
TG, TR, TT = "tests/test_rna_genes.py", "tests/test_rna_rates.py", "tests/test_rna_tests.py"

MUTATIONS = [
    # --- rna_genes.py
    ("a duplicated gene takes its last copy", G, "        if f.name not in seqs:", "        if True:", TG),
    ("an unresolved Leu or Ser copy is not set apart", G, "        if f.name in PLAIN_ISO:", "        if False:", TG),
    ("a gene is read in the record's strand, not its own sense", G, "seqs[f.name] = gm.extract(g.record.sequence, f.segments).upper()", "seqs[f.name] = g.record.sequence[f.start - 1:f.end].upper()", TG),
    ("the rRNA features are not read", G, 'if f.type not in ("tRNA", "rRNA"):', 'if f.type not in ("tRNA",):', TG),
    ("the gene order is not used to place Leu and Ser", G, 'gm.fallback_table(phylum), phylum == "Chordata")', "gm.fallback_table(phylum), False)", TG),
    ("an unresolved gene counted as not annotated", G, '"n_not_annotated": len(lacking) - len(unres),', '"n_not_annotated": len(lacking),', TG),
    ("any second feature counted as a duplicate", G, 'r["copies"].get(gene, 0) > 1', 'r["copies"].get(gene, 0) > 0', TG),
    ("no species ever counted as unresolved", G, 'unres = [s for s in lacking if aa and records[s]["unresolved"].get(aa, 0) > 0]', "unres = []", TG),
    ("a column with a gap in most rows is occupied", G, "if sum(counts.values()) < INFORMATIVE_OCCUPANCY * len(rows):", "if False:", TG),
    ("a base seen once makes a column informative", G, "if sum(1 for v in counts.values() if v >= 2) >= 2:", "if sum(1 for v in counts.values() if v >= 1) >= 2:", TG),
    ("gaps counted in a gene's length and bases", G, 'seq = seq.upper().replace("-", "")', "seq = seq.upper()", TG),
    ("shares of an unresolved sequence read as zero", G, "(b[x.lower()] / n if n else None)", "(b[x.lower()] / max(n, 1))", TG),
    ("an unpinned release read", G, "    if md5_of(release) != pins.REFSEQ_MD5:", "    if False:", TG),
    ("a missing record passed over", G, "    if missing:", "    if False:", TG),
    ("a retired species loaded", G, '    sp = sp[sp.record_status.ne("retired")]\n', "", TG),
    ("a CONTIG-only record read without its download", G, "            if not p.is_file():", "            if False:", TG),
    ("species not in slug order", G, 'for r in sp.sort_values("slug").itertuples(index=False):', "for r in sp.itertuples(index=False):", TG),
    ("rows aligned in the order given", G, 'for slug, seq in sorted(seqs.items())]', "for slug, seq in seqs.items()]", TG),
    ("two sequences aligned", G, "    if len(rows) < MIN_ALIGNED:", "    if len(rows) < 1:", TG),
    # --- rna_rates.py
    ("no floor of shared columns", R, "    if len(pairs) < MIN_SITES:", "    if len(pairs) < 0:", TR),
    ("no Jukes-Cantor correction", R, '"dist": rates.jukes_cantor(p), "sites"', '"dist": p, "sites"', TR),
    ("the two sisters measured on their own columns", R, "    keep = [i for i in range(len(s1)) if all(r[i] in BASES for r in rows)]", "    keep = [i for i in range(len(s1)) if s1[i] in BASES]", TR),
    ("rows of unequal length accepted", R, "    if len({len(s1), len(s2), len(outgroup)}) != 1:", "    if False:", TR),
    ("saturation at distance two", R, '> rates.SATURATED_DS\n    return c\n\n\n# --- the tests', "> 2\n    return c\n\n\n# --- the tests", TR),
    ("every gene on the species of the whole frame", R, "        usable = set(ali) & set(species.slug)", "        usable = set(species.slug)", TR),
    ("the cherry's depth is its longer branch", R, '"time_my": float(np.mean([k.branch_length or 0.0 for k in node.clades])), "sites"', '"time_my": float(max([k.branch_length or 0.0 for k in node.clades])), "sites"', TR),
    ("the contrast reversed", R, 'c["d_dist"] = c.dist_a - c.dist_b', 'c["d_dist"] = c.dist_b - c.dist_a', TR),
    ("the effect size is not a partial correlation", R, "return t / math.sqrt(t * t + dof)", "return t / math.sqrt(t * t + dof + 1)", TR),
    ("Benjamini-Hochberg where the bar names Holm", R, 'multipletests(p, method="holm")', 'multipletests(p, method="fdr_bh")', TR),
    ("q over saturated and unsaturated pairs together", R, 'grp = ["class", "normalisation", "includes_saturated"]', 'grp = ["class", "normalisation"]', TR),
    ("the control ignores differing cells", R, "bad = [c for c in cols if c in m and m[c][mk[k]] != ref[c][rk[k]]]", "bad = []", TR),
    ("the control compares at two digits", R, 'f"{v:.5g}"', 'f"{v:.2g}"', TR),
    ("a table with differing rows is reproduced", R, '"reproduced": not (only_mine or only_ref or differing or missing_cols)}', '"reproduced": not (only_mine or only_ref or missing_cols)}', TR),
    ("a comparison that always reproduces", R, '"reproduced": not (only_mine or only_ref or differing or missing_cols)}', '"reproduced": True}', TR),
    ("a table with an extra row is reproduced", R, '"reproduced": not (only_mine or only_ref or differing or missing_cols)}', '"reproduced": not (only_ref or differing or missing_cols)}', TR),
    ("a table with a missing row is reproduced", R, '"reproduced": not (only_mine or only_ref or differing or missing_cols)}', '"reproduced": not (only_mine or differing or missing_cols)}', TR),
    ("a table with a missing column is reproduced", R, '"reproduced": not (only_mine or only_ref or differing or missing_cols)}', '"reproduced": not (only_mine or only_ref or differing)}', TR),
    ("a control of another build accepted", R, '    if v.get("species_sha256") != _sha(data / "species.parquet"):', "    if False:", TR),
    ("a failed control accepted", R, '    if v.get("reproduced") is not True:', "    if False:", TR),
    ("no control file needed", R, "    if not p.is_file():\n        raise SystemExit(f\"no {p}", "    if False:\n        raise SystemExit(f\"no {p}", TR),
    ("the environment restriction accepted", R, "    if species_sets.active() is not None:", "    if False:", TR),
    ("the high set not applied", R, "        sp = sp[sp.slug.isin(species_sets.slugs(data, species_set))].reset_index(drop=True)", "        sp = sp.reset_index(drop=True)", TR),
    ("the cherries of the genes some species lack", R, "    usable = set.intersection(*[set(a) for a in alignments.values()]) & set(species.slug)", "    usable = set.union(*[set(a) for a in alignments.values()]) & set(species.slug)", TR),
    ("a failed control exits as a passed one", Path("src/mitoage_pipelines/cli.py"), "        raise typer.Exit(1)          # the table and its json are written", "        pass          # the table and its json are written", TR),
    # --- rna_tests.py
    ("the PIC model not required", T, "d.p_holm < ALPHA and d.p_holm_pic < ALPHA and", "d.p_holm < ALPHA and", TT),
    ("the sign across models not required", T, " and np.sign(d.slope) == np.sign(d.slope_pic))", ")", TT),
    ("the threshold inclusive", T, "d.p_holm < ALPHA and d.p_holm_pic", "d.p_holm <= ALPHA and d.p_holm_pic", TT),
    ("a strictly larger effect required", T, "c2 = bool(abs(d.r) >= median)", "c2 = bool(abs(d.r) > median)", TT),
    ("the effect's sign counted", T, "c2 = bool(abs(d.r) >= median)", "c2 = bool(d.r >= median)", TT),
    ("the high set not required", T, "c3 = bool(h is not None and pd.notna(h.p_holm) and h.p_holm < ALPHA and np.sign(h.slope) == np.sign(d.slope))", "c3 = True", TT),
    ("the high set's sign not required", T, " and h.p_holm < ALPHA and np.sign(h.slope) == np.sign(d.slope))", " and h.p_holm < ALPHA)", TT),
    ("the high set's threshold not required", T, " and h.p_holm < ALPHA and np.sign(h.slope)", " and np.sign(h.slope)", TT),
    ("few informative columns inclusive", T, "d.informative_columns < WEAK_COLUMNS", "d.informative_columns <= WEAK_COLUMNS", TT),
    ("an unknown column count read as informative", T, "weak = pd.isna(d.informative_columns) or d.informative_columns < WEAK_COLUMNS", "weak = d.informative_columns < WEAK_COLUMNS", TT),
    ("a median from one protein gene", T, "if len(prot) >= MIN_PROTEIN else np.nan", "if len(prot) >= 1 else np.nan", TT),
    ("the median over both sets", T, '& (protein.set == "default")].dropna(subset=["r"])', "].dropna(subset=[\"r\"])", TT),
    ("the median over every group", T, '& (protein["class"] == cls) & (protein.set', "& (protein.set", TT),
    ("the median over every test", T, "(protein.test == test) & ", "", TT),
    ("a gene with no row is left out of the verdicts", T, "            if gene not in have:", "            if False:", TT),
    ("Holm replaced by Benjamini-Hochberg", T, 'res.loc[g.index, "p_holm"] = multipletests(g.p, method="holm")[1]', 'res.loc[g.index, "p_holm"] = multipletests(g.p, method="fdr_bh")[1]', TT),
    ("the correction over the features together", T, 'rna.groupby(["set", "class", "feature", "model"])', 'rna.groupby(["set", "class", "model"])', TT),
    ("the protein genes corrected with the RNA genes", T, 'rna = res[(res.group == "rna") & res.p.notna()]', "rna = res[res.p.notna()]", TT),
    ("the feature not standardised", T, "z = (y - y.mean()) / y.std(ddof=1)", "z = y", TT),
    ("body mass left out of the model", T, "np.column_stack([ll, lm]), lam)", "np.column_stack([ll]), lam)", TT),
    ("lambda fixed at one not applied", T, '(("lambda_ml", None), ("lambda_1", 1.0))', '(("lambda_ml", None), ("lambda_1", None))', TT),
    ("the high set not applied to the species", T, '& (sp0.high if name == "high" else True)]', "& True]", TT),
    ("informative columns over the whole alignment", T, "informative_columns(ali[gene], set(f.slug))", "informative_columns(ali[gene])", TT),
    ("a gene of too few species tested", T, "if ok.sum() < MIN_SPECIES or", "if ok.sum() < 3 or", TT),
    ("a constant feature tested", T, "np.nanstd(y[ok]) == 0", "False", TT),
    ("lambda one result read as the primary", T, 'a = tests[tests.model == "lambda_ml"]', 'a = tests[tests.model == "lambda_1"]', TT),
    ("the rate with the body-mass contrast not chosen", T, '(rna.normalisation == "d_mass_adj")', '(rna.normalisation == "d_raw")', TT),
    ("saturated pairs in the rate", T, '& ~rna.includes_saturated].copy()', "].copy()", TT),
    ("the protein rate not dS", T, '(protein.normalisation == "dS_mass_adj")', '(protein.normalisation == "dN_raw")', TT),
    ("the summary's best cell by the smallest effect", T, "b = ok.loc[ratio.idxmax()]", "b = ok.loc[ratio.idxmin()]", TT),
    ("the summary counts a cell without a verdict", T, 'row["n_not_assessed"] = int(g.verdict.str.startswith("not assessed").sum())', 'row["n_not_assessed"] = 0', TT),
    ("the summary's survivors counted from the verdicts that count", T, 'row["n_survives_controls"] = int(g.survives_controls.fillna(False).astype(bool).sum())', 'row["n_survives_controls"] = int((g.verdict == "counts").sum())', TT),
    ("the effect over the median taken as the effect", T, "ratio = (ok.r.abs() / ok.median_protein_r)", "ratio = (ok.r.abs() / 1.0)", TT),
    ("the gate before the composition tests removed", T, "    require_control(data, out)\n    res = run_class", "    res = run_class", TT),
]


def run(tests: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", "-rs", tests], capture_output=True, text=True)


originals = {p: p.read_text() for p in {m[1] for m in MUTATIONS}}
try:
    for tests in (TG, TR, TT):
        base = run(tests)
        last = base.stdout.strip().splitlines()[-1] if base.stdout.strip() else ""
        if base.returncode != 0 or "skipped" in last:
            print(f"the unbroken tests of {tests} are not all passing with none skipped:\n{base.stdout[-1500:]}")
            sys.exit(2)
        print(f"baseline {tests}: {last}")
    caught = 0
    chosen = [m for m in MUTATIONS if re.search(sys.argv[1], m[0])] if len(sys.argv) > 1 else MUTATIONS
    for what, f, old, new, tests in chosen:
        text = originals[f]
        if text.count(old) != 1:
            raise SystemExit(f"{what}: the piece to break is in {f} {text.count(old)} times, not once")
        f.write_text(text.replace(old, new))
        r = run(tests)
        failed = r.returncode != 0
        caught += failed
        print(("caught: " if failed else "SURVIVED: ") + what + f"  [{f.name}]", flush=True)
        f.write_text(text)
finally:
    for p, t in originals.items():
        p.write_text(t)
print(f"{caught} of {len(chosen)} mutations caught; sources restored: {all(p.read_text() == t for p, t in originals.items())}")
sys.exit(0 if caught == len(chosen) else 1)
