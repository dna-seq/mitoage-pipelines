"""Break each guard of H12's code in turn and check that a test fails.

Run by the H12 Slurm jobs after the tests have passed on the unbroken code: `cpu` (the default) for
`longevity structure-pgls`, `mitoport partners` and `structure embed`'s helpers; `embed` for the embedding
distance, which needs the table `structure embed` writes on a GPU. Each mutation replaces one line of a source
file, runs tests/test_structure_pgls.py and tests/test_partners.py, and puts the file back byte for byte (from a
copy: git is not installed on the compute nodes). Prints one line per mutation, "caught" or "SURVIVED", and exits
non-zero if any survived or a file did not come back.
"""
import hashlib
import subprocess
import sys
from pathlib import Path

FILE = Path("src/mitoage_pipelines/longevity/structure_pgls.py")
PARTNERS = Path("src/mitoage_pipelines/mitoport/partners.py")
EMBED = Path("src/mitoage_pipelines/structure/embed.py")
TESTS = ["tests/test_structure_pgls.py", "tests/test_partners.py"]
MUTATIONS = [
    ("genes not weighted by their columns",
     "num += share * len(cols[grp][g]); den += len(cols[grp][g])", "num += share; den += 1"),
    ("contact plus control instead of minus",
     'res["contact minus control"] = res.contact - res.control', 'res["contact minus control"] = res.contact + res.control'),
    ("retired records kept",
     "pgf = pgf[pgf.slug.isin(live) & pgf.region.isin(GRAVY_REGIONS)]", "pgf = pgf[pgf.region.isin(GRAVY_REGIONS)]"),
    ("another region read as total_pc",
     'return wide[[r for r in GRAVY_REGIONS if r in wide.columns]].add_prefix("gravy:")',
     'return wide[[r for r in GRAVY_REGIONS if r in wide.columns]].rename(columns={"total_pc": "x", "ND1": "total_pc"})'
     '.rename(columns={"x": "ND1"}).add_prefix("gravy:")'),
    ("feature not standardised", "z = (x[ok] - x[ok].mean()) / x[ok].std()", "z = x[ok]"),
    ("normal instead of t interval", 'h = stats.t.ppf(0.975, r["df"]) * r["se"]', 'h = 1.96 * r["se"]'),
    ("fitted on fewer than MIN_SPECIES values",
     "if ok.sum() < MIN_SPECIES or np.std(x[ok]) == 0:", "if ok.sum() < 3 or np.std(x[ok]) == 0:"),
    ("q over the whole class, not within family and model",
     'for _, idx in res[tested].groupby(["class", "family", "model"]).groups.items():',
     'for _, idx in res[tested].groupby(["class"]).groups.items():'),
    ("a class too small is dropped",
     'return [{"class": cls, "family": fam, "feature": "*", "model": m, "n": len(sp),',
     'return []; [{"class": cls, "family": fam, "feature": "*", "model": m, "n": len(sp),'),
    ("one model reports the class size as its n",
     'rows.append({**base, "model": model, "n": int(ok.sum()), "status": "tested",',
     'rows.append({**base, "model": model, "n": int(ok.sum()) if model != "mass" else len(sp), "status": "tested",'),
    ("a missing contact map not refused", '        if not (work / "sites.tsv").is_file():', "        if False:"),
]


# the second round of H12: (file, name, line, replacement)
MORE = [
    (FILE, "own composition read from the 12-gene covariate",
     "r = fit(C[np.ix_(ok, ok)], sp.ll.to_numpy()[ok], np.column_stack([z, sp.lm.to_numpy()[ok], cov[ok]]))",
     'r = fit(C[np.ix_(ok, ok)], sp.ll.to_numpy()[ok], np.column_stack([z, sp.lm.to_numpy()[ok], sp.loc[ok, ["A4", "C4", "T4"]].to_numpy()]))'),
    (FILE, "U share taken from G", 'out[r] = pd.DataFrame({"a": d.ff_a / n, "c": d.ff_c / n, "u": d.ff_u / n}).set_axis(d.slug)',
     'out[r] = pd.DataFrame({"a": d.ff_a / n, "c": d.ff_c / n, "u": d.ff_g / n}).set_axis(d.slug)'),
    (PARTNERS, "either side off its consensus instead of exactly one",
     "n += 1; d += (a != cons[(r.mt_gene, c)]) != (b != cons[(r.nuc_gene, r.nuc_col)])",
     "n += 1; d += (a != cons[(r.mt_gene, c)]) or (b != cons[(r.nuc_gene, r.nuc_col)])"),
    (PARTNERS, "consensus over every species, not the class",
     "cons[key] = consensus(aln, slugs, col)", "cons[key] = consensus(aln, list(aln), col)"),
    (PARTNERS, "species filling few pairs kept",
     'res.loc[res.contact_pairs < MIN_SHARE * len(pr), ["contact", "control"]] = np.nan', "pass"),
    (PARTNERS, "a pairing that disagrees with sites.tsv accepted",
     "        if (not set(m.contact_aln_col) <= contacts or m.contact_aln_col.duplicated().any()",
     "        if False and (not set(m.contact_aln_col) <= contacts or m.contact_aln_col.duplicated().any()"),
    (PARTNERS, "control pair is the contact itself",
     '"ctl_col": ctl - 1,', '"ctl_col": int(s.aln_col) - 1,'),
    (EMBED, "residue index off by one", "return [int(pos[c]) for c in cols if row[c] != \"-\"]",
     "return [int(pos[c]) + 1 for c in cols if row[c] != \"-\"]"),
]
EMBED_MUTATIONS = [
    (FILE, "genes not weighted in the embedding distance",
     "num = num.add(pd.Series(dist * weights[g], index=d.index), fill_value=np.nan)",
     "num = num.add(pd.Series(dist * 100, index=d.index), fill_value=np.nan)"),
    (FILE, "distance to the class mean taken as a dot product",
     "dist = 1 - X @ c / (np.linalg.norm(X, axis=1) * np.linalg.norm(c))", "dist = 1 - X @ c"),
]


def run(mutations) -> int:
    originals = {f: f.read_bytes() for f in {m[0] for m in mutations}}
    digests = {f: hashlib.sha256(b).hexdigest() for f, b in originals.items()}
    survived = 0
    try:
        for path, name, old, new in mutations:
            text = originals[path].decode()
            if old not in text:                 # a line that recurs further down: its first occurrence is the one meant
                print(f"NOT APPLIED (line not found): {name}"); survived += 1; continue
            path.write_text(text.replace(old, new, 1))
            r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "--show-capture=no", *TESTS],
                               capture_output=True, text=True)
            failed = [l.split(" - ")[0].removeprefix("FAILED ") for l in r.stdout.splitlines() if l.startswith(("FAILED", "ERROR"))]
            if r.returncode == 0:
                print(f"SURVIVED: {name}"); survived += 1
            else:
                print(f"caught: {name}: {', '.join(failed) or r.stdout.strip().splitlines()[-1]}")
            path.write_bytes(originals[path])
    finally:
        for f, b in originals.items():
            f.write_bytes(b)
    back = all(hashlib.sha256(f.read_bytes()).hexdigest() == d for f, d in digests.items())
    print(f"{len(mutations) - survived} of {len(mutations)} mutations caught; sources restored: {back}")
    return 0 if survived == 0 and back else 1


def main() -> int:
    group = sys.argv[1] if len(sys.argv) > 1 else "cpu"
    if group == "embed":
        return run(EMBED_MUTATIONS)
    return run([(FILE, *m) for m in MUTATIONS] + MORE)


if __name__ == "__main__":
    sys.exit(main())
