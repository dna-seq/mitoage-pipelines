"""Break each guard of `longevity structure-pgls` in turn and check that a test of tests/test_structure_pgls.py fails.

Run by slurm/structure-pgls.sbatch after the tests have passed on the unbroken code. Each mutation replaces one line of
src/mitoage_pipelines/longevity/structure_pgls.py, runs the tests, and puts the file back byte for byte (from a
copy: git is not installed on the compute nodes). Prints one line per mutation, "caught" or "SURVIVED", and exits
non-zero if any survived or the file did not come back.
"""
import hashlib
import subprocess
import sys
from pathlib import Path

FILE = Path("src/mitoage_pipelines/longevity/structure_pgls.py")
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


def main() -> int:
    original = FILE.read_bytes()
    digest = hashlib.sha256(original).hexdigest()
    survived = 0
    try:
        for name, old, new in MUTATIONS:
            text = original.decode()
            if text.count(old) != 1:
                print(f"NOT APPLIED (line not found once): {name}"); survived += 1; continue
            FILE.write_text(text.replace(old, new))
            r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "--show-capture=no",
                                "tests/test_structure_pgls.py"], capture_output=True, text=True)
            failed = [l.split(" - ")[0].removeprefix("FAILED ") for l in r.stdout.splitlines() if l.startswith(("FAILED", "ERROR"))]
            if r.returncode == 0:
                print(f"SURVIVED: {name}"); survived += 1
            else:
                print(f"caught: {name}: {', '.join(failed) or r.stdout.strip().splitlines()[-1]}")
            FILE.write_bytes(original)
    finally:
        FILE.write_bytes(original)
    back = hashlib.sha256(FILE.read_bytes()).hexdigest() == digest
    print(f"{len(MUTATIONS) - survived} of {len(MUTATIONS)} mutations caught; source restored: {back}")
    return 0 if survived == 0 and back else 1


if __name__ == "__main__":
    sys.exit(main())
