"""Break each guard of H20's code (`longevity thermal`) in turn and check that a test fails.

Run by slurm/h20-thermal.sbatch after the tests have passed on the unbroken code. Each mutation replaces one line of
the source, runs tests/test_thermal.py (stopping at the first failure), and puts the file back byte for byte from a
copy (git is not installed on the compute nodes). Prints one line per mutation, "caught" or "SURVIVED", and exits
non-zero if any survived, was not applied, or the file did not come back.
"""
import hashlib
import signal
import subprocess
import sys
from pathlib import Path

TH = Path("src/mitoage_pipelines/longevity/thermal.py")
TESTS = ["tests/test_thermal.py"]
MUTATIONS = [
    (TH, "a body temperature out of range is kept",
     "body = species.body_temperature_k.where(species.body_temperature_k.between(*BODY_K))",
     "body = species.body_temperature_k"),
    (TH, "a water temperature out of range is kept", "ok = lambda s: s.where(s.between(*WATER_C))", "ok = lambda s: s"),
    (TH, "several stocks give their mean, not their median",
     'temps = st.groupby("SpecCode")[["tpref", "trange"]].median()', 'temps = st.groupby("SpecCode")[["tpref", "trange"]].mean()'),
    (TH, "1/kT is computed with the wrong constant", "    return 1.0 / (K_EV * kelvin)", "    return 1.0 / (8.617e-5 * 1.01 * kelvin)"),
    (TH, "fishes take their size from mass", 't["size"] = t.lm.where(~fish, t.llen)', 't["size"] = t.lm'),
    (TH, "a synonym of two species takes the first", "        if len(codes) > 1:", "        if False:"),
    (TH, "a misapplied name counts as a synonym", "    syn = syn[syn.Status.isin(SYNONYM_STATUS)]", "    syn = syn"),
    (TH, "a match that disagrees with the build's record is kept", "        if m and int(m.group(1)) != code:", "        if False:"),
    (TH, "a species goes into two pairs",
     "            (x, dx), (y, dy) = up.pop(0), up.pop(0)", "            (x, dx), (y, dy) = up.pop(0), up[0]"),
    (TH, "the members left over below a node are not passed up", "        return up\n", "        return []\n"),
    (TH, "a pair's age is the sum of its two branches", "            out.append((x, y, (dx + dy) / 2))", "            out.append((x, y, dx + dy))"),
    (TH, "the rate is per pair, not per branch", "    per_branch = 2 * p.age_my.where(p.age_my > 0)", "    per_branch = p.age_my.where(p.age_my > 0)"),
    (TH, "a pair takes the one sister's trait when the other lacks it",
     "**{c: (t.at[a, c] + t.at[b, c]) / 2 for c in terms}})", "**{c: np.nanmean([t.at[a, c], t.at[b, c]]) for c in terms}})"),
    (TH, "saturated pairs are used", '    p["usable"] = ~p.saturated & (p.dS > 0)', '    p["usable"] = (p.dS > 0)'),
    (TH, "E has the wrong sign", '                        row |= {"E_eV": -math.log(10) * b,', '                        row |= {"E_eV": math.log(10) * b,'),
    (TH, "the primary tests are corrected by false discovery rate, not Holm",
     '    out.loc[prim, "p_holm"] = multipletests(out.loc[prim, "p"], method="holm")[1]',
     '    out.loc[prim, "p_holm"] = multipletests(out.loc[prim, "p"], method="fdr_bh")[1]'),
    (TH, "a class whose size control failed is read anyway",
     "    failed = set(out[out.kind.eq(\"control\")", "    failed = set() and set(out[out.kind.eq(\"control\")"),
    (TH, "Lehmann's R^2 is computed on the tree's residuals",
     "        r2 = 1 - ((y - X1 @ b_ols) ** 2).sum() / ((y - y.mean()) ** 2).sum()",
     "        r2 = 1 - ((y - X1 @ b_ols) ** 2).sum() / (y ** 2).sum()"),
]


def main() -> int:
    # Slurm's time limit (and kill) send SIGTERM: raise, so that `finally` puts the source back
    signal.signal(signal.SIGTERM, lambda *_: sys.exit("terminated: sources put back"))
    originals = {f: f.read_bytes() for f in {m[0] for m in MUTATIONS}}
    digests = {f: hashlib.sha256(b).hexdigest() for f, b in originals.items()}
    bad = 0
    try:
        for path, name, old, new in MUTATIONS:
            text = originals[path].decode()
            if text.count(old) != 1:
                print(f"NOT APPLIED (line found {text.count(old)} times): {name}"); bad += 1; continue
            path.write_text(text.replace(old, new, 1))
            r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", "--show-capture=no", *TESTS],
                               capture_output=True, text=True)
            failed = [l.split(" - ")[0].removeprefix("FAILED ") for l in r.stdout.splitlines() if l.startswith(("FAILED", "ERROR"))]
            if r.returncode == 0:
                print(f"SURVIVED: {name}"); bad += 1
            else:
                print(f"caught: {name}: {', '.join(failed) or r.stdout.strip().splitlines()[-1]}", flush=True)
            path.write_bytes(originals[path])
    finally:
        for f, b in originals.items():
            f.write_bytes(b)
    back = all(hashlib.sha256(f.read_bytes()).hexdigest() == d for f, d in digests.items())
    print(f"{len(MUTATIONS) - bad} of {len(MUTATIONS)} mutations caught; sources restored: {back}")
    return 0 if bad == 0 and back else 1


if __name__ == "__main__":
    sys.exit(main())
