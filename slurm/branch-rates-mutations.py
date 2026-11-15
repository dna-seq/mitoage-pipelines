"""Break each piece of H3's two-step check (`longevity branch-rates`, `branch-pgls`) in turn and check that a test
fails. Each mutation replaces one line of src/mitoage_pipelines/longevity/branch_rates.py, runs
tests/test_branch_rates.py, and puts the file back byte for byte. Prints "caught" or "SURVIVED" per mutation and
exits non-zero if any survived or the file did not come back."""
import hashlib
import subprocess
import sys
from pathlib import Path

FILE = Path("src/mitoage_pipelines/longevity/branch_rates.py")
TESTS = ["tests/test_branch_rates.py"]
MUTATIONS = [
    ("dN/dS inverted", "df.beta / df.alpha.where(df.alpha > 0, np.nan)", "df.alpha.where(df.alpha > 0, np.nan) / df.beta"),
    ("dS from S alone", '    df["dS"] = df.alpha * (S + N)\n', '    df["dS"] = df.alpha * S\n'),
    ("dN from N alone", '    df["dN"] = df.beta * (S + N)\n', '    df["dN"] = df.beta * N\n'),
    ("lengths not checked against the model", "    if worst > 1e-6 * max(1.0, df.length.max()):\n", "    if False:\n"),
    ("root branches not merged", "            if s in fitted_sets and other not in fitted_sets:\n", "            if False:\n"),
    ("missing species not refused", "    if missing:\n", "    if False:\n"),
    ("dS per Myr not divided by the duration", "df.dS / df.myr.where(df.myr > 0, np.nan)", "df.dS"),
    ("lifespan read as mass", '    df["ll"] = np.log10(df.longevity.astype(float))\n', '    df["ll"] = np.log10(df.mass.astype(float))\n'),
    ("lambda not estimated", "            lam = float(optimize.minimize_scalar(nll, bounds=(0.0, 1.0), args=(full, y, C), method=\"bounded\").x)\n",
     "            lam = 1.0\n"),
    ("tree rows not matched to species", "            idx = np.array([order[s] for s in d.species])\n", "            idx = np.arange(len(d))\n"),
    ("saturated branches kept", "                d = d[d.dS <= SATURATED_DS]\n", "                pass\n"),
    ("no-log branches not counted", "            dropped = int(d[resp].isna().sum())\n", "            dropped = 0\n"),
    ("Holm over every row", '    q = (res.branches == "all") & (res.model == "lifespan+mass") & (res.term == "ll")\n',
     '    q = res.branches == "all"\n'),
    ("control not required", "    control = c.slope > 0 and c.p < 0.05\n", "    control = True\n"),
    ("Coevol's sign ignored", '            verdict = "agrees with Coevol\'s sign" if np.sign(r.slope) == sign else "DISAGREES with Coevol\'s sign"\n',
     '            verdict = "agrees with Coevol\'s sign"\n'),
]


def main() -> int:
    orig = FILE.read_bytes(); digest = hashlib.sha256(orig).hexdigest(); text = orig.decode()
    survived = 0
    try:
        for name, old, new in MUTATIONS:
            if text.count(old) != 1:
                print(f"BAD MUTATION {name}: the line occurs {text.count(old)} times"); survived += 1; continue
            FILE.write_text(text.replace(old, new))
            try:
                r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", *TESTS],
                                   capture_output=True, text=True, timeout=600)
                caught = r.returncode != 0
            except subprocess.TimeoutExpired:
                caught = True
            survived += not caught
            print(f"{'caught' if caught else 'SURVIVED'}: {name}")
    finally:
        FILE.write_bytes(orig)
    back = hashlib.sha256(FILE.read_bytes()).hexdigest() == digest
    print(f"{len(MUTATIONS) - survived} of {len(MUTATIONS)} caught; source restored: {back}")
    return 0 if survived == 0 and back else 1


if __name__ == "__main__":
    sys.exit(main())
