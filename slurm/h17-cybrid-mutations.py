"""Break each guard of H17's code (`mitoage mitoport cybrid-coev`, `cybrid-esm`, `cybrid-readout`) in turn and check
that a test of tests/test_cybrid_h17.py fails.

Run by slurm/h17-cybrid.sbatch after the tests have passed on the unbroken code. Each mutation replaces one line
of a source file, runs the tests, and puts the file back byte for byte (from a copy: git is not installed on the
compute nodes). Prints one line per mutation, "caught" or "SURVIVED", and exits non-zero if any survived, was not
applied, or a file did not come back.
"""
import hashlib
import subprocess
import sys
from pathlib import Path

P = Path("src/mitoage_pipelines/mitoport")
TESTS = ["tests/test_cybrid_h17.py"]
MUTATIONS = [
    (P / "cybrid_coev.py", "the identity threshold strict", "ident >= identity", "ident > identity"),
    (P / "cybrid_coev.py", "host and donor counted in their own pairs", "        if s in skip or a not in IDX", "        if a not in IDX"),
    (P / "cybrid_coev.py", "sequence weights ignored", "f[IDX[a], IDX[b]] += w[s]", "f[IDX[a], IDX[b]] += 1.0"),
    (P / "cybrid_coev.py", "the loss with its sign turned", "    return pmi(p, host_mt, host_nuc) - pmi(p, donor_mt, host_nuc)",
     "    return pmi(p, donor_mt, host_nuc) - pmi(p, host_mt, host_nuc)"),
    (P / "cybrid_coev.py", "a pair with too few species scored", "            if n < MIN_SPECIES:", "            if False:"),
    (P / "cybrid_coev.py", "an unscored pair counted as scored", '            out["pairs_unscored"] += 1', '            out["pairs_scored"] += 1'),
    (P / "cybrid_esm.py", "a pair over the model's length run", "        if len(seq) > MAX_RESIDUES:", "        if False:"),
    (P / "cybrid_esm.py", "the partner's effect not taken against the protein alone",
     "            vals[name] = -(float(np.mean(with_p)) - alone)", "            vals[name] = -float(np.mean(with_p))"),
    (P / "cybrid_esm.py", "a partner that was not run read as assessed", "            if alone is None or not with_p:", "            if alone is None:"),
    (P / "cybrid_esm.py", "every position masked", "if i is not None and aa in IDX and seq[i] in IDX and aa != seq[i]:", "if i is not None:"),
    (P / "cybrid.py", "a missing score table passed over in silence",
     '            print(f"NOT ASSESSED: {out / f} is missing; its scores are left out of the summary")', "            pass"),
]


def main() -> int:
    originals = {f: f.read_bytes() for f in {m[0] for m in MUTATIONS}}
    digests = {f: hashlib.sha256(b).hexdigest() for f, b in originals.items()}
    bad = 0
    try:
        for path, name, old, new in MUTATIONS:
            text = originals[path].decode()
            if old not in text:
                print(f"NOT APPLIED (line not found): {name}"); bad += 1; continue
            path.write_text(text.replace(old, new, 1))
            r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", "--show-capture=no",
                                *TESTS], capture_output=True, text=True)
            failed = [l.split(" - ")[0].removeprefix("FAILED ") for l in r.stdout.splitlines() if l.startswith(("FAILED", "ERROR"))]
            if r.returncode == 0:
                print(f"SURVIVED: {name}"); bad += 1
            else:
                print(f"caught: {name}: {', '.join(failed) or r.stdout.strip().splitlines()[-1]}")
            path.write_bytes(originals[path])
    finally:
        for f, b in originals.items():
            f.write_bytes(b)
    back = all(hashlib.sha256(f.read_bytes()).hexdigest() == d for f, d in digests.items())
    print(f"{len(MUTATIONS) - bad} of {len(MUTATIONS)} mutations caught; sources restored: {back}")
    return 0 if bad == 0 and back else 1


if __name__ == "__main__":
    sys.exit(main())
