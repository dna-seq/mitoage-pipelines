"""Break each guard of H18's code (`mitoage mitoport cybrid-energy`) in turn and check that a test of
tests/test_cybrid_energy.py fails.

Run by slurm/h18-energy-check.sbatch after the tests have passed on the unbroken code, in the environment with OpenMM.
Each mutation replaces one line of a source file, runs the tests, and puts the file back byte for byte. Prints one
line per mutation, "caught" or "SURVIVED", and exits non-zero if any survived, was not applied, or a file did not
come back.
"""
import hashlib
import subprocess
import sys
from pathlib import Path

F = Path("src/mitoage_pipelines/mitoport/cybrid_energy.py")
TESTS = ["tests/test_cybrid_energy.py"]
MUTATIONS = [
    (F, "a chain matched at any identity", "        if ident >= MIN_IDENTITY and (best is None or ident > best[0]):", "        if (best is None or ident > best[0]):"),
    (F, "a short piece matched", "    if len(seq) < 20:", "    if len(seq) < 2:"),
    (F, "a deletion not counted", "        else:\n            unmodelled += 1", "        else:\n            unmodelled += 0"),
    (F, "an unresolved difference not counted", "    unmodelled += sum(1 for p, d in enumerate(donor_on_host) if p not in resolved and d != host[p])", "    unmodelled += 0"),
    (F, "the structure's own residue substituted again", "            if d != struct_seq[i]:\n                subs[i] = d", "            subs[i] = d"),
    (F, "a chain not cut at a break", "            if not joined:", "            if cur is None:"),
    (F, "backbone atoms left free", "    free = hv & side & in_free", "    free = hv & in_free"),
    (F, "every residue free", "    in_free = np.array([tuple(k) in free_res for k in key])", "    in_free = np.ones(len(key), dtype=bool)"),
    (F, "a missing task passed over in silence", '        print(f"NOT ASSESSED: {len(missing)} of {len(want)} tasks have no table: " + "; ".join(missing))', "        pass"),
    (F, "an mtDNA-encoded chain of another complex kept", "            drop.add(c.index)", "            pass"),
    (F, "the parts' energies not taken off", "    return total - parts[0] - parts[1]", "    return total - parts[0]"),
    (F, "substitutions not made", "        fx.applyMutations(muts, chains[c].id)", "        pass"),
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
