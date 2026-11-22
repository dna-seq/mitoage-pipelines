"""Break each guard of H28's code (`mitoage mitoport cybrid-stab`, its collection and the rule's verdict) in turn
and check that a test of tests/test_cybrid_stab.py fails.

Run by slurm/h28-stab-check.sbatch after the tests have passed on the unbroken code, in the environment with the model,
on a clone that holds the structures, the records and the sites (one break is caught only by the test on a real
structure). Each mutation replaces one line of a source file, runs the tests, and puts the file back byte for byte.
Prints one line per mutation, "caught" or "SURVIVED", and exits non-zero if any survived, was not applied, or a file did
not come back.
"""
import hashlib
import subprocess
import sys
from pathlib import Path

F = Path("src/mitoage_pipelines/mitoport/cybrid_stab.py")
C = Path("src/mitoage_pipelines/mitoport/cybrid.py")
TESTS = ["tests/test_cybrid_stab.py"]
MUTATIONS = [
    (F, "a far chain kept", '        elif hits[i] is None and np.isfinite(tree.query(c["heavy"], distance_upper_bound=CONTACT_A)[0]).any():', "        elif hits[i] is None:"),
    (F, "another complex's mtDNA-encoded chain kept", "        elif hits[i] is None and np.isfinite(", "        elif np.isfinite("),
    (F, "a gene without a chain passed over", "    if not set(mt_proteins) <= found:", "    if False:"),
    (F, "a break not kept", '            gap = min(max(n - chain["num"][i - 1] - 1, 0), GAP_CAP)', "            gap = 0"),
    (F, "a break not capped", '            gap = min(max(n - chain["num"][i - 1] - 1, 0), GAP_CAP)', '            gap = max(n - chain["num"][i - 1] - 1, 0)'),
    (F, "a later chain's places not shifted", "        total.append(seq); offset += len(seq)", "        total.append(seq)"),
    (F, "a substitution at a break accepted", "        if s[i] not in THREE or aa not in THREE:", "        if False:"),
    (F, "a site's residue not checked", "            if human[g][r.human_pos - 1] != r.human_aa:", "            if False:"),
    (F, "sites not carried by the alignment", "        out[g] = {to_host[r.human_pos - 1]: r.group for r in s.itertuples() if r.human_pos - 1 in to_host}",
     "        out[g] = {r.human_pos - 1: r.group for r in s.itertuples()}"),
    (F, "every substitution a contact", "            grp = groups[g].get(to_protein[i])", '            grp = "contact"'),
    (F, "positions not modelled not counted", "        unmodelled += u\n        for i, aa in subs.items():", "        for i, aa in subs.items():"),
    (F, "the sign not turned", "    whole, one, two = (-np.mean(vals, axis=0)).tolist()", "    whole, one, two = np.mean(vals, axis=0).tolist()"),
    (F, "the partners' terms not taken off", '"ddG": whole - (one + two),', '"ddG": whole - one,'),
    (F, "a protein's set not filled", "            sets[GENE + g][(c, i)] = aa\n", ""),
    (F, "a protein's set scored with every weights", 'and sets[k] and name == "stabddg")]:', "and sets[k])]:"),
    (F, "an empty set whose stability term is not 0 accepted", '    bad = long[(long.set == "none") & ((long.ddG != 0) | (long.fold != 0))]', '    bad = long[(long.set == "none") & (long.ddG != 0)]'),
    (F, "a protein's set entered as a score of the complex", '    long = long[(long.set != "none") & ~long.set.str.startswith(GENE)]', '    long = long[long.set != "none"]'),
    (F, "a missing task passed over in silence", '        print(f"NOT ASSESSED: {len(missing)} of {len(want)} tasks have no table: " + "; ".join(missing))', "        pass"),
    (F, "an empty set that is not 0 accepted", "    if len(bad):", "    if False:"),
    (F, "the empty set kept as a score", '    long = long[(long.set != "none") & ~long.set.str.startswith(GENE)]', "    long = long[~long.set.str.startswith(GENE)]"),
    (F, "the spread over the seeds as a population's", 'g.std(ddof=1).unstack("score")', 'g.std(ddof=0).unstack("score")'),
    (C, "T1 equal to the baseline's beats it", "        ok1 = bool(t1 < b1)", "        ok1 = bool(t1 <= b1)"),
    (C, "T2 equal to the baseline's beats it", "        ok2 = bool(t2 > b2)", "        ok2 = bool(t2 >= b2)"),
    (C, "four pairs are enough for T3", "H28_T3_MIN = 5", "H28_T3_MIN = 4"),
    (C, "a missing readout does not stop a verdict", '"verdict": "not assessed" if None in (ok1, ok2, ok3) else', '"verdict":'),
    (C, "a pair without a T3 value counted as not named", "        ok3 = bool(t3.sum() >= H28_T3_MIN) if len(t3) and t3.notna().all() else None", "        ok3 = bool(t3.sum() >= H28_T3_MIN)"),
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
