"""Break each guard of the audit's code (H13) in turn and check that a test of tests/test_audit.py fails.

Run by slurm/job3-audit.sbatch after the tests have passed on the unbroken code. Each mutation replaces one line of a
source file, runs the tests, and puts the file back byte for byte from a copy (git is not installed on the compute
nodes). Prints one line per mutation, "caught" or "SURVIVED", and exits non-zero if any survived, a line was not
found, or a file did not come back.
"""
import hashlib
import subprocess
import sys
from pathlib import Path

A = Path("src/mitoage_pipelines/audit")
TESTS = ["tests/test_audit.py"]
MUTATIONS = [
    (A / "strand.py", "ND6 not on the other strand", 'return -1 if gene == "ND6" else 1', "return 1"),
    (A / "strand.py", "too few genes still classified", "if len(present) < MIN_GENES:", "if len(present) < 1:"),
    (A / "strand.py", "the convention applied outside Vertebrata", "    if not vertebrate:", "    if False:"),
    (A / "strand.py", "no reversed class", "    if not agree:", "    if False:"),
    (A / "records.py", "every record called vertebrate", '"vertebrate": "Vertebrata" in h.lineage.split("; "),', '"vertebrate": True,'),
    (A / "records.py", "reverse complement keeps the strand", "gm.Segment(n - s.end + 1, n - s.start + 1, -s.strand)",
     "gm.Segment(n - s.end + 1, n - s.start + 1, s.strand)"),
    (A / "order.py", "the other strand read without reversing the order",
     'return [("-" if t[0] == "+" else "+") + t[1:] for t in reversed(tok)]', 'return [("-" if t[0] == "+" else "+") + t[1:] for t in tok]'),
    (A / "order.py", "not canonical by rotation", '        return " ".join(tok[i:] + tok[:i])', '        return " ".join(tok)'),
    (A / "order.py", "families of two flagged", "MIN_FAMILY = 3", "MIN_FAMILY = 2"),
    (A / "order.py", "a tRNA-less record given a full arrangement",
     'n_trna = sum(1 for t in full if t[1:].startswith("tRNA"))', "n_trna = len(full)"),
    (A / "order.py", "a reversed deposit called a rearrangement",
     '    if own[f"order_{tier}_unoriented"] == other[f"order_{tier}_unoriented"]:', "    if False:"),
    (A / "order.py", "strand differences not told apart",
     '    if canonical(["+" + x for x in _names(a)]) == canonical(["+" + x for x in _names(b)]):', "    if False:"),
    (A / "code.py", "the terminal codon counted", "    for i in range(n - 1):", "    for i in range(n):"),
    (A / "code.py", "one clean code is not enough", "    if YES in (clean_declared, clean_expected):", "    if YES == clean_declared:"),
    (A / "code.py", "an unknown code read as clean", "    return UNKNOWN if stops is None else YES if stops == 0 else NO",
     "    return YES if not stops else NO"),
    (A / "code.py", "merged taxids not followed", "        cur = self.merged.get(t, t)", "        cur = t"),
    (A / "code.py", "codon_start ignored", '"nt": gm.extract(seq, f.segments)[cs - 1:].upper(),', '"nt": gm.extract(seq, f.segments).upper(),'),
    (A / "similar.py", "compared columns counted as matches", "    compared = valid @ valid.T", "    compared = match"),
    (A / "similar.py", "a looser identity threshold", "NEAR = 0.99", "NEAR = 0.5"),
    (A / "similar.py", "own order not required", "            same = (tax[lv].to_numpy() == own) & (ident[i] >= 0)",
     "            same = ident[i] >= 0"),
    (A / "common.py", "tables not sorted", '    out = df.sort_values(keys, ascending=ascending, kind="mergesort") if keys else df', "    out = df"),
    (A / "common.py", "warnings in the order raised", "for reason, n in sorted(self.count.items())]", "for reason, n in self.count.items()]"),
]


def main() -> int:
    originals = {f: f.read_bytes() for f in {m[0] for m in MUTATIONS}}
    digests = {f: hashlib.sha256(b).hexdigest() for f, b in originals.items()}
    survived = 0
    try:
        for path, name, old, new in MUTATIONS:
            text = originals[path].decode()
            if old not in text:
                print(f"NOT APPLIED (line not found): {name}", flush=True); survived += 1; continue
            path.write_text(text.replace(old, new, 1))
            r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "--show-capture=no", *TESTS],
                               capture_output=True, text=True)
            failed = [l.split(" - ")[0].removeprefix("FAILED ").removeprefix("ERROR ") for l in r.stdout.splitlines()
                      if l.startswith(("FAILED", "ERROR"))]
            if r.returncode == 0:
                print(f"SURVIVED: {name}", flush=True); survived += 1
            else:
                print(f"caught: {name}: {', '.join(failed) or r.stdout.strip().splitlines()[-1]}", flush=True)
            path.write_bytes(originals[path])
    finally:
        for f, b in originals.items():
            f.write_bytes(b)
    back = all(hashlib.sha256(f.read_bytes()).hexdigest() == d for f, d in digests.items())
    print(f"{len(MUTATIONS) - survived} of {len(MUTATIONS)} mutations caught; sources restored: {back}")
    return 0 if survived == 0 and back else 1


if __name__ == "__main__":
    sys.exit(main())
