"""Break each guard of H15's code (`mitoage peptides ...`) in turn and check that a test of tests/test_peptides.py fails.

Run by slurm/job4-peptides.sbatch after the tests have passed on the unbroken code. Each mutation replaces one line
of a source file, runs the tests, and puts the file back byte for byte (from a copy: git is not installed on the
compute nodes). Prints one line per mutation, "caught" or "SURVIVED", and exits non-zero if any survived, was not
applied, or a file did not come back.
"""
import hashlib
import subprocess
import sys
from pathlib import Path

P = Path("src/mitoage_pipelines/peptides")
TESTS = ["tests/test_peptides.py"]
MUTATIONS = [
    (P / "frames.py", "open over half the codons, not 90 %", "OPEN_SHARE = 0.9 ", "OPEN_SHARE = 0.5 "),
    (P / "frames.py", "a gapped window assessed", "gap = (cov < COVERAGE_MIN) | ~inside", "gap = ~inside"),
    (P / "frames.py", "an early stop with a start not counted as lost",
     "state[decided & (early | ~start)] = LOST", "state[decided & ~start] = LOST"),
    (P / "frames.py", "an unknown base read as a sense codon",
     "amb = ~gap & ~early & (k_bad < m)", "amb = ~gap & ~early & (k_bad < 0)"),
    (P / "frames.py", "the window read from the column after the start", "p = host.cum[:, c0]", "p = host.cum[:, c0 + 1]"),
    (P / "codes.py", "AGA and AGG kept as stops in the AGR reading", 'if reading.endswith("agr_r"):', "if False:"),
    (P / "codes.py", "the standard reading starts at TTG and CTG too", 'starts = {"ATG"}', "starts = set(t.start_codons)"),
    (P / "shuffle.py", "codons shuffled across amino acids", "pos = np.flatnonzero(aa == a)", "pos = np.flatnonzero(aa >= 0)"),
    (P / "shuffle.py", "one shuffle repeated in every replicate",
     'perm = np.argsort(rng.random((replicates, len(pos))), axis=1, kind="stable")',
     'perm = np.argsort(rng.random((1, len(pos))), axis=1, kind="stable").repeat(replicates, 0)'),
    (P / "shuffle.py", "the last exit not walked last", "lst.append(tail)", "lst.insert(0, tail)"),
    (P / "shuffle.py", "the frame read one base off", "b = bases[..., offset:offset + 3 * m].astype(np.int16)",
     "b = bases[..., offset + 1:offset + 1 + 3 * m].astype(np.int16)"),
    (P / "nulls.py", "the peptide's own window in its null", "if (q - pc.q0) % 3 == 0 and abs(q - pc.q0) < 3 * n:",
     "if False:"),
    (P / "nulls.py", "null windows in any frame of a protein gene", "if protein and q % 3 != pc.q0 % 3:", "if False:"),
    (P / "nulls.py", "p without the observation counted", 'p": float((1 + (v >= o - 1e-12).sum()) / (1 + replicates))',
     'p": float((v >= o - 1e-12).sum() / replicates)'),
    (P / "nulls.py", "q over everything, not within null, statistic, level and reading",
     'keys = ["null", "statistic", "level", "reading"]', 'keys = ["status"]'),
    (P / "tree.py", "a missing tip read as lost", "s = states.get(name)", "s = states.get(name, np.zeros(width, dtype=int))"),
    (P / "tree.py", "a change preferred on a tie", "t = s if stay <= move else 1 - s", "t = s if stay < move else 1 - s"),
    (P / "lifespan.py", "species not assessed fitted as lost",
     '(states.state != "not assessed")]', '(states.state != "-")]'),
    (P / "catalogue.py", "the published nucleotides not checked",
     "if pep.nucleotides and not nt.startswith(pep.nucleotides[:3 * n]):", "if False:"),
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
