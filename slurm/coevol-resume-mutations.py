"""Break each guard of `longevity coevol-resume` in turn and check that a test fails.

Each mutation replaces one line of src/mitoage_pipelines/longevity/coevol.py, runs tests/test_coevol_resume.py,
and puts the file back byte for byte (from a copy: git is not installed on the compute nodes). Prints one line
per mutation, "caught" or "SURVIVED", and exits non-zero if any survived or the file did not come back.
"""
import hashlib
import subprocess
import sys
from pathlib import Path

FILE = Path("src/mitoage_pipelines/longevity/coevol.py")
TESTS = ["tests/test_coevol_resume.py"]
MUTATIONS = [
    ("param size not compared with the trace",
     "    if p != n:\n", "    if False:\n"),
    ("chain lines not compared with the points",
     "    if n <= 0 or lines % n:\n", "    if n <= 0:\n"),
    ("checkpoint of an inconsistent chain",
     "    if not ok:\n        raise RuntimeError(f\"chain{chain} is not consistent", "    if False:\n        raise RuntimeError(f\"chain{chain} is not consistent"),
    ("restore without truncating the chain",
     "        f.truncate(state[\"chain_bytes\"])\n", "        pass\n"),
    ("restore without checking the prefix",
     "    if cf.stat().st_size < state[\"chain_bytes\"] or _sha256_prefix(cf, state[\"chain_bytes\"]) != state[\"chain_sha256\"]:\n",
     "    if cf.stat().st_size < state[\"chain_bytes\"]:\n"),
    ("restore from the oldest checkpoint", "    d = cps[-1]\n", "    d = cps[0]\n"),
    ("old checkpoints not pruned", "    for old in sorted(p for p in root.iterdir() if p.is_dir())[:-keep]:\n",
     "    for old in []:\n"),
    ("resume without the consistency check",
     "    ok, why, _ = consistency(out, chain)\n    if not ok:\n        print(", "    ok, why = True, ''\n    if not ok:\n        print("),
    ("no checkpoint after a segment", "        run_segment(cmd, out, chain, hours, until, grace_s, poll_s)\n        checkpoint(out, chain)\n",
     "        run_segment(cmd, out, chain, hours, until, grace_s, poll_s)\n"),
    ("a run that stopped without exiting is not ended",
     "                elif (asked or n > until) and now - since >= grace_s and consistency(out, chain)[0]:\n",
     "                elif False:\n"),
    ("a run ended while its files disagree",
     "                elif (asked or n > until) and now - since >= grace_s and consistency(out, chain)[0]:\n",
     "                elif (asked or n > until) and now - since >= grace_s:\n"),
    ("the stop not written again after Coevol starts",
     "                if asked and run_file.read_text().strip() != \"0\":\n", "                if False:\n"),
    ("time budget ignored", "            if hours < min_segment_hours:\n", "            if False:\n"),
    ("until not written before a segment",
     "        param.write_text(set_until(param.read_text(), size, until))\n", "        pass\n"),
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
            except subprocess.TimeoutExpired:     # the tests' own 60 s alarm should come first; a hang is a catch
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
