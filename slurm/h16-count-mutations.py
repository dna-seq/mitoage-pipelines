"""Breaks of the guards added to H16's count after the review of 2026-11-22: each is applied to the source, the
named test must then fail, and the source is restored.

    uv run --frozen --group dev python slurm/h16-count-mutations.py
"""
import subprocess
import sys
from pathlib import Path

SRC = Path("src/mitoage_pipelines/mitoport/cybrid.py")
T = "tests/test_cybrid.py::"
MUTATIONS = [
    ("the count aligns the donor to human by itself again",
     "ds = donor_through_host(human[g], hp, dp, dpart, al)", "ds = on_human(human[g], dp, dpart, al)",
     T + "TestDonorThroughHost::test_the_count_reads_the_donor_through_the_host"),
    ("the donor through the host is the donor's own alignment",
     "return [own[p] if i is None else direct[i] for p, i in enumerate(index_on_human(human, host, al))]", "return own",
     T + "TestDonorThroughHost::test_mouse_and_caroli_atp8_differ_at_five_residues"),
    ("the first residue of a partial record is compared",
     "itertools.chain(range(0, first + 1), range(last - 1, len(human)))", "itertools.chain(range(0, first), range(last - 1, len(human)))",
     T + "TestOnHuman::test_the_first_residue_of_a_partial_sequence_is_not_a_substitution"),
    ("the last residue of a partial record is compared",
     "itertools.chain(range(0, first + 1), range(last - 1, len(human)))", "itertools.chain(range(0, first + 1), range(last, len(human)))",
     T + "TestOnHuman::test_the_ends_of_a_partial_sequence_are_unknown_not_deleted"),
    ("the index on human is off by one",
     "out[h0 + k] = o0 + k", "out[h0 + k] = o0 + k + 1",
     T + "TestDonorThroughHost::test_index_on_human"),
]
original = SRC.read_text()
caught = 0
try:
    for what, old, new, test in MUTATIONS:
        if original.count(old) != 1:
            raise SystemExit(f"{what}: the line to break is not there exactly once")
        SRC.write_text(original.replace(old, new))
        failed = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", test], capture_output=True).returncode != 0
        caught += failed
        print(("caught: " if failed else "NOT CAUGHT: ") + what + ": " + test)
finally:
    SRC.write_text(original)
print(f"{caught} of {len(MUTATIONS)} mutations caught; sources restored: {SRC.read_text() == original}")
sys.exit(0 if caught == len(MUTATIONS) else 1)
