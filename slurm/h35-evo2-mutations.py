"""Break each guard of H35's table logic (`mitoport/evo_table.py`) in turn and check that a test of
tests/test_evo_table.py fails. Prints one line per break, "caught" or "SURVIVED", and exits non-zero if any survived
or was not applied."""
import subprocess
import sys
from pathlib import Path

T = Path("src/mitoage_pipelines/mitoport/evo_table.py")
breaks = [
 ("always the reading from the start", "return np.where((before_first >= before_second)[:, None], first, second)", "return first"),
 ("the reading with fewer bases before the place", "(before_first >= before_second)[:, None], first, second)", "(before_first < before_second)[:, None], first, second)"),
 ("the turned reading not turned back", "second = np.roll(turned, half, axis=0)", "second = turned"),
 ("a wrong shape taken", "        if got.shape != (n, 4):", "        if False:"),
 ("the minus strand's places not reversed", "    return minus[::-1, ::-1]", "    return minus[:, ::-1]"),
 ("the minus strand's bases not complemented", "    return minus[::-1, ::-1]", "    return minus[::-1, :]"),
 ("the reverse complement without the reverse", "    return seq.translate(COMPLEMENT)[::-1]", "    return seq.translate(COMPLEMENT)"),
 ("the plus strand alone as both", "\"both\": normalised((plus + minus) / 2)}", "\"both\": normalised(plus)}"),
 ("the mean not normalised", "\"both\": normalised((plus + minus) / 2)}", "\"both\": (plus + minus) / 2}"),
 ("an unknown base scored", "        t[unknown] = np.nan\n", "        pass\n"),
 ("lower case not read", "    seq = seq.upper()\n", ""),
]
caught = 0
for name, a, b in breaks:
    orig = T.read_text()
    if orig.count(a) != 1:
        print("NOT APPLIED", name); continue
    try:
        T.write_text(orig.replace(a, b))
        r = subprocess.run(["uv", "run", "--frozen", "pytest", "-q", "-x", "-p", "no:cacheprovider", "tests/test_evo_table.py"], capture_output=True, text=True)
    finally:
        T.write_text(orig)
    ok = r.returncode != 0
    caught += ok
    print("caught " if ok else "SURVIVED", name)
print(caught, "of", len(breaks))
sys.exit(0 if caught == len(breaks) else 1)
