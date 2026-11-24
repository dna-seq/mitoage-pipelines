"""Break each guard of `mitoport where` in turn and check that a test of tests/test_where.py fails. Prints one line
per break, "caught" or "SURVIVED", and exits non-zero if any survived or was not applied."""
import subprocess
import sys
from pathlib import Path

T = Path("src/mitoage_pipelines/mitoport/where.py")
breaks = [
 ("an unknown flag counted as no", "            d = v[(v.label == label) & v[flag].notna()]", "            d = v[(v.label == label)]"),
 ("the labels not kept apart", "            d = v[(v.label == label) & v[flag].notna()]", "            d = v[v[flag].notna()]"),
 ("the bins of the site's amino acids shifted", "k = int(((d.n_aa >= lo) & (d.n_aa <= hi)).sum())", "k = int(((d.n_aa > lo) & (d.n_aa <= hi)).sum())"),
 ("equal residues listed as differences", "                if \"?\" in (a, b) or a == b:", "                if \"?\" in (a, b):"),
 ("host and donor exchanged", "\"host_aa\": a, \"donor_aa\": b,", "\"host_aa\": b, \"donor_aa\": a,"),
 ("every difference counted as conserved", "        cons = x[x.conserved == 1]", "        cons = x"),
 ("buried not asked of the conserved", "int((cons.buried == 1).sum())", "int((x.buried == 1).sum())"),
]
caught = 0
for name, a, b in breaks:
    orig = T.read_text()
    if orig.count(a) != 1:
        print("NOT APPLIED", name); continue
    try:
        T.write_text(orig.replace(a, b))
        r = subprocess.run(["uv", "run", "--frozen", "pytest", "-q", "-x", "-p", "no:cacheprovider", "tests/test_where.py"], capture_output=True, text=True)
    finally:
        T.write_text(orig)
    ok = r.returncode != 0
    caught += ok
    print("caught " if ok else "SURVIVED", name)
print(caught, "of", len(breaks))
sys.exit(0 if caught == len(breaks) else 1)
