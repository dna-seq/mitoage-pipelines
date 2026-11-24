"""Break each guard of `mitoport nuclear-where` in turn and check that a test of tests/test_nuclear_where.py fails.
Prints one line per break, "caught" or "SURVIVED", and exits non-zero if any survived or was not applied."""
import subprocess
import sys
from pathlib import Path

T = Path("src/mitoage_pipelines/mitoport/nuclear_where.py")
breaks = [
 ("a mitochondrial partner counted as nuclear", "    c = contacts[contacts.partner_is_mt == 0].copy()", "    c = contacts.copy()"),
 ("a residue counted once per difference", "    return out.drop_duplicates([\"host\", \"donor\", \"complex\", \"nuclear_gene\", \"partner_resnum\"])", "    return out"),
 ("disagreeing papers kept", "        if s != \"disagree\":", "        if True:"),
 ("reduced and not reduced exchanged", "        red = (g.status == \"reduced\").to_numpy()", "        red = (g.status != \"reduced\").to_numpy()"),
 ("the share of reduced over the wrong total", "g[\"share_reduced\"] = [r.reduced / n_red.get(r.complex, np.nan) for r in g.itertuples()]", "g[\"share_reduced\"] = [r.reduced / n_not.get(r.complex, np.nan) for r in g.itertuples()]"),
 ("the ranking not by the difference", "    return g.sort_values([\"difference\", \"reduced\"], ascending=False).reset_index(drop=True)", "    return g.sort_values([\"reduced\"], ascending=False).reset_index(drop=True)"),
 ("a pair without a difference left out of the counts", "            x = g.merge(n, on=key, how=\"left\").fillna({\"n\": 0})", "            x = g.merge(n, on=key, how=\"inner\")"),
]
caught = 0
for name, a, b in breaks:
    orig = T.read_text()
    if orig.count(a) != 1:
        print("NOT APPLIED", name); continue
    try:
        T.write_text(orig.replace(a, b))
        r = subprocess.run(["uv", "run", "--frozen", "pytest", "-q", "-x", "-p", "no:cacheprovider", "tests/test_nuclear_where.py"], capture_output=True, text=True)
    finally:
        T.write_text(orig)
    ok = r.returncode != 0
    caught += ok
    print("caught " if ok else "SURVIVED", name)
print(caught, "of", len(breaks))
sys.exit(0 if caught == len(breaks) else 1)
