"""Break each guard of the nucleotide baselines (`mitoage mitoport cybrid-dna`) in turn and check that a test of
tests/test_cybrid_dna.py fails. Run on a clone that holds the records of `download cybrid` (four breaks are caught
only on them). Each mutation replaces one piece of the source, runs the tests, and the file is put back at the end.
Prints one line per mutation, "caught" or "SURVIVED", and exits non-zero if any survived or was not applied."""
import subprocess, sys
from pathlib import Path
p = Path("src/mitoage_pipelines/mitoport/cybrid_dna.py"); orig = p.read_text()
breaks = [
 ("a gap not counted", "return sum(x != y for x, y in cols), len(cols)", "return sum(x != y and '-' not in (x, y) for x, y in cols), len(cols)"),
 ("the place one early", "else start < i <= end ", "else start <= i < end "),
 ("an insertion inside dropped", 'inside = start < i < end if x == "-" else', 'inside = False if x == "-" else'),
 ("an insertion after the last place kept", 'inside = start < i < end if x == "-" else', 'inside = start < i <= end if x == "-" else'),
 ("the own line ignored", "acc = cybrid.OWN_LINE.get((slug, gene), cybrid.GENOMES[slug])", "acc = cybrid.GENOMES[slug]"),
 ("humanin's frame moved", "HUMANIN = (2633, 2707)", "HUMANIN = (2634, 2708)"),
 ("the donor compared with itself", "h, d = features[cybrid.GENOMES[host]], features[cybrid.GENOMES[donor]]", "h = d = features[cybrid.GENOMES[donor]]"),
 ("every changed codon silent", "if pa[a0 + k] == pb[b0 + k]:", "if True:"),
 ("the standard code", "pa, pb = (str(Seq(\"\".join(c)).translate(table=2)) for c in (ca, cb))", "pa, pb = (str(Seq(\"\".join(c)).translate(table=1)) for c in (ca, cb))"),
 ("an unknown base compared", "if set(x + y) - set(\"ACGT\"):", "if False:"),
 ("a final stop kept", "    if out and str(Seq(out[-1]).translate(table=2)) == \"*\":\n        out.pop()\n", ""),
 ("the blocks without the host's start", "carried(hs, ds, a - start, b - start, al)", "carried(hs, ds, a, b, al)"),
 ("a missing control region scored", "    if h[\"dloop\"] is None or d[\"dloop\"] is None:\n        return out\n", "    if h[\"dloop\"] is None or d[\"dloop\"] is None:\n        return dict(out, nt_dloop_d=0)\n"),
 ("the control region in the record's order of parts", "key=lambda q: int(q.start) == 0", "key=lambda q: 0"),
 ("tRNA genes not numbered", "i = seen.get(product, 0); seen[product] = i + 1", "i = 0"),
]
caught = 0
try:
    for name, a, b in breaks:
        assert orig.count(a) == 1, name
        p.write_text(orig.replace(a, b))
        r = subprocess.run(["uv", "run", "--frozen", "pytest", "-q", "-x", "tests/test_cybrid_dna.py"], capture_output=True, text=True)
        ok = r.returncode != 0
        caught += ok
        print("caught " if ok else "SURVIVED", name)
finally:
    p.write_text(orig)
print(caught, "of", len(breaks))
sys.exit(0 if caught == len(breaks) and p.read_text() == orig else 1)
