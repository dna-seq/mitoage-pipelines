"""Break each guard of H32's code (`mitoage mitoport cybrid-blocks`) in turn and check that a test of
tests/test_cybrid_blocks.py fails. Run on a clone that holds the records of `download cybrid`. Each mutation replaces
one piece of the source, runs the tests, and the file is put back at the end. Prints one line per mutation, "caught"
or "SURVIVED", and exits non-zero if any survived or was not applied."""
import subprocess
import sys
from pathlib import Path

p = Path("src/mitoage_pipelines/mitoport/cybrid_blocks.py"); orig = p.read_text()
breaks = [
 ("a deletion not counted", "out[i] += x != y", "out[i] += x != y and y != '-'"),
 ("an insertion not counted", "out[max(i, 0)] += 1", "pass"),
 ("windows may overlap", "if all(a + s <= b for (a, s), (b, _) in zip(spans, spans[1:])):", "if True:"),
 ("a window past the end", "starts = [int(rng.integers(0, length - s + 1)) for s in sizes]", "starts = [int(rng.integers(0, length)) for s in sizes]"),
 ("a window one short", "c[starts[:, j] + s] - c[starts[:, j]]", "c[starts[:, j] + s - 1] - c[starts[:, j]]"),
 ("the area upside down", "red = (use.status == \"reduced\").to_numpy()", "red = (use.status != \"reduced\").to_numpy()"),
 ("a disagreeing point kept", "use = points[points.status != \"disagree\"]", "use = points"),
 ("the spearman on the scores, not their ranks", "rt = stats.rankdata(t, axis=0); rp", "rt = t; rp"),
 ("the mouse blocks as places of the record", "return {k: (a - start, b - start) for k, (a, b) in sorted(h[\"blocks\"].items())}", "return {k: (a, b) for k, (a, b) in sorted(h[\"blocks\"].items())}"),
 ("the human blocks from the whole mouse region", "block = mseq[a - mstart:b - mstart]", "block = mseq[a - mstart:b - mstart + 40]"),
 ("the donor compared with the host's own", "h, d = features[cybrid.GENOMES[host]], features[cybrid.GENOMES[donor]]", "h = d = features[cybrid.GENOMES[host]]"),
]
caught = 0
try:
    for name, a, b in breaks:
        assert orig.count(a) == 1, name
        p.write_text(orig.replace(a, b))
        r = subprocess.run(["uv", "run", "--frozen", "pytest", "-q", "-x", "tests/test_cybrid_blocks.py"], capture_output=True, text=True)
        ok = r.returncode != 0
        caught += ok
        print("caught " if ok else "SURVIVED", name)
finally:
    p.write_text(orig)
print(caught, "of", len(breaks))
sys.exit(0 if caught == len(breaks) and p.read_text() == orig else 1)
