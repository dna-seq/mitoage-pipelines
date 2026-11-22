"""Break each guard of H34's comparison with the existing predictors (`mitoage mitoport variant-compare`) in turn
and check that a test of tests/test_variant_compare.py fails. Each mutation replaces one piece of the source, runs
the tests, and the file is put back. Prints one line per mutation and exits non-zero if any survived or was not
applied.

    uv run --frozen --group dev python slurm/h34-compare-mutations.py
"""
import subprocess
import sys
from pathlib import Path

F = Path("src/mitoage_pipelines/mitoport/variant_compare.py")
MUTATIONS = [
    ("ESM1b's sign not turned", '"esm1b": ("ESM1b_score", -1,', '"esm1b": ("ESM1b_score", 1,'),
    ("another file than the pinned one read", "    if digest is not None and sha256(csv) != digest:", "    if False:"),
    ("changes that are not missense kept", '    d = d[d.functional_effect_general == "missense_variant"]\n', ""),
    ("an unknown gene passed over", "    if out.gene.isna().any():", "    if False:"),
    ("a change twice for one gene passed over", '    if out.duplicated(["pos", "ref", "alt", "gene"]).any():', "    if False:"),
    ("the join without the gene", 'v = variants.merge(mitimpact, on=["pos", "ref", "alt", "gene"], how="left", validate="one_to_one")',
     'v = variants.merge(mitimpact.drop_duplicates(["pos", "ref", "alt"]).drop(columns="gene"), on=["pos", "ref", "alt"], how="left")'),
    ("a missing score read as zero", '    v["harmful"] = v.label == "harmful"', '    v = v.fillna({k: 0 for k in list(PREDICTORS) + ["kron"]}); v["harmful"] = v.label == "harmful"'),
    ("the baseline upside down", "    v[BASELINE] = -v.n_aa", "    v[BASELINE] = v.n_aa"),
    ("the area over all variants, covered or not", "        d = v[v[s].notna() & v[BASELINE].notna()]", "        d = v.fillna({s: -1e9})"),
    ("one kind only still gets an area", "        if h.any() and not h.all():", "        if len(h):"),
]
original = F.read_text()
caught = 0
try:
    for what, old, new in MUTATIONS:
        if original.count(old) != 1:
            raise SystemExit(f"{what}: the piece to break is not there exactly once")
        F.write_text(original.replace(old, new))
        failed = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", "tests/test_variant_compare.py"], capture_output=True).returncode != 0
        caught += failed
        print(("caught: " if failed else "SURVIVED: ") + what)
finally:
    F.write_text(original)
print(f"{caught} of {len(MUTATIONS)} mutations caught; sources restored: {F.read_text() == original}")
sys.exit(0 if caught == len(MUTATIONS) else 1)
