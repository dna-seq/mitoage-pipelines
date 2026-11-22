"""Break each guard of H31's scoring (`mitoage mitoport kron-cybrids`, `kron-variants`, and the likelihood table's
logic) in turn and check that a test of tests/test_kron.py or tests/test_kron_table.py fails. Run on a clone that
holds the records of `download cybrid`. Each mutation replaces one piece of a source file, runs the tests, and the
file is put back. Prints one line per mutation, "caught" or "SURVIVED", and exits non-zero if any survived or was
not applied. Not among them: the order in which overlapping kinds of place are named, since the human
record has no place in both a protein gene and another annotated stretch, so no test could tell."""
import subprocess
import sys
from pathlib import Path

K = Path("src/mitoage_pipelines/mitoport/kron.py")
T = Path("src/mitoage_pipelines/mitoport/kron_table.py")
breaks = [
 (K, "the drop upside down", "return float(lp[place, BASES.index(ref)] - lp[place, BASES.index(alt)])", "return float(lp[place, BASES.index(alt)] - lp[place, BASES.index(ref)])"),
 (K, "an unknown base scored", "    if ref not in BASES or alt not in BASES:\n        return float(\"nan\")\n", ""),
 (K, "the minus strand not complemented", "ref, alt = (x, y) if strand == 1 else (x.translate(COMPLEMENT), y.translate(COMPLEMENT))", "ref, alt = x, y"),
 (K, "the host's base not checked", "        if genome[places[i]] != ref:", "        if False:"),
 (K, "a gap scored as a change", "        if x == \"-\" or y == \"-\":\n            gaps += 1\n            continue\n", "        if x == \"-\" or y == \"-\":\n            gaps += 1; n += 1\n            continue\n"),
 (K, "the place not advanced past a deletion in the donor", "        if x != \"-\":\n            i += 1\n        if x == \"-\" or y == \"-\":", "        if x != \"-\" and y != \"-\":\n            i += 1\n        if x == \"-\" or y == \"-\":"),
 (K, "the rule passes on two", "\"counts\" if all(beats) else", "\"counts\" if sum(beats) >= 2 else"),
 (K, "a tie with the bar passes", "beats = [t1 < bar[\"T1\"], t2 > bar[\"T2\"], t3 > bar[\"T3\"]]", "beats = [t1 <= bar[\"T1\"], t2 >= bar[\"T2\"], t3 >= bar[\"T3\"]]"),
 (K, "T3 without its floor", "\"T3\": max(T3_MIN - 1, max(value(\"T3\", b) for b in BASELINE))", "\"T3\": max(value(\"T3\", b) for b in BASELINE) - 1"),
 (K, "any reviewer's pathogenic is harmful", "        if r[\"review\"] == variant_stab.PANEL:\n            rows[k] = (\"harmful\", np.nan)", "        if True:\n            rows[k] = (\"harmful\", np.nan)"),
 (K, "a pathogenic variant among the harmless", "        if k in any_pathogenic:", "        if False:"),
 (K, "a rare variant among the harmless", "& (g.AF_hom >= variant_stab.AF_HOM_MIN)]", "& (g.AF_hom >= 0)]"),
 (K, "a filtered variant kept", "g = g[(g.filters == \"PASS\") & ", "g = g[(g.filters != \"\") & "),
 (K, "the control region in the record's order of parts", "parts = sorted(f.location.parts, key=lambda q: int(q.start) == 0)", "parts = list(f.location.parts)"),
 (T, "two places hidden side by side", "return [np.arange(k, n, stride) for k in range(min(stride, n))]", "return [np.arange(k, n, max(stride // 2, 1)) for k in range(min(stride, n))]"),
 (T, "an unknown place scored", "    out[codes == N] = np.nan\n", ""),
 (T, "always the record's reading", "return np.where((edge_first >= edge_second)[:, None], first, second)", "return first"),
 (T, "the turned reading not turned back", "second = np.roll(turned, half, axis=0)", "second = turned"),
 (T, "no padding", "    out = np.full(-(-len(codes) // multiple) * multiple, N, dtype=np.int64)", "    out = np.full(len(codes), N, dtype=np.int64)"),
]
caught = 0
for f, name, a, b in breaks:
    orig = f.read_text()
    if orig.count(a) != 1:
        print("NOT APPLIED", name); continue
    try:
        f.write_text(orig.replace(a, b))
        r = subprocess.run(["uv", "run", "--frozen", "pytest", "-q", "-x", "tests/test_kron.py", "tests/test_kron_table.py"], capture_output=True, text=True)
    finally:
        f.write_text(orig)
    ok = r.returncode != 0
    caught += ok
    print("caught " if ok else "SURVIVED", name)
print(caught, "of", len(breaks))
sys.exit(0 if caught == len(breaks) else 1)
