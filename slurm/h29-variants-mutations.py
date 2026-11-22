"""Break each guard of H29's code (`mitoage mitoport variant-stab`, its variants and its readout) in turn and check
that a test of tests/test_variant_stab.py fails.

Run by slurm/h29-variants-check.sbatch after the tests have passed on the unbroken code. Each mutation replaces one
line of a source file, runs the tests, and puts the file back byte for byte. Prints one line per mutation, "caught" or
"SURVIVED", and exits non-zero if any survived, was not applied, or a file did not come back.
"""
import hashlib
import subprocess
import sys
from pathlib import Path

F = Path("src/mitoage_pipelines/mitoport/variant_stab.py")
TESTS = ["tests/test_variant_stab.py"]
MUTATIONS = [
    (F, "the reference base not checked", "    if seq[pos - 1] != ref:", "    if False:"),
    (F, "a gene on the other strand read forwards", "        off = pos - 1 - s if strand == 1 else e - pos", "        off = pos - 1 - s"),
    (F, "the new base not complemented on the other strand", "        a = alt if strand == 1 else alt.translate(COMPLEMENT)", "        a = alt"),
    (F, "a change of the start codon taken as an amino acid change", 'kind = "start" if k == 0 else "synonymous"', 'kind = "synonymous"'),
    (F, "the second codon taken as the start", 'kind = "start" if k == 0 else', 'kind = "start" if k <= 1 else'),
    (F, "a lost stop called nonsense", '"stop lost" if x == "*" else "nonsense" if y == "*"', '"nonsense" if "*" in (x, y)'),
    (F, "the standard genetic code", "CodonTable.unambiguous_dna_by_id[2]", "CodonTable.unambiguous_dna_by_id[1]"),
    (F, "conflicting classifications taken as pathogenic", 'classification.lower().startswith(("pathogenic", "likely pathogenic"))', '"pathogenic" in classification.lower()'),
    (F, "harmful at any review level", '        if r["review"] == PANEL:\n            rows +=', "        if True:\n            rows +="),
    (F, "a pathogenic variant kept among the harmless", "        if key in any_pathogenic:", "        if False:"),
    (F, "gnomAD's filter not applied", '    g = g[(g.filters == "PASS") & ', "    g = g["),
    (F, "rare variants kept as harmless", "(g.AF_hom >= AF_HOM_MIN)]", "(g.AF_hom >= 0)]"),
    (F, "a change that is not missense kept", '                        for k, label, af, c in rows if c[0] == "missense"])', "                        for k, label, af, c in rows])"),
    (F, "a share dealt twice", "    return mine[mine.index % CHUNKS[complex_] == chunk]", "    return mine[mine.index % CHUNKS[complex_] <= chunk]"),
    (F, "the structure's residue not checked", '                if chains[c]["seq"][i] != v.ref_aa:', "                if False:"),
    (F, "an unresolved position said as scored", '    return subs, "scored" if subs else "not assessed: the position is not resolved"', '    return subs, "scored"'),
    (F, "the sign not turned", "    fw, f1 = -np.mean(fw, axis=0), -np.mean(f1, axis=0)", "    fw, f1 = np.mean(fw, axis=0), np.mean(f1, axis=0)"),
    (F, "the partner's term not taken off", "    return fw, fw - f1", "    return fw, fw"),
    (F, "a weak area passes", "        passed.append(bool(own >= AUC_MIN))", "        passed.append(True)"),
    (F, "a baseline need not be beaten", "            passed.append(bool(lo > 0))", "            passed.append(True)"),
    (F, "the baselines read the wrong way round", '    ok["base_blosum"], ok["base_n_aa"] = -ok.blosum, -ok.n_aa', '    ok["base_blosum"], ok["base_n_aa"] = ok.blosum, ok.n_aa'),
    (F, "one kind only still gets a verdict", "    if h.any() and not h.all():", "    if True:"),
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
