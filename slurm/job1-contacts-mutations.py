"""Break each guard of H14's code (the contact maps of complexes I, III and V, and `mitoport readouts`) in turn and
check that a test fails.

Run by slurm/job1-contacts.sbatch after the tests have passed on the unbroken code. Each mutation replaces one line
of a source file, runs tests/test_complexes.py and tests/test_partners.py (stopping at the first failure), and puts
the file back byte for byte from a copy (git is not installed on the compute nodes). Prints one line per mutation,
"caught" or "SURVIVED", and exits non-zero if any survived, was not applied, or a file did not come back.
"""
import hashlib
import signal
import subprocess
import sys
from pathlib import Path

SRC = Path("src/mitoage_pipelines/mitoport")
CX, CT, BU, IF, PT, RO = (SRC / f for f in ("complexes.py", "contacts.py", "burial.py", "iface.py", "partners.py", "readouts.py"))
TESTS = ["tests/test_complexes.py", "tests/test_partners.py"]
MUTATIONS = [
    (CX, "an entity the table does not name is guessed",
     "if len(codes) != 1 or next(iter(codes)) not in cx.uniprot:", "if False:"),
    (CX, "a reference chain is not checked against its gene", "        if out.get(c) != g:", "        if False:"),
    (CT, "only the reference chains are mitochondrial partners",
     "is_mt = {c for c, g in gene_of_chain.items() if g in MT_GENES}", "is_mt = set(mt_chain)"),
    (BU, "chains.tsv marks only the reference chains as mitochondrial",
     '"is_mt": int(gene_of_chain[c] in MT_GENES)', '"is_mt": int(c in MTCHAIN)'),
    (BU, "species from the first alignment only",
     "species = sorted(set.intersection(*(set(alns[g]) for g in cx.mt)) & tree_tips",
     "species = sorted(set(alns[next(iter(cx.mt))]) & tree_tips"),
    (BU, "contacts counted against mitochondrial partners too",
     "nuc = ct[ct.partner_is_mt == 0].copy()", "nuc = ct.copy()"),
    (IF, "pair species from COX1, as the complex-IV-only code had them",
     "            sp = pair_species(aln, mt_aln, tips)", '            sp = sorted(set(aln) & set(mt_aln["COX1"]) & set(tips))'),
    (IF, "pair species from the first mitochondrial alignment only",
     "return sorted(set(aln) & set.intersection(*(set(a) for a in mt_aln.values())) & set(tips))",
     "return sorted(set(aln) & set(next(iter(mt_aln.values()))) & set(tips))"),
    (IF, "a gene in several chains keeps the last chain's interface", "            cols |= mine", "            cols = mine"),
    (IF, "a contact without a column is read as the last column",
     '    if len(unplaced):', '    if False:'),
    (IF, "every alignment outside the named groups is a cytosolic control",
     'groups["cytosolic control"] = sorted(g for g in nuc_aln if g in CYTOSOLIC_CONTROLS)',
     'groups["cytosolic control"] = sorted(g for g in nuc_aln if g not in {x for v in groups.values() for x in v})'),
    (IF, "subunits that touch nothing set the usual interface size",
     "sizes = [len(v) for v in iface_cols.values() if len(v) > 0]", "sizes = [len(v) for v in iface_cols.values()]"),
    (PT, "a contact without a control is refused again",
     "        if (not set(m.contact_aln_col) <= contacts or m.contact_aln_col.duplicated().any()",
     "        if (set(m.contact_aln_col) != contacts or m.contact_aln_col.duplicated().any()"),
    (RO, "one-sided group test", 'alternative="two-sided").pvalue)', 'alternative="greater").pvalue)'),
    (RO, "R1b on every subunit, with an interface or not", "touch = part[part.iface_cols.ge(MIN_IFACE)]", "touch = part"),
    (RO, "normal instead of t interval of the intercept", "h = stats.t.ppf(0.975, df) * se[0]", "h = 1.96 * se[0]"),
    (RO, "lambda fixed at zero in the intercept",
     'lam = float(optimize.minimize_scalar(nll, bounds=(0.0, 1.0), args=(X, y, C), method="bounded").x)', "lam = 0.0"),
    (RO, "Holm over every class of R3",
     'fam = res.status.eq("tested") & res.readout.isin(HOLM) & (res.readout.ne("R3") | res.set.eq(PRIMARY_CLASS))',
     'fam = res.status.eq("tested") & res.readout.isin(HOLM)'),
    (RO, "a contact paired with the wrong control",
     "d = np.array([by[(r.gene, r.contact_aln_col)] - by[(r.gene, r.control_aln_col)] for r in matched.itertuples()])",
     "d = np.array([by[(r.gene, r.contact_aln_col)] - by[(r.gene, r.contact_aln_col)] for r in matched.itertuples()])"),
    (RO, "too few proteins scored anyway", "    if len(a) < MIN_PROTEINS or len(b) < MIN_PROTEINS:", "    if len(a) < 1 or len(b) < 1:"),
]


def main() -> int:
    # Slurm's time limit (and kill) send SIGTERM: raise, so that `finally` puts the sources back
    signal.signal(signal.SIGTERM, lambda *_: sys.exit("terminated: sources put back"))
    originals = {f: f.read_bytes() for f in {m[0] for m in MUTATIONS}}
    digests = {f: hashlib.sha256(b).hexdigest() for f, b in originals.items()}
    bad = 0
    try:
        for path, name, old, new in MUTATIONS:
            text = originals[path].decode()
            if old not in text:
                print(f"NOT APPLIED (line not found): {name}"); bad += 1; continue
            path.write_text(text.replace(old, new, 1))
            r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", "--show-capture=no", *TESTS],
                               capture_output=True, text=True)
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
