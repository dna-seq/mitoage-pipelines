"""Break each guard of H21-H25's code (src/mitoage_pipelines/adaptations/) in turn and check that a test fails.

Run by slurm/adapt-*.sbatch after the tests have passed on the unbroken code. Each mutation replaces one line of the
source, runs tests/test_adapt.py and tests/test_hibernation.py (stopping at the first failure), and puts the file back
byte for byte from a copy (git is not installed on the compute nodes). Prints one line per mutation, "caught" or
"SURVIVED", and exits non-zero if any survived, was not applied, or a file did not come back.
"""
import hashlib
import signal
import subprocess
import sys
from pathlib import Path

A = Path("src/mitoage_pipelines/adaptations")
CL, SE, PA, CO, FI, RU, TR = (A / f"{m}.py" for m in ("clades", "selection", "pairs", "coupling", "fits", "run", "traits"))
TESTS = ["tests/test_adapt.py", "tests/test_hibernation.py"]
MUTATIONS = [
    (CL, "the inside of a trait clade is searched again", "        if not p:\n            for c in n.children: rec(c, False)", "        if True:\n            for c in n.children: rec(c, False)"),
    (CL, "a merged unary node loses its branch length", "            kids[0].length += n.length", "            pass"),
    (CL, "two pairs may share a branch", "        if edges & used:", "        if False:"),
    (CL, "the other sister may carry the trait", "        found = nearest_outside(c, lambda s: state[s] == 0)", "        found = nearest_outside(c, lambda s: True)"),
    (CL, "the outgroup may sit inside the pair's clade", "        og = nearest_outside(anc, lambda s: True)", "        og = (r, anc)"),
    (CL, "control clades of one draw may overlap", "            free = [z for z in p if not (z[4] & covered)]", "            free = list(p)"),
    (CL, "control clades are not matched on size", "        cands = sorted((z for z in zs if n / 2 <= z[1] <= 2 * n),", "        cands = sorted((z for z in zs),"),
    (CL, "every draw uses the same generator", "        rng = np.random.default_rng([0, k])", "        rng = np.random.default_rng(0)"),
    (SE, "a reference tip may carry the trait", "                           if t.is_tip and t.name not in inside and state[t.name] == 0 and t.name not in got)",
     "                           if t.is_tip and t.name not in inside and t.name not in got)"),
    (SE, "only the stems are test branches in RELAX", "            if s == cs or (not stems_only and s <= cs):", "            if s == cs:"),
    (SE, "stop codons are kept", '                if c in stops or ("-" in c and c != "---") or any(x not in "ACGT-" for x in c):',
     '                if ("-" in c and c != "---") or any(x not in "ACGT-" for x in c):'),
    (SE, "RELAX's K is read from the LRT", '"K": float(t["relaxation or intensification parameter"])', '"K": float(t["LRT"])'),
    (SE, "R1 is read with control draws missing", "    if len(ctl) < n_draws:", "    if False:"),
    (PA, "saturated pairs count as unsaturated", "max(vals[\"dS_trait\"], vals[\"dS_other\"]) > SATURATED_DS)", "max(vals[\"dS_trait\"], vals[\"dS_other\"]) > 99)"),
    (PA, "d_omega has the wrong sign", '            row["d_omega"] = row["d_dN"] - row["d_dS"]', '            row["d_omega"] = row["d_dS"] - row["d_dN"]'),
    (PA, "the minimum number of pairs is ignored", "    if len(x) < MIN_PAIRS:", "    if False:"),
    (CO, "a control's background keeps its own length", "        bg = (tot - scaled[g]) / (len(ctl) - 1) if g in ctl else tot / len(ctl)", "        bg = tot / len(ctl)"),
    (CO, "every branch of a clade counts as its stem", "        if stems_only and inside[id(n.parent)]:", "        if False:"),
    (CO, "C1 compares the foreground with itself", "    s = {g: logr[g][fg].median() - logr[g][~fg].median() for g in logr}", "    s = {g: logr[g][fg].median() - logr[g][fg].median() for g in logr}"),
    (CO, "an empty species set is run", "    if not sp:      # `prepare` kept no species", "    if False:      # `prepare` kept no species"),
    (FI, "L2 calls it hibernation's without the count halved", "    reading = (\"largely hibernation's\" if ratio < L2_LARGELY and halved", "    reading = (\"largely hibernation's\" if ratio < L2_LARGELY"),
    (FI, "L2 compares models on different species", "    if set(a.n) != set(b.n) or not a.index.equals(b.index):", "    if False:"),
    (RU, "R1 is read without beating the controls", "        if ctl < 19:", "        if ctl < 0:"),
    (RU, "C is read although the positive control failed", '    if r.readout in ("C1", "C2") and not positive_control_passed:', "    if False:"),
    (TR, "an intermediate depth counts as shallow", "    st = s.copy(); st[depth >= DEEP_M] = 1; st[depth < SHALLOW_M] = 0", "    st = s.copy(); st[depth >= DEEP_M] = 1; st[depth < DEEP_M] = 0"),
    (TR, "an empty AnaCat counts as not diadromous", "    st = s.copy(); st[ana.isin(DIADROMOUS)] = 1; st[ana.isin(OTHER_ANACAT)] = 0", "    st = s.copy(); st[ana.isin(DIADROMOUS)] = 1; st[~ana.isin(DIADROMOUS)] = 0"),
]


def main() -> int:
    # Slurm's time limit (and kill) send SIGTERM: raise, so that `finally` puts the source back
    signal.signal(signal.SIGTERM, lambda *_: sys.exit("terminated: sources put back"))
    originals = {f: f.read_bytes() for f in {m[0] for m in MUTATIONS}}
    digests = {f: hashlib.sha256(b).hexdigest() for f, b in originals.items()}
    bad = 0
    try:
        for path, name, old, new in MUTATIONS:
            text = originals[path].decode()
            if text.count(old) != 1:
                print(f"NOT APPLIED (line found {text.count(old)} times): {name}"); bad += 1; continue
            path.write_text(text.replace(old, new, 1))
            r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", "--show-capture=no", *TESTS],
                               capture_output=True, text=True)
            failed = [l.split(" - ")[0].removeprefix("FAILED ") for l in r.stdout.splitlines() if l.startswith(("FAILED", "ERROR"))]
            if r.returncode == 0:
                print(f"SURVIVED: {name}"); bad += 1
            else:
                print(f"caught: {name}: {', '.join(failed) or r.stdout.strip().splitlines()[-1]}", flush=True)
            path.write_bytes(originals[path])
    finally:
        for f, b in originals.items():
            f.write_bytes(b)
    back = all(hashlib.sha256(f.read_bytes()).hexdigest() == d for f, d in digests.items())
    print(f"{len(MUTATIONS) - bad} of {len(MUTATIONS)} mutations caught; sources restored: {back}")
    return 0 if bad == 0 and back else 1


if __name__ == "__main__":
    sys.exit(main())
