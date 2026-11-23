"""Break each guard of H34's part B (`mitoage mitoport variant-mpnn`, its reading, its collection and its readout)
in turn and check that a test of tests/test_variant_mpnn.py fails. Each mutation replaces one piece of the source,
runs the tests, and the file is put back byte for byte. Prints one line per mutation, "caught" or "SURVIVED", and
exits non-zero if any survived, was not applied, or the file did not come back.

    UV_PROJECT_ENVIRONMENT=.venv-stab uv run --frozen --group stab --group dev python slurm/h34-mpnn-mutations.py
"""
import hashlib
import subprocess
import sys
from pathlib import Path

F = Path("src/mitoage_pipelines/mitoport/variant_mpnn.py")
TESTS = ["tests/test_variant_mpnn.py"]
MUTATIONS = [
    ("the scored place decoded first", "            rank[torch.arange(b), p] = L ", "            rank[torch.arange(b), p] = -1 "),
    ("the others not in the seed's order", "        base[torch.argsort(noise.to(E_idx.device), stable=True)] = torch.arange(L, device=E_idx.device)",
     "        base[:] = torch.arange(L, device=E_idx.device)"),
    ("the scored place sees its own residue", ".view(b, L, K) < rank[:, :, None])", ".view(b, L, K) <= rank[:, :, None])"),
    ("a residue without its backbone not masked", "mask.append(ok.all(axis=1))", "mask.append(np.ones(len(bb), bool))"),
    ("a missing atom not set to zero", "X.append(np.where(ok[:, :, None], bb, 0.0))", "X.append(bb)"),
    ("an insertion code not renumbered", '"R_idx": torch.tensor(renumber(nums), dtype', '"R_idx": torch.tensor(nums, dtype'),
    ("every chain given one label", "chain += [k] * len(bb)", "chain += [0] * len(bb)"),
    ("the membrane label ignored", "torch.full((1, L), int(label), dtype", "torch.full((1, L), 1, dtype"),
    ("a checkpoint not checked", "    if not w.is_file() or sha256(w) != CHECKPOINTS[file]:", "    if not w.is_file():"),
    ("the repository's code not checked", "    if not f.is_file() or sha256(f) != MODEL_UTILS_SHA256:", "    if not f.is_file():"),
    ("the variants of other complexes kept", "    v = v[v.gene.isin(list(genes))].reset_index(drop=True)\n", ""),
    ("one seed's order for every seed", "lp = reader(model, feats, positions, order_noise(L, seed), batch)", "lp = reader(model, feats, positions, order_noise(L, 0), batch)"),
    ("the residues read at the wrong place", "                    x = lp[col[at[(c, i)]]]", "                    x = lp[0]"),
    ("the sign turned", '"lp_ref": a, "lp_alt": b, "score": a - b}', '"lp_ref": a, "lp_alt": b, "score": b - a}'),
    ("the variants not placed left out", "        if status != \"scored\":\n            rows.append", "        if False:\n            rows.append"),
    ("the first copy for all copies", 'per_seed = sc.groupby(["arm", "id", "seed"]).score.mean()', 'per_seed = sc.groupby(["arm", "id", "seed"]).score.first()'),
    ("a change of the start kept", "    if (long.protein_pos == 1).any():", "    if False:"),
    ("a weak area passes", "verdict = \"counts\" if cybrid.auc(ok[MAIN], ok.harmful) >= AUC_MIN and lo > 0", "verdict = \"counts\" if lo > 0"),
    ("the baseline need not be beaten", "cybrid.auc(ok[MAIN], ok.harmful) >= AUC_MIN and lo > 0", "cybrid.auc(ok[MAIN], ok.harmful) >= AUC_MIN"),
    ("the baseline read the wrong way round", 'ok["base_blosum"], ok[BASELINE] = -ok.blosum, -ok.n_aa', 'ok["base_blosum"], ok[BASELINE] = -ok.blosum, ok.n_aa'),
    ("the soluble weights' difference the wrong way round", 'difference(ok, MAIN, "smpnn", "all"', 'difference(ok, "smpnn", MAIN, "all"'),
    ("one kind only still gets a verdict", "    if h.any() and not h.all():\n        area(ok, MAIN", "    if True:\n        area(ok, MAIN"),
]


def main() -> int:
    original = F.read_bytes()
    digest = hashlib.sha256(original).hexdigest()
    bad = 0
    try:
        for name, old, new in MUTATIONS:
            text = original.decode()
            if text.count(old) != 1:
                print(f"NOT APPLIED (the piece is not there exactly once): {name}"); bad += 1; continue
            F.write_text(text.replace(old, new))
            r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-x", "-p", "no:cacheprovider", "--show-capture=no", *TESTS],
                               capture_output=True, text=True)
            failed = [line.split(" - ")[0].removeprefix("FAILED ") for line in r.stdout.splitlines() if line.startswith(("FAILED", "ERROR"))]
            if r.returncode == 0:
                print(f"SURVIVED: {name}", flush=True); bad += 1
            else:
                print(f"caught: {name}: {', '.join(failed) or r.stdout.strip().splitlines()[-1]}", flush=True)
            F.write_bytes(original)
    finally:
        F.write_bytes(original)
    back = hashlib.sha256(F.read_bytes()).hexdigest() == digest
    print(f"{len(MUTATIONS) - bad} of {len(MUTATIONS)} mutations caught; sources restored: {back}")
    return 0 if bad == 0 and back else 1


if __name__ == "__main__":
    sys.exit(main())
