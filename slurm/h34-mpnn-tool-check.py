"""H34, part B: our per-position reading (`variant_mpnn.read`, on the structure as `variant_mpnn.features` gives it)
against the LigandMPNN repository's own `score.py`, end to end from a file, for each of the four arms.

Two structures: the repository's example 1BC8, and two chains of complex IV (COX2 and the smallest nuclear chain of
H29's system around it, from 5Z62, written to a PDB file). For each, `score.py --single_aa_score 1 --use_sequence 0`
(the reading that sees the other residues) with three draws; then our reading of the same file at every place
where the repository decoded it last, with the order of the other places its output reports. Prints, per structure
and arm, the places compared, the largest difference of a log-probability and the native residue's mean
log-probability, and exits non-zero if a sequence or mask differs or a difference exceeds TOLERANCE.

Needs ProDy for the repository's reader, which the project does not carry:
    UV_PROJECT_ENVIRONMENT=.venv-stab uv run --frozen --group stab --group dev --with prody python slurm/h34-mpnn-tool-check.py <folder>
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import torch

from mitoage_pipelines import paths
from mitoage_pipelines.mitoport import cybrid
from mitoage_pipelines.mitoport import cybrid_stab as cs
from mitoage_pipelines.mitoport import variant_mpnn as vm

TOLERANCE = 1e-4
SEED = 1                 # score.py takes a seed of 0 as none and draws one
FLAG = {"protein_mpnn": "--checkpoint_protein_mpnn", "soluble_mpnn": "--checkpoint_soluble_mpnn",
        "global_label_membrane_mpnn": "--checkpoint_global_label_membrane_mpnn"}


def complex_iv_pair(folder: Path) -> Path:
    from Bio.PDB import PDBIO, MMCIFParser, Select
    from mitoage_pipelines.mitoport.complexes import COMPLEXES
    hp = cybrid.species_proteins(cybrid.read_proteins(paths.CYBRID_GB / "records.gb"), cybrid.HUMAN)
    genes = tuple(COMPLEXES["IV"].mt)
    chains, mt = cs.system(cs.read_chains(paths.PDB / "5Z62.cif"), {g: hp[g][0] for g in genes}, {g: p for g, (p, _, _) in hp.items() if g not in genes})
    cox2 = next(k for k, (g, _) in mt.items() if g == "COX2")
    small = min((k for k in range(len(chains)) if k not in mt), key=lambda k: len(chains[k]["seq"]))
    ids = {chains[cox2]["id"], chains[small]["id"]}
    if any(len(i) != 1 for i in ids):
        raise SystemExit(f"chain ids {ids} do not fit a PDB file")

    class Two(Select):
        def accept_chain(self, chain):
            return chain.id in ids

        def accept_residue(self, res):        # protein residues only, as both readers take them
            return res.id[0] == " "

    s = MMCIFParser(QUIET=True).get_structure("5Z62", str(paths.PDB / "5Z62.cif"))
    io = PDBIO(); io.set_structure(next(iter(s))); f = folder / "5Z62_COX2_pair.pdb"
    io.save(str(f), Two())
    return f


def tool_output(pdb: Path, arm: str, folder: Path) -> dict:
    file, model_type, label = vm.ARMS[arm]
    out = folder / f"{pdb.stem}.{arm}"
    cmd = [sys.executable, "score.py", "--model_type", model_type, FLAG[model_type], str((paths.LIGANDMPNN / "model_params" / file).resolve()),
           "--seed", str(SEED), "--single_aa_score", "1", "--use_sequence", "0", "--batch_size", "1", "--number_of_batches", "3",
           "--pdb_path", str(pdb.resolve()), "--out_folder", str(out.resolve())]
    if label is not None:
        cmd += ["--global_transmembrane_label", str(label)]
    subprocess.run(cmd, cwd=paths.LIGANDMPNN, check=True, stdout=subprocess.DEVNULL)
    return torch.load(out / f"{pdb.stem}.pt", weights_only=False)


def main(folder: Path) -> int:
    folder.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(int(__import__("os").environ.get("OMP_NUM_THREADS", "4")))
    bad = 0
    for pdb in (paths.LIGANDMPNN / "inputs" / "1BC8.pdb", complex_iv_pair(folder)):
        chains = cs.read_chains(pdb)
        for arm in vm.ARMS:
            t = tool_output(pdb, arm, folder)
            model = vm.load(paths.LIGANDMPNN, arm)
            feats, _ = vm.features(chains, list(range(len(chains))), vm.ARMS[arm][2])
            L = feats["S"].shape[1]
            if L != len(t["native_sequence"]) or not np.array_equal(feats["S"][0].numpy(), t["native_sequence"]) or not np.array_equal(feats["mask"][0].numpy(), t["mask"]):
                print(f"{pdb.name} {arm}: the structure differs from the repository's reading ({L} places, the repository's {len(t['native_sequence'])})"); bad += 1
                continue
            diffs, native, n = [], [], 0
            for b in range(t["log_probs"].shape[0]):
                for idx in range(L):
                    row = t["decoding_order"][b, idx].astype(int)
                    if row[-1] != idx:
                        continue
                    noise = torch.empty(L); noise[torch.from_numpy(row)] = torch.arange(L, dtype=torch.float32)
                    mine = vm.read(model, feats, [idx], noise)[0]
                    diffs.append(np.abs(mine - t["log_probs"][b, idx]).max()); native.append(mine[t["native_sequence"][idx]]); n += 1
            worst = max(diffs)
            bad += worst > TOLERANCE
            print(f"{pdb.name} {arm}: {n} of {3 * L} places compared (the rest the repository did not decode last), largest difference "
                  f"{worst:.2e}; native residue's mean log-probability {np.mean(native):.3f}; {'agrees' if worst <= TOLERANCE else 'DIFFERS'}", flush=True)
    print("tool check:", "every reading agrees" if bad == 0 else f"{bad} readings differ")
    return 0 if bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
