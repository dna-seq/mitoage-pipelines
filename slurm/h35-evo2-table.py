# /// script
# requires-python = "==3.12.*"
# dependencies = [
#   "numpy>=2.0",
#   "biopython>=1.85",
#   "torch==2.7.1",
#   "evo2==0.6.0",
#   "flash-attn @ https://github.com/Dao-AILab/flash-attention/releases/download/v2.8.0.post2/flash_attn-2.8.0.post2+cu12torch2.7cxx11abiTRUE-cp312-cp312-linux_x86_64.whl",
# ]
# ///
"""H35, part 1: Evo 2's likelihoods along a host's mitochondrial genome, one table per (genome, arm).

Run with `uv run --script slurm/h35-evo2-table.py` from the repository's root, on a GPU (the model's code does not
import without one). The table's logic is the package's `mitoport/evo_table.py`, read from `src/`. The weights are
`evo2_7b.pt` of arcinstitute/evo2_7b at revision bda0089f92582d5baabf0f22d9fc85f3588f6b58, used only with the
SHA-256 below.

  --weights FILE  evo2_7b.pt
  --gb FILE       the GenBank records of `mitoage download cybrid`
  --out DIR       <accession>.<arm>.tsv are written here (arms: evo2, evo2plus, evo2minus)
  --accession A   (repeatable) the genomes; default the two hosts
"""
import argparse
import hashlib
import sys
import time
from pathlib import Path

import numpy as np
import torch
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "mitoage_pipelines" / "mitoport"))
import evo_table as et  # noqa: E402  (the one file and kron_table beside it, without the package's other imports)
import kron_table as kt  # noqa: E402

SHA = "c66645929dc1b9c631f5be656da8726f38946315dc9167000a615dd626fcecf4"
HOSTS = ("NC_012920.1", "NC_005089.1")
ARMS = {"evo2": "both", "evo2plus": "plus", "evo2minus": "minus"}
TOKENS = [ord(b) for b in et.BASES]          # the model's tokens are the bytes of the letters


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def predictor(model):
    def predict(seq: str) -> np.ndarray:
        ids = torch.tensor(model.tokenizer.tokenize(seq), dtype=torch.int).unsqueeze(0).to("cuda:0")
        with torch.inference_mode():
            outputs, _ = model(ids)
        logits = outputs[0][0] if outputs[0].dim() == 3 else outputs[0]          # [L, vocabulary]
        if logits.shape[0] != len(seq):
            raise SystemExit(f"the model returned {tuple(logits.shape)} for {len(seq)} bases")
        lp = torch.log_softmax(logits[:, TOKENS].float(), dim=-1).cpu().numpy().astype(np.float64)
        out = np.zeros((len(seq), 4))
        out[1:] = lp[:-1]                      # the logits at place i give the base at place i + 1
        return out
    return predict


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--weights", type=Path, required=True)
    ap.add_argument("--gb", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--accession", action="append")
    a = ap.parse_args()
    got = sha256(a.weights)
    if got != SHA:
        raise SystemExit(f"{a.weights}: SHA-256 {got}, not the pinned {SHA}")
    a.out.mkdir(parents=True, exist_ok=True)
    torch.manual_seed(0)
    from evo2 import Evo2
    model = Evo2("evo2_7b", local_path=str(a.weights))
    print(f"model loaded on {torch.cuda.get_device_name(0)}")
    records = {r.id: r for r in SeqIO.parse(str(a.gb), "genbank")}
    for acc in a.accession or HOSTS:
        seq = str(records[acc].seq).upper()
        start = time.time()
        torch.cuda.reset_peak_memory_stats()
        tabs = et.tables(seq, predictor(model))
        codes = kt.encode(seq)
        known = codes < 4
        for arm, key in ARMS.items():
            tab = tabs[key]
            kt.write(a.out / f"{acc}.{arm}.tsv", seq, tab)
            own = tab[np.arange(len(seq)), np.clip(codes, 0, 3)][known]
            top = (np.nanargmax(np.where(np.isnan(tab), -np.inf, tab), axis=1) == codes)[known].mean()
            print(f"{acc} {arm}: {len(seq)} places, mean log-probability of the genome's own base {own.mean():.4f}, "
                  f"the likeliest base is the genome's own at {100 * top:.1f}% of the places")
        print(f"{acc}: {time.time() - start:.0f} s, peak GPU memory {torch.cuda.max_memory_allocated() / 2**30:.2f} GiB")


if __name__ == "__main__":
    main()
