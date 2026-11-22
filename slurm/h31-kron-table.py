# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "numpy>=2.0",
#   "biopython>=1.85",
#   "safetensors>=0.4",
#   "kron-fm @ git+https://github.com/antonkulaga/kron-fm@58dfc5c39c52380082ad01ed34c8b65116237ffc",
# ]
# ///
"""H31: KronGenome's masked likelihoods along a host's mitochondrial genome, one table per (genome, context, arm).

Run with `uv run --script slurm/h31-kron-table.py` from the repository's root (the table's logic is the package's
`mitoport/kron_table.py`, read from `src/`; the model's code is kron-fm at the commit the model card names, which
cannot share an environment with the repository's other torch groups). The weights are the files of
just-dna-seq/KronGenome-F-1Mb at revision 3420801c763bda2b0d6a33c1e21d2d7276e3ec7c; each is used only with the
SHA-256 below.

  --weights DIR   config.json and model.safetensors
  --gb FILE       the GenBank records of `mitoage download cybrid`
  --out DIR       <accession>.<arm>.tsv are written here
  --accession A   (repeatable) the genomes; default the two hosts
  --device D      cuda or cpu
Arms: `mito` (context: genetic code 2), `nuclear` (context 0), `shuffled` (context 2, the bases shuffled within each
annotated feature and each stretch between features, seed 0: what composition alone gives).
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
from Bio import SeqIO
from safetensors.torch import load_file

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "mitoage_pipelines" / "mitoport"))
import kron_table as kt  # noqa: E402  (the one file, without the package's other imports)

SHA = {"config.json": "9f52afd382ebbae31f67d27aa13d4a3b241f65f2685e1bbd40d112eaf7beacab",
       "model.safetensors": "a294e18402ebde01508d39f720247bc25de1ebd5997d1302984875905a6aa4c2"}
HOSTS = ("NC_012920.1", "NC_005089.1")
ARMS = {"mito": 2, "nuclear": 0, "shuffled": 2}


def load(weights: Path, device: str):
    from kron_dna import genome as G
    for name, want in SHA.items():
        got = hashlib.sha256((weights / name).read_bytes()).hexdigest()
        if got != want:
            raise SystemExit(f"{weights / name}: SHA-256 {got}, not the pinned {want}")
    config = G.KronGenomeConfig(**json.loads((weights / "config.json").read_text()))
    model = G.KronGenome(config)
    model.load_state_dict(load_file(str(weights / "model.safetensors")), strict=True)
    return model.eval().to(device), 2 ** config.n_down


def predictor(model, context: int, device: str):
    def predict(rows: np.ndarray) -> np.ndarray:
        seq = torch.from_numpy(rows).to(device)
        ctx = torch.full((len(rows),), context, dtype=torch.long, device=device)
        with torch.no_grad():
            logits = model(seq, context=ctx)[0]                      # [B, L, channels]; the first four are A, C, G, T
        return torch.log_softmax(logits[..., :4].float(), dim=-1).cpu().numpy().astype(np.float64)
    return predict


def shuffled(rec, seed: int = 0) -> str:
    """The record's bases shuffled within each annotated gene-level feature and each stretch between them."""
    rng = np.random.default_rng(seed)
    n = len(rec.seq)
    cuts = {0, n}
    for f in rec.features:
        if f.type in ("CDS", "rRNA", "tRNA", "D-loop"):
            cuts |= {int(f.location.start), int(f.location.end)}
    cuts = sorted(c for c in cuts if 0 <= c <= n)
    s = np.array(list(str(rec.seq).upper()))
    for a, b in zip(cuts, cuts[1:]):
        s[a:b] = rng.permutation(s[a:b])
    return "".join(s)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--weights", type=Path, required=True)
    ap.add_argument("--gb", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--accession", action="append")
    ap.add_argument("--arm", action="append", choices=sorted(ARMS))
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--batch", type=int, default=8)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    torch.manual_seed(0)
    model, multiple = load(a.weights, a.device)
    print(f"model loaded on {a.device}: {sum(p.numel() for p in model.parameters())} parameters, lengths a multiple of {multiple}")
    records = {r.id: r for r in SeqIO.parse(str(a.gb), "genbank")}
    for acc in a.accession or HOSTS:
        rec = records[acc]
        for arm in a.arm or sorted(ARMS):
            seq = shuffled(rec) if arm == "shuffled" else str(rec.seq).upper()
            start = time.time()
            if a.device == "cuda":
                torch.cuda.reset_peak_memory_stats()
            tab = kt.table(kt.encode(seq), predictor(model, ARMS[arm], a.device), multiple=multiple, batch=a.batch)
            kt.write(a.out / f"{acc}.{arm}.tsv", seq, tab)
            own = tab[np.arange(len(seq)), np.clip(kt.encode(seq), 0, 3)][kt.encode(seq) < 4]
            mem = f", peak GPU memory {torch.cuda.max_memory_allocated() / 2**30:.2f} GiB" if a.device == "cuda" else ""
            print(f"{acc} {arm}: {len(seq)} places, mean log-probability of the genome's own base {own.mean():.4f}, "
                  f"own base the likeliest at {(tab.argmax(1)[kt.encode(seq) < 4] == kt.encode(seq)[kt.encode(seq) < 4]).mean():.4f} "
                  f"of the places, {time.time() - start:.0f} s{mem}", flush=True)


if __name__ == "__main__":
    main()
