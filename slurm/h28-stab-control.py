"""H28's positive control, added after the result on a reviewer's advice: the pipeline's own input path
(`read_chains`, `structure`, `encode`, `predict`) on complexes of the model's benchmark, beside the predictions the
model's authors publish for the same mutations (their `baselines/StaB-ddG.csv`: `pred_1`) and the measured values
(`ddG`). Agreement says the model is called as its authors call it; it says nothing about the cybrids.

    uv run --frozen --group stab python slurm/h28-stab-control.py --benchmark <StaB-ddG.csv> --cif-dir <folder of mmCIF> \
        --weights <folder of the weights> --out <table.tsv> [--entries 1AO7,1A4Y,...] [--device cpu]

The mmCIF files are the Protein Data Bank's entries (the authors used SKEMPI's cleaned PDB files of the same
entries). A mutation whose residue number carries an insertion code, or whose stated residue is not the one the
entry holds at that number, is left out and counted.
"""
import argparse
import ast
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from mitoage_pipelines.mitoport import cybrid_stab as cs

p = argparse.ArgumentParser()
p.add_argument("--benchmark", type=Path, required=True)
p.add_argument("--cif-dir", type=Path, required=True)
p.add_argument("--weights", type=Path, required=True)
p.add_argument("--out", type=Path, required=True)
p.add_argument("--entries", default="1AO7,1A4Y,2G2U,1LFD,4G0N")
p.add_argument("--device", default="cpu")
a = p.parse_args()

file, digest = cs.WEIGHTS["stabddg"]
if cs.sha256(a.weights / file) != digest:
    raise SystemExit("the weights are not the pinned ones")
model = cs.load_model(a.weights / file, a.device)
bench = pd.read_csv(a.benchmark)
rows, left = [], {}
for entry in a.entries.split(","):
    sub = bench[bench["#Pdb"].str.startswith(entry + "_")]
    if sub.empty:
        raise SystemExit(f"{entry} is not in the benchmark")
    name = sub["#Pdb"].iloc[0]
    _, one, two = name.split("_")
    chains = cs.read_chains(a.cif_dir / f"{entry}.cif")
    index = {c["id"]: k for k, c in enumerate(chains)}
    missing = [x for x in one + two if x not in index]
    if missing:
        raise SystemExit(f"{entry}: no chain {missing}")
    groups = ([index[x] for x in one + two], [index[x] for x in one], [index[x] for x in two])
    structs, ats = zip(*(cs.structure(entry, chains, g) for g in groups))
    for _, r in sub.iterrows():
        subs, why = [{}, {}, {}], None
        for m in ast.literal_eval(r["Mutation"]):
            got = re.fullmatch(r"([A-Z])([A-Za-z])(-?\d+)([A-Z])", m)
            if not got:
                why = "an insertion code or another form"; break
            wt, ch, num, mut = got.group(1), got.group(2), int(got.group(3)), got.group(4)
            c = index.get(ch)
            where = [i for i, n in enumerate(chains[c]["num"]) if n == num] if c is not None else []
            if len(where) != 1 or chains[c]["seq"][where[0]] != wt:
                why = "the entry holds another residue there, or none"; break
            for k, at in enumerate(ats):
                if (c, where[0]) in at:
                    subs[k][at[(c, where[0])]] = mut
        if why:
            left[why] = left.get(why, 0) + 1
            continue
        seqs = tuple(cs.encode(s["seq"], x) for s, x in zip(structs, subs))
        out = cs.predict(model, structs, seqs)
        rows.append({"entry": name, "mutation": r["Mutation"], "measured": r["ddG"], "authors": r["pred_1"], "ours": out["ddG"], "ours_fold": out["fold"]})
        print(rows[-1], flush=True)
t = pd.DataFrame(rows)
a.out.parent.mkdir(parents=True, exist_ok=True)
t.to_csv(a.out, sep="\t", index=False, float_format="%.5f")
print("left out:", left or "none")
for name, g in list(t.groupby("entry")) + [("all", t)]:
    print(f"{name}: n {len(g)}; ours against the authors' prediction: Pearson {pearsonr(g.ours, g.authors)[0]:.3f}, "
          f"mean absolute difference {np.abs(g.ours - g.authors).mean():.3f}; against the measured: Spearman ours "
          f"{spearmanr(g.ours, g.measured)[0]:.3f}, the authors' {spearmanr(g.authors, g.measured)[0]:.3f}")
