"""H34, part B: does ProteinMPNN's membrane checkpoint tell harmful from harmless substitutions within a
mitochondrial protein better than the site's variability among mammals?

The model is the membrane variant of ProteinMPNN in the LigandMPNN repository (github.com/dauparas/LigandMPNN at
COMMIT; its weights by the repository's `get_model_params.sh`). Each variant of H29 (`variant_stab.variants`, the set
corrected on 2026-11-22: no changes of a start codon) is read on the human structure and system of H28 and H29
(`cybrid_stab.system`, every copy of the chain) as

    score = log p(reference residue) - log p(new residue)

at the variant's position, decoded last, after every other residue of the system in the order the seed gives, so
that it sees the identity of every other residue and not its own: the repository's `single_aa_score` with
`--use_sequence 0` (which, against its name, is the reading that sees the sequence; H34's note of 2026-11-22). Higher
is the prediction of harm. The repository's own function builds an L by L by L product for each of the L positions;
`read` computes the same numbers for the variants' positions only, with the model's own encoder and decoder layers,
and the tests hold it equal to the repository's function on a real structure.

Fixed in the entry before any score: the checkpoint `global_label_membrane_mpnn_v_48_020` with the global label 1
(transmembrane) is the score; the controls are the same checkpoint with label 0, `proteinmpnn_v_48_020` and
`solublempnn_v_48_020`; seeds SEEDS, the score the mean over them and over the copies of the chain. The rule is
H29's against the site's variability alone: the area under the ROC curve for harmful against harmless is at least
AUC_MIN and above that of `n_aa` on the same variants, the 95% bootstrap interval of the difference above 0.

The structure is given to the model as the repository's reader gives it: backbone N, CA, C, O as deposited (zeros
and masked where one of the four is missing), the residues' own numbers (renumbered past an insertion code as its
`featurize` does), one label per chain. `variant-mpnn` runs one complex; `variant-mpnn-collect` and
`variant-mpnn-readout` follow. Needs `uv sync --group stab` and the repository with its weights in
data/input/ligandmpnn.
"""
from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import cybrid
from . import cybrid_stab as cs
from . import variant_stab as vs
from .cybrid_energy import STRUCTURES

TOOL = "https://github.com/dauparas/LigandMPNN"
COMMIT = "26ec57ac976ade5379920dbd43c7f97a91cf82de"
MODEL_UTILS_SHA256 = "5ad8f6992f83cdd1c97103e19c80833afff2bff201662a23bcb06eae669f0f95"
CHECKPOINTS = {          # file in the repository's model_params -> SHA-256
    "global_label_membrane_mpnn_v_48_020.pt": "89ef0abddbfb956c4c7c02dcba6523c0b5152733055f7509292365eaecc38b21",
    "proteinmpnn_v_48_020.pt": "c9cb4a671d79604111231f8dbfc7c590e06f1197453b7a6854ac6661a642f5bd",
    "solublempnn_v_48_020.pt": "7af52d090172c230c7f0e9d21e02203f6b3a38b16db58d3c7a3960e0a9a6e31a",
}
ARMS = {                 # name -> (checkpoint, the repository's model type, global membrane label or None)
    "memb": ("global_label_membrane_mpnn_v_48_020.pt", "global_label_membrane_mpnn", 1),
    "memb0": ("global_label_membrane_mpnn_v_48_020.pt", "global_label_membrane_mpnn", 0),
    "pmpnn": ("proteinmpnn_v_48_020.pt", "protein_mpnn", None),
    "smpnn": ("solublempnn_v_48_020.pt", "soluble_mpnn", None),
}
MAIN = "memb"
SEEDS = (0, 1, 2)
BATCH = 4                # positions decoded in one call; it changes the memory, not the numbers
AUC_MIN = 0.70
BASELINE = "base_n_aa"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def tool_module(tool: Path):
    """The repository's model_utils, refused unless it is the file of the pinned commit."""
    f = tool / "model_utils.py"
    if not f.is_file() or sha256(f) != MODEL_UTILS_SHA256:
        raise SystemExit(f"{f} is missing or is not the file of {TOOL} at {COMMIT}")
    if "ligandmpnn_model_utils" not in sys.modules:
        spec = importlib.util.spec_from_file_location("ligandmpnn_model_utils", f)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        sys.modules["ligandmpnn_model_utils"] = mod
    return sys.modules["ligandmpnn_model_utils"]


def load(tool: Path, arm: str, device: str = "cpu"):
    """The model of an arm, built as the repository's score.py builds it, in eval mode. Refuses a checkpoint whose
    SHA-256 is not the pinned one."""
    import torch
    file, model_type, _ = ARMS[arm]
    w = tool / "model_params" / file
    if not w.is_file() or sha256(w) != CHECKPOINTS[file]:
        raise SystemExit(f"{w} is missing or is not the pinned checkpoint")
    mu = tool_module(tool)
    ck = torch.load(w, map_location=device, weights_only=True)
    model = mu.ProteinMPNN(node_features=128, edge_features=128, hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3,
                           k_neighbors=ck["num_edges"], device=device, atom_context_num=1, model_type=model_type,
                           ligand_mpnn_use_side_chain_context=0)
    model.load_state_dict(ck["model_state_dict"])
    return model.to(device).eval()


def renumber(nums: list[int]) -> list[int]:
    """The residue numbers as the repository's `featurize` gives them: a number equal to the one before (an insertion
    code) moves it and every later number on by one."""
    out, count, prev = [], 0, -100000
    for n in nums:
        if n == prev:
            count += 1
        out.append(n + count)
        prev = n
    return out


def features(chains: list[dict], which: list[int], label: int | None, device: str = "cpu") -> tuple[dict, dict[tuple[int, int], int]]:
    """(the model's input for the chains `which`, as the repository's reader and `featurize` make it, with a batch
    of one; {(chain index, residue index in the chain): index in the input}). A residue without all four backbone
    atoms has zeros there and is masked. `label` is the global membrane label of every residue (1 transmembrane, 0
    soluble), or None for a checkpoint without one."""
    import torch
    X, mask, nums, chain, S, at = [], [], [], [], [], {}
    for k, c in enumerate(which):
        bb = chains[c]["bb"]
        ok = np.isfinite(bb).all(axis=2)
        at.update({(c, i): len(S) + i for i in range(len(bb))})
        X.append(np.where(ok[:, :, None], bb, 0.0)); mask.append(ok.all(axis=1))
        nums += list(chains[c]["num"]); chain += [k] * len(bb)
        S += [cs.ALPHABET.index(a) if a in cs.ALPHABET else cs.ALPHABET.index("X") for a in chains[c]["seq"]]
    L = len(S)
    f = {"X": torch.tensor(np.concatenate(X), dtype=torch.float32)[None], "mask": torch.tensor(np.concatenate(mask), dtype=torch.int32)[None],
         "R_idx": torch.tensor(renumber(nums), dtype=torch.int64)[None], "chain_labels": torch.tensor(chain, dtype=torch.int64)[None],
         "S": torch.tensor(S, dtype=torch.int64)[None], "chain_mask": torch.ones(1, L, dtype=torch.int32)}
    if label is not None:
        f["membrane_per_residue_labels"] = torch.full((1, L), int(label), dtype=torch.int64)
    return {k: v.to(device) for k, v in f.items()}, at


def order_noise(L: int, seed: int):
    """The draw that orders the decoding for one seed: |N(0, 1)| for each place, drawn on the CPU (the same numbers
    on any device), as the repository draws `randn` for a batch of one."""
    import torch
    g = torch.Generator().manual_seed(seed)
    return torch.randn(L, generator=g).abs()


def read(model, feats: dict, positions: list[int], noise, batch: int = BATCH) -> np.ndarray:
    """[len(positions), 21] log-probabilities of the model's alphabet at each position, the position decoded last
    and every other place before it in the order of `noise` (smallest first), so that it sees every other residue's
    identity and not its own. The repository's `single_aa_score` with `use_sequence` False, at these positions only."""
    import torch
    mu = sys.modules[type(model).__module__]
    out = []
    with torch.no_grad():
        h_V0, h_E, E_idx = model.encode(feats)
        S, mask = feats["S"], feats["mask"]
        L, K = E_idx.shape[1], E_idx.shape[2]
        h_S = model.W_s(S)
        h_ES = mu.cat_neighbors_nodes(h_S, h_E, E_idx)
        h_EXV = mu.cat_neighbors_nodes(h_V0, mu.cat_neighbors_nodes(torch.zeros_like(h_S), h_E, E_idx), E_idx)
        base = torch.empty(L, dtype=torch.int64, device=E_idx.device)
        base[torch.argsort(noise.to(E_idx.device), stable=True)] = torch.arange(L, device=E_idx.device)
        mask_1D = mask.view(1, L, 1, 1)
        for b0 in range(0, len(positions), batch):
            p = torch.tensor(positions[b0:b0 + batch], dtype=torch.int64, device=E_idx.device)
            b = len(p)
            rank = base.repeat(b, 1)
            rank[torch.arange(b), p] = L                      # each scored position after every other place
            E = E_idx.expand(b, -1, -1)
            earlier = (torch.gather(rank, 1, E.reshape(b, -1)).view(b, L, K) < rank[:, :, None]).float().unsqueeze(-1)
            mask_bw, mask_fw = mask_1D * earlier, mask_1D * (1.0 - earlier)
            fw = mask_fw * h_EXV
            h_V, hES = h_V0.expand(b, -1, -1), h_ES.expand(b, -1, -1, -1)
            for layer in model.decoder_layers:
                h_ESV = mask_bw * mu.cat_neighbors_nodes(h_V, hES, E) + fw
                h_V = layer(h_V, h_ESV, mask)
            logits = model.W_out(h_V[torch.arange(b), p])
            out.append(torch.nn.functional.log_softmax(logits, dim=-1).cpu().numpy())
    return np.concatenate(out) if out else np.zeros((0, len(cs.ALPHABET)))


def run(complex_: str, gb: Path, pdb_dir: Path, source: Path, tool: Path, out: Path, device: str = "cpu", batch: int = BATCH,
        seeds: tuple[int, ...] = SEEDS, arms: tuple[str, ...] = tuple(ARMS), model_for=None, reader=read) -> pd.DataFrame:
    """Every variant of the complex's genes, each arm and seed: one row per variant, arm, seed and copy of its chain
    with the log-probabilities of the reference and the new residue there; a variant that cannot be placed has one
    row saying why. `model_for(arm)` and `reader` replace the model and its reading (for tests)."""
    from .complexes import COMPLEXES
    v, _ = vs.variants(gb, source / "clinvar_mt.json", source / "gnomad.chrM.tsv")
    genes = tuple(COMPLEXES[complex_].mt)
    v = v[v.gene.isin(list(genes))].reset_index(drop=True)
    pid = STRUCTURES[cybrid.HUMAN][complex_]
    cif = pdb_dir / f"{pid}.cif"
    if not cif.is_file():
        raise SystemExit(f"{cif} is missing: `mitoage download pdb {pid}`")
    hp = cybrid.species_proteins(cybrid.read_proteins(gb), cybrid.HUMAN)
    chains, mt_chain = cs.system(cs.read_chains(cif), {g: hp[g][0] for g in genes}, {g: p for g, (p, _, _) in hp.items() if g not in genes})
    which = list(range(len(chains)))
    placed = [vs.place(r, chains, mt_chain) for r in v.itertuples()]
    todo = [k for k, (_, status) in enumerate(placed) if status == "scored"]
    _, at = features(chains, which, None)
    positions = sorted({at[ci] for k in todo for ci in placed[k][0]})
    col = {p: j for j, p in enumerate(positions)}
    L = len(at)
    print(f"{complex_} ({pid}): {len(chains)} chains, {L} residues; {len(todo)} of {len(v)} variants placed, {len(positions)} positions, on {device}", flush=True)
    rows = []
    for arm in arms:
        model = model_for(arm) if model_for else load(tool, arm, device)
        feats, _ = features(chains, which, ARMS[arm][2], device)
        for seed in seeds:
            lp = reader(model, feats, positions, order_noise(L, seed), batch)
            for k in todo:
                r = v.iloc[k]
                for (c, i), alt in sorted(placed[k][0].items()):
                    x = lp[col[at[(c, i)]]]
                    a, b = float(x[cs.ALPHABET.index(r.ref_aa)]), float(x[cs.ALPHABET.index(alt)])
                    rows.append({**r.to_dict(), "complex": complex_, "arm": arm, "seed": seed, "chain": c, "copies": len(placed[k][0]),
                                 "status": "scored", "lp_ref": a, "lp_alt": b, "score": a - b})
            print(f"{complex_} {arm} seed {seed}: {len(positions)} positions read", flush=True)
    for k, (_, status) in enumerate(placed):
        if status != "scored":
            rows.append({**v.iloc[k].to_dict(), "complex": complex_, "arm": "", "seed": -1, "chain": -1, "copies": 0, "status": status,
                         "lp_ref": np.nan, "lp_alt": np.nan, "score": np.nan})
    res = pd.DataFrame(rows, columns=[*v.columns, "complex", "arm", "seed", "chain", "copies", "status", "lp_ref", "lp_alt", "score"])
    d = out / "mpnn"; d.mkdir(parents=True, exist_ok=True)
    f = d / f"{complex_}.tsv"
    res.to_csv(f, sep="\t", index=False, float_format="%.6f")
    print("wrote", f, f"({len(todo)} scored, {len(v) - len(todo)} not assessed)")
    return res


def scores(long: pd.DataFrame) -> pd.DataFrame:
    """One row per variant from the tasks' rows: each arm's score, the mean over its seeds of the mean over the
    copies, and for the score arm the spread over the seeds and the number of seeds. Refuses a change in a first
    codon, which is a change of the start and no missense variant (the correction of 2026-11-22)."""
    if (long.protein_pos == 1).any():
        raise SystemExit("a change of a start codon is in the tables: they predate the correction of 2026-11-22")
    info =long.drop_duplicates("id").set_index("id")[["gene", "protein_pos", "ref_aa", "alt_aa", "label", "pos", "ref", "alt", "af_hom", "complex", "copies", "status"]]
    sc = long[long.status == "scored"]
    per_seed = sc.groupby(["arm", "id", "seed"]).score.mean()
    cols = {}
    for arm in ARMS:
        if arm not in per_seed.index.get_level_values(0):
            continue
        s = per_seed.loc[arm].groupby("id")
        cols[arm] = s.mean()
        if arm == MAIN:
            cols[f"{arm}_sd"], cols["seeds"] = s.std(ddof=1), s.size()
    return info.join(pd.DataFrame(cols)).reset_index()


def collect(out: Path, data: Path, sites_dir) -> pd.DataFrame:
    """mpnn_scores.tsv: one row per variant with each arm's score, BLOSUM62, the site's variability among the build's
    mammals (`n_aa`) and whether the site is a mito-nuclear contact, as H29's collection. A complex without a table is
    said, never passed over; its variants are absent, not zero."""
    from Bio.Align import substitution_matrices
    blosum = substitution_matrices.load("BLOSUM62")
    tabs, missing = [], []
    for c in STRUCTURES[cybrid.HUMAN]:
        f = out / "mpnn" / f"{c}.tsv"
        (tabs.append(pd.read_csv(f, sep="\t", keep_default_na=False, na_values=[""])) if f.is_file() else missing.append(c))
    if missing:
        print(f"NOT ASSESSED: {len(missing)} of {len(STRUCTURES[cybrid.HUMAN])} complexes have no table: " + ", ".join(missing))
    if not tabs:
        raise SystemExit("no table to collect")
    res = scores(pd.concat(tabs, ignore_index=True))
    res["blosum"] = [float(blosum[a][b]) for a, b in zip(res.ref_aa, res.alt_aa)]
    var = {g: vs.site_variability(data, g) for g in sorted(set(res.gene))}
    res["n_aa"] = [var[g].get(p, np.nan) for g, p in zip(res.gene, res.protein_pos)]
    contact = set()
    for c in sorted(set(res["complex"])):
        s = pd.read_csv(sites_dir(c) / "sites.tsv", sep="\t", dtype={"chain": str})
        contact |= {(r.gene, r.human_pos) for r in s[s.group == "contact"].itertuples()}
    res["contact"] = [(g, p) in contact for g, p in zip(res.gene, res.protein_pos)]
    res = res.sort_values(["gene", "protein_pos", "alt_aa", "pos"]).reset_index(drop=True)
    res.to_csv(out / "mpnn_scores.tsv", sep="\t", index=False, float_format="%.5f")
    print(res.groupby(["label", "status"]).size().to_string()); print("wrote", out / "mpnn_scores.tsv")
    return res


def readout(out: Path) -> pd.DataFrame:
    """The rule on mpnn_scores.tsv, and beside it, deciding nothing: each control's area, BLOSUM62's, the score minus
    the soluble weights with its interval, the split by contact, and each complex alone. Writes h34_mpnn_summary.tsv."""
    f = out / "mpnn_scores.tsv"
    if not f.is_file():
        raise SystemExit(f"{f} is missing: `mitoage mitoport variant-mpnn-collect`")
    v = pd.read_csv(f, sep="\t")
    ok = v[(v.status == "scored") & v[list(ARMS)].notna().all(axis=1) & v.n_aa.notna()].copy()
    ok["harmful"] = ok.label == "harmful"
    ok["base_blosum"], ok[BASELINE] = -ok.blosum, -ok.n_aa        # higher is the prediction of harm, as for the scores
    rows = [{"what": "variants", "group": "all", "score": "", "n_harmful": int((v.label == "harmful").sum()), "n_harmless": int((v.label == "harmless").sum()),
             "value": np.nan, "reads": f"{len(ok)} scored by every arm with the baseline; {len(v) - len(ok)} not assessed"}]

    def area(d, s, group, what="area (reported, claimed as nothing)"):
        h = d.harmful
        rows.append({"what": what, "group": group, "score": s, "n_harmful": int(h.sum()), "n_harmless": int((~h).sum()),
                     "value": cybrid.auc(d[s], h) if h.any() and not h.all() else np.nan, "reads": "area under the ROC curve, harmful against harmless"})

    def difference(d, s, b, group, what):
        h = d.harmful.to_numpy()
        a, ab, lo, hi = vs.compare(h, d[s].to_numpy(float), d[b].to_numpy(float))
        rows.append({"what": what, "group": group, "score": f"{s} - {b}", "n_harmful": int(h.sum()), "n_harmless": int((~h).sum()),
                     "value": a - ab, "reads": f"95% bootstrap interval {lo:.3f} to {hi:.3f}"})
        return lo

    h = ok.harmful.to_numpy()
    if h.any() and not h.all():
        area(ok, MAIN, "all", "area")
        area(ok, BASELINE, "all", "area")
        lo = difference(ok, MAIN, BASELINE, "all", "difference")
        verdict = "counts" if cybrid.auc(ok[MAIN], ok.harmful) >= AUC_MIN and lo > 0 else "does not count"
        for s in ("memb0", "pmpnn", "smpnn", "base_blosum"):
            area(ok, s, "all")
        difference(ok, MAIN, "smpnn", "all", "difference (reported, claimed as nothing)")
        for name, d in (("contact", ok[ok.contact]), ("not contact", ok[~ok.contact])):
            for s in (MAIN, BASELINE):
                area(d, s, name)
        for c, d in ok.groupby("complex"):
            for s in (MAIN, BASELINE):
                area(d, s, f"complex {c}")
    else:
        verdict = "not assessed"
    rows.append({"what": "rule", "group": "all", "score": MAIN, "n_harmful": int(h.sum()), "n_harmless": int((~h).sum()), "value": np.nan,
                 "reads": f"{verdict}: an area of at least {AUC_MIN} and above the site's variability with the interval above 0"})
    res = pd.DataFrame(rows)
    res.to_csv(out / "h34_mpnn_summary.tsv", sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.width", 220, "display.max_colwidth", 110, "display.max_rows", 200):
        print(res.to_string(index=False))
    print("wrote", out / "h34_mpnn_summary.tsv")
    return res
