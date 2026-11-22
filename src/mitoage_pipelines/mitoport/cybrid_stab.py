"""H28, a learned binding energy: does StaB-ddG predict that the donor's mitochondrial subunits bind the host's nuclear
subunits worse where the measured complexes fail?

StaB-ddG (Deng et al., ICML 2025; github.com/LDeng0205/StaB-ddG at COMMIT) is ProteinMPNN fine-tuned on folding
stabilities and then on measured binding energy changes. For a set of substitutions it predicts

    ddG = ddG_fold(complex) - ddG_fold(partner one alone) - ddG_fold(partner two alone)

in kcal/mol, positive for weaker binding, each term from the model's log-likelihood of the whole sequence on the given
backbone, averaged over MC_SAMPLES draws of backbone noise and decoding order. Here partner one is the complex's
mtDNA-encoded chains and partner two the nuclear chains of the system, on a solved structure of the host's complex.

Fixed here, before any score was computed:
- Structures: H18's (`cybrid_energy.STRUCTURES`). The system is H18's too: the complex's mtDNA-encoded chains (every
  copy, found by sequence as in H18) and each protein chain with a heavy atom within 4.5 A of one of them; a chain
  that is an mtDNA-encoded protein of another complex is dropped. Only the backbone (N, CA, C, O) as deposited is
  read; a residue that lacks one of the four is known to the model by its name and masked, as the model's own reader
  does. A non-standard residue is read as its standard parent (Biopython's table), any other as unknown.
- A break in a chain's numbering is kept as that many masked places, at most GAP_CAP (the model's positional
  encoding does not tell longer distances apart).
- Substitutions only, H18's rule (`cybrid_energy.substitutions`); deletions, insertions and unresolved positions
  are counted (`unmodelled`), never modelled.
- Sets, each scored on its own: `all`, every substitution of the complex's mtDNA-encoded chains; `contact` and
  `control`, those at the contact residues and at their burial-matched controls. Contacts and controls are H17's
  (`mitoport burial`, sites.tsv of the human structures), carried to the host's protein by a global alignment with
  the human one. `none`, no substitution, must give exactly 0.
- Weights: `stabddg` (the score), and as controls `proteinmpnn` (plain ProteinMPNN) and `stability_finetuned`. Seeds
  SEEDS; the model's defaults otherwise (backbone noise 0.1, 20 samples, antithetic variates).

`cybrid-stab` runs one (host, donor, complex) and writes stab/<host>.<donor>.<complex>.tsv, one row per weights, seed
and set; `cybrid-stab-collect` joins them into stab_scores.tsv (the mean over the seeds and its spread) and
stab_by_seed.tsv. Needs `uv sync --group stab` and the weights (`mitoage download stabddg`).
"""
from pathlib import Path

import numpy as np
import pandas as pd

from . import cybrid
from .cybrid_energy import STRUCTURES, THREE, insertions, match_protein, substitutions, triples

REPO = "https://github.com/LDeng0205/StaB-ddG"
COMMIT = "b19a5d72637a6b581b16f38af39383658e5ae1a5"
WEIGHTS = {       # name -> (file in the repository's model_ckpts, SHA-256)
    "stabddg": ("stabddg.pt", "51a0feb6f9e4f296dbd31db2470c2b8d917d10f64f1302ae63650f43095e348f"),
    "proteinmpnn": ("proteinmpnn.pt", "7af52d090172c230c7f0e9d21e02203f6b3a38b16db58d3c7a3960e0a9a6e31a"),
    "stability_finetuned": ("stability_finetuned.pt", "dc645acee165a7d6072699e9ff6278f695fa5f88a0f85cfb42c28cd2380a32eb"),
}
PREFIX = {"stabddg": "stab", "proteinmpnn": "mpnn", "stability_finetuned": "stabfold"}
SEEDS = (0, 1, 2)
MC_SAMPLES = 20
NOISE = 0.1
CONTACT_A = 4.5
GAP_CAP = 33
SETS = ("all", "contact", "control", "none")
ALPHABET = "ACDEFGHIKLMNPQRSTVWYX"
BACKBONE = ("N", "CA", "C", "O")


def read_chains(cif: Path) -> list[dict]:
    """The protein chains of the entry's first model, in the file's order: {id, seq (one letter per residue, X when
    unknown), num (the residues' numbers), bb [n, 4, 3] (N, CA, C, O; NaN where an atom is missing), heavy [m, 3]}.
    A residue counts as protein when it is a standard amino acid or has a standard parent, and has a CA."""
    from Bio.Data.PDBData import protein_letters_3to1_extended as LETTER
    from Bio.PDB.MMCIFParser import MMCIFParser
    model = next(iter(MMCIFParser(QUIET=True).get_structure(cif.stem, str(cif))))
    out = []
    for chain in model:
        seq, num, bb, hv = [], [], [], []
        for res in chain:
            name = res.get_resname().strip()
            if name not in LETTER or "CA" not in res:
                continue
            seq.append(LETTER[name] if LETTER[name] in THREE else "X"); num.append(res.id[1])
            bb.append([res[a].coord if a in res else np.full(3, np.nan) for a in BACKBONE])
            hv += [a.coord for a in res if a.element != "H"]
        if seq:
            out.append({"id": chain.id, "seq": "".join(seq), "num": num, "bb": np.array(bb, dtype=float), "heavy": np.array(hv, dtype=float)})
    return out


def system(chains: list[dict], mt_proteins: dict[str, str], other_mt: dict[str, str]) -> tuple[list[dict], dict[int, tuple[str, dict[int, int]]]]:
    """(the chains kept, {index among them: (gene, {residue index in the chain: 0-based protein position})} for the
    mtDNA-encoded ones): the complex's mtDNA-encoded chains and every chain within CONTACT_A of one, without the
    chains that are an mtDNA-encoded protein of another complex. Refuses when a gene of the complex has no chain."""
    from scipy.spatial import cKDTree
    al = cybrid.aligner()
    hits = [match_protein(c["seq"], {**mt_proteins, **other_mt}, al) for c in chains]
    mine = [i for i, h in enumerate(hits) if h and h[0] in mt_proteins]
    found = {hits[i][0] for i in mine}
    if not set(mt_proteins) <= found:
        raise SystemExit(f"no chain of {sorted(set(mt_proteins) - found)}")
    tree = cKDTree(np.concatenate([chains[i]["heavy"] for i in mine]))
    keep = []
    for i, c in enumerate(chains):
        if i in mine:
            keep.append(i)
        elif hits[i] is None and np.isfinite(tree.query(c["heavy"], distance_upper_bound=CONTACT_A)[0]).any():
            keep.append(i)
    return [chains[i] for i in keep], {k: hits[i] for k, i in enumerate(keep) if i in mine}


def gapped(chain: dict) -> tuple[str, np.ndarray, list[int]]:
    """(sequence with '-' at the breaks, backbone with NaN there, {residue index -> index in the gapped sequence} as
    a list). A break is a step of more than one in the residues' numbers; it becomes that many places less one, at
    most GAP_CAP. A step of one or less (an insertion code, a numbering that runs back) adds nothing."""
    seq, bb, where = [], [], []
    for i, (aa, n) in enumerate(zip(chain["seq"], chain["num"])):
        if i:
            gap = min(max(n - chain["num"][i - 1] - 1, 0), GAP_CAP)
            seq += ["-"] * gap; bb += [np.full((4, 3), np.nan)] * gap
        where.append(len(seq)); seq.append(aa); bb.append(chain["bb"][i])
    return "".join(seq), np.array(bb), where


def structure(name: str, chains: list[dict], which: list[int]) -> tuple[dict, dict[tuple[int, int], int]]:
    """The model's own description of the chains `which` (the dictionary its reader makes from a PDB file), and
    {(chain index, residue index in the chain): index in its concatenated sequence}."""
    d = {"name": name, "num_of_chains": len(which), "masked_list": [], "visible_list": []}
    at, total, offset = {}, [], 0
    for c in which:
        seq, bb, where = gapped(chains[c])
        cid = f"c{c}"
        d[f"seq_chain_{cid}"] = seq
        d[f"coords_chain_{cid}"] = {f"{a}_chain_{cid}": bb[:, k, :].tolist() for k, a in enumerate(BACKBONE)}
        d["masked_list"].append(cid)
        at.update({(c, i): offset + w for i, w in enumerate(where)})
        total.append(seq); offset += len(seq)
    d["seq"] = "".join(total)
    return d, at


def encode(seq: str, subs: dict[int, str]) -> np.ndarray:
    """The sequence as the model's indices, with the residues at `subs` replaced."""
    s = list(seq)
    for i, aa in subs.items():
        if s[i] not in THREE or aa not in THREE:
            raise ValueError(f"substitution at a place that holds no amino acid: {s[i]} to {aa} at {i}")
        s[i] = aa
    return np.array([ALPHABET.index(a) if a in ALPHABET else ALPHABET.index("X") for a in s], dtype=np.int64)


def site_groups(sites: pd.DataFrame, human: dict[str, str], host: dict[str, str]) -> dict[str, dict[int, str]]:
    """{gene: {0-based position in the host's protein: "contact" or "control"}} from the human structures' sites,
    carried by a global alignment of the human protein with the host's. Refuses a site whose residue is not the
    human protein's."""
    al = cybrid.aligner(); out = {}
    for g in host:
        s = sites[(sites.gene == g) & sites.group.isin(["contact", "control"])]
        for r in s.itertuples():
            if human[g][r.human_pos - 1] != r.human_aa:
                raise SystemExit(f"sites.tsv: {g} {r.human_pos} is {r.human_aa}, the human protein has {human[g][r.human_pos - 1]}")
        a = al.align(human[g], host[g])[0]
        to_host = {h0 + k: o0 + k for (h0, h1), (o0, o1) in zip(*a.aligned) for k in range(h1 - h0)}
        out[g] = {to_host[r.human_pos - 1]: r.group for r in s.itertuples() if r.human_pos - 1 in to_host}
    return out


def plan(chains: list[dict], mt_chain: dict[int, tuple[str, dict[int, int]]], host: dict[str, str],
         donor_on: dict[str, list[str]], groups: dict[str, dict[int, str]]) -> tuple[dict[str, dict[tuple[int, int], str]], int]:
    """({set: {(chain index, residue index): the donor's residue}}, positions not modelled) for one (host, donor,
    complex); every copy of a chain carries its substitutions."""
    sets = {s: {} for s in SETS}; unmodelled = 0
    for c, (g, to_protein) in sorted(mt_chain.items()):
        subs, u = substitutions(chains[c]["seq"], to_protein, host[g], donor_on[g])
        unmodelled += u
        for i, aa in subs.items():
            sets["all"][(c, i)] = aa
            grp = groups[g].get(to_protein[i])
            if grp:
                sets[grp][(c, i)] = aa
    return sets, unmodelled


def load_model(weights: Path, device: str, noise: float = NOISE):
    import torch
    from stabddg.model import StaBddG
    from stabddg.mpnn_utils import ProteinMPNN
    pm = ProteinMPNN(node_features=128, edge_features=128, hidden_dim=128, num_encoder_layers=3, num_decoder_layers=3,
                     k_neighbors=48, dropout=0.0, augment_eps=0.0)
    ck = torch.load(weights, map_location="cpu")
    pm.load_state_dict(ck["model_state_dict"] if "model_state_dict" in ck else ck)
    return StaBddG(pmpnn=pm, use_antithetic_variates=True, noise_level=noise, device=device).to(device).eval()


def predict(model, structs: tuple[dict, dict, dict], seqs: tuple[np.ndarray, np.ndarray, np.ndarray], samples: int = MC_SAMPLES) -> float:
    """The predicted change of binding energy for one set of substitutions, as the repository's own script computes
    it: the mean over `samples` draws, sign turned so that positive is weaker binding."""
    import torch
    vals = []
    with torch.no_grad():
        t = [torch.from_numpy(s)[None, :].to(model.device) for s in seqs]
        for _ in range(samples):
            vals.append(float(model(*structs, *t).squeeze().cpu()))
    return -float(np.mean(vals))


def sha256(path: Path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(host: str, donor: str, complex_: str, gb: Path, pdb_dir: Path, sites_dir, weights_dir: Path, out: Path,
        device: str = "cpu", seeds: tuple[int, ...] = SEEDS, samples: int = MC_SAMPLES, model_for=None) -> pd.DataFrame:
    """`sites_dir(complex)` is the folder of the complex's sites.tsv. `model_for(weights name)` replaces the model
    (for tests); a model is a callable of (complex, partner one, partner two, three index tensors) with `.device`."""
    from .complexes import COMPLEXES
    pid = STRUCTURES[host][complex_]
    cif = pdb_dir / f"{pid}.cif"
    if not cif.is_file():
        raise SystemExit(f"{cif} is missing: `mitoage download pdb {pid}`")
    records = cybrid.read_proteins(gb)
    hp, dp, hum = (cybrid.species_proteins(records, s) for s in (host, donor, cybrid.HUMAN))
    genes = tuple(COMPLEXES[complex_].mt)
    mt = {g: hp[g][0] for g in genes}
    chains, mt_chain = system(read_chains(cif), mt, {g: p for g, (p, _, _) in hp.items() if g not in genes})
    al = cybrid.aligner()
    donor_on = {g: cybrid.on_human(mt[g], dp[g][0], dp[g][1], al) for g in genes}
    sites = pd.read_csv(sites_dir(complex_) / "sites.tsv", sep="\t", dtype={"chain": str})
    sets, unmodelled = plan(chains, mt_chain, mt, donor_on, site_groups(sites, {g: hum[g][0] for g in genes}, mt))
    unmodelled += sum(insertions(mt[g], dp[g][0], al) for g in genes if not dp[g][1])
    one, two = sorted(mt_chain), [i for i in range(len(chains)) if i not in mt_chain]
    if not two:
        raise SystemExit(f"{pid}: no nuclear chain touches the mtDNA-encoded chains")
    parts = [structure(pid, chains, sorted(one + two)), structure(f"{pid}_mt", chains, one), structure(f"{pid}_nuc", chains, two)]
    base = {"host": host, "donor": donor, "complex": complex_, "pdb": pid, "chains": len(chains), "mt_chains": len(one),
            "places": len(parts[0][0]["seq"]), "unmodelled": unmodelled, "device": device, "samples": samples}
    rows = []
    for name, (file, digest) in WEIGHTS.items():
        if model_for is None:
            w = weights_dir / file
            if not w.is_file() or sha256(w) != digest:
                raise SystemExit(f"{w} is missing or is not the pinned file: `mitoage download stabddg`")
            model = load_model(w, device)
        else:
            model = model_for(name)
        for seed in seeds:
            if model_for is None:
                import torch
                torch.manual_seed(seed)
            for s in SETS:
                seqs = tuple(encode(d["seq"], {at[k]: aa for k, aa in sets[s].items() if k in at}) for d, at in parts)
                v = predict(model, tuple(d for d, _ in parts), seqs, samples)
                rows.append({**base, "weights": name, "seed": seed, "set": s, "substituted": len(sets[s]), "ddG": v})
                print(f"{host} {donor} {complex_} {name} seed {seed} {s}: {len(sets[s])} substitutions, ddG {v:.4f}", flush=True)
    res = pd.DataFrame(rows)
    d = out / "stab"; d.mkdir(parents=True, exist_ok=True)
    f = d / f"{host}.{donor}.{complex_}.tsv"
    res.to_csv(f, sep="\t", index=False, float_format="%.6f")
    print("wrote", f)
    return res


def collect(out: Path, measurements: Path = cybrid.MEASUREMENTS) -> pd.DataFrame:
    """stab_scores.tsv (one row per host, donor and complex: each weights' sets as the mean over the seeds, with the
    standard deviation over them) and stab_by_seed.tsv (the same per seed) from the per-task tables. A task without a
    table is listed as missing, never as zero; a table whose empty set is not exactly 0 is refused."""
    want = triples(cybrid.read_measurements(measurements))
    tabs, missing = [], []
    for h, d, c in want:
        f = out / "stab" / f"{h}.{d}.{c}.tsv"
        if f.is_file():
            tabs.append(pd.read_csv(f, sep="\t"))
        else:
            missing.append(f"{h} / {d} / {c}")
    if missing:
        print(f"NOT ASSESSED: {len(missing)} of {len(want)} tasks have no table: " + "; ".join(missing))
    if not tabs:
        raise SystemExit("no table to collect")
    long = pd.concat(tabs, ignore_index=True)
    bad = long[(long.set == "none") & (long.ddG != 0)]
    if len(bad):
        raise SystemExit(f"the empty set is not 0 in {len(bad)} rows: the two sides of a prediction do not share their draws")
    long = long[long.set != "none"]
    long["score"] = long.weights.map(PREFIX) + "_" + long.set
    key = ["host", "donor", "complex"]
    by_seed = long.pivot_table(index=key + ["seed"], columns="score", values="ddG", sort=False).reset_index()
    by_seed.columns.name = None
    g = long.groupby(key + ["score"], sort=False).ddG
    mean = g.mean().unstack("score"); sd = g.std(ddof=1).unstack("score").add_suffix("_sd")
    n = long[long.weights == "stabddg"].pivot_table(index=key, columns="set", values="substituted", aggfunc="first", sort=False).add_prefix("n_")
    info = long.groupby(key, sort=False)[["pdb", "chains", "mt_chains", "places", "unmodelled", "device", "samples"]].first()
    seeds = long.groupby(key, sort=False).seed.nunique().rename("seeds")
    res = pd.concat([info, seeds, n, mean, sd], axis=1).reset_index()
    res.columns.name = None
    res.to_csv(out / "stab_scores.tsv", sep="\t", index=False, float_format="%.5f")
    by_seed.to_csv(out / "stab_by_seed.tsv", sep="\t", index=False, float_format="%.5f")
    with pd.option_context("display.width", 250, "display.max_rows", 200, "display.max_columns", 40):
        print(res.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("wrote", out / "stab_scores.tsv", out / "stab_by_seed.tsv")
    return res
