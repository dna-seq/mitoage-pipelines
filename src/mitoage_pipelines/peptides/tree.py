"""RM16's third question, first part: how a frame was kept and lost on the build's dated species tree.

Per peptide, reading and group (each class with at least SPECIES_MIN assessed species in the tree, and all species):
the minimum number of changes between kept and lost (Sankoff with unit costs, which handles the tree's polytomies
exactly; a species not assessed, or outside the group, is a missing tip and costs nothing), and one most
parsimonious reconstruction: the root takes the cheaper state (kept on a tie, the published state), and each child
the state that keeps the cost minimal, its parent's on a tie (changes as close to the tips as parsimony allows).
Losses and gains are counted on that reconstruction; their sum is the minimum.
Clustering: the minimum number of changes against PERMUTATIONS random reassignments of the same states to the same
tips (Maddison & Slatkin 1991): fewer changes than random means kept and lost are clustered on the tree.
p_clustered = (1 + permutations with at most the observed changes) / (1 + permutations).
Writes tree.tsv (from states.tsv).
"""
import hashlib
import io
import sys
from pathlib import Path

import numpy as np
import polars as pl
from Bio import Phylo

from . import frames

PERMUTATIONS = 1000
SEED = 20261115


class Tree:
    """The tree as arrays: nodes in postorder, each node's children, the tip names."""

    def __init__(self, newick: Path | str):
        sys.setrecursionlimit(100000)
        src = newick if isinstance(newick, Path) else io.StringIO(newick)
        nodes = list(Phylo.read(src, "newick").find_clades(order="postorder"))
        index = {id(c): i for i, c in enumerate(nodes)}
        self.children = [[index[id(k)] for k in c.clades] for c in nodes]
        self.tip = {c.name: i for i, c in enumerate(nodes) if not c.clades}
        self.root = len(nodes) - 1


def sankoff(tree: Tree, states: dict[str, np.ndarray]) -> np.ndarray:
    """Minimum number of changes for each column of `states` ({tip: (columns,) of 1, 0, -1})."""
    return np.min(_costs(tree, states)[tree.root], axis=0)


def _costs(tree: Tree, states: dict[str, np.ndarray]) -> list[np.ndarray]:
    width = len(next(iter(states.values())))
    inf = 1e9
    cost: list = [None] * len(tree.children)
    for name, i in tree.tip.items():
        s = states.get(name)
        c = np.zeros((2, width))
        if s is not None:
            c[0, s == 1] = inf
            c[1, s == 0] = inf
        cost[i] = c
    for i, ch in enumerate(tree.children):
        if not ch:
            continue
        c = np.zeros((2, width))
        for k in ch:
            ck = cost[k]
            c[0] += np.minimum(ck[0], ck[1] + 1)
            c[1] += np.minimum(ck[1], ck[0] + 1)
        cost[i] = c
    return cost


def reconstruct(tree: Tree, states: dict[str, int]) -> tuple[int, int, int, int]:
    """One most parsimonious reconstruction of a single column. -> (minimum changes, losses, gains, root state)."""
    cost = _costs(tree, {k: np.array([v]) for k, v in states.items()})
    c = [x[:, 0] for x in cost]
    root = 1 if c[tree.root][1] <= c[tree.root][0] else 0
    losses = gains = 0
    stack = [(tree.root, root)]
    while stack:
        i, s = stack.pop()
        for k in tree.children[i]:
            stay, move = c[k][s], c[k][1 - s] + 1
            t = s if stay <= move else 1 - s
            if t != s:
                losses += s == 1
                gains += s == 0
            stack.append((k, t))
    return int(min(c[tree.root])), losses, gains, root


def seed_for(seed: int, *parts: str) -> np.random.Generator:
    h = hashlib.sha256(":".join([str(seed), *parts]).encode()).digest()
    return np.random.default_rng(int.from_bytes(h[:8], "little"))


def test_group(tree: Tree, tips: dict[str, int], permutations: int, rng: np.random.Generator) -> dict:
    """Changes, losses, gains and clustering for one group's tips ({slug: 1 kept / 0 lost})."""
    names = sorted(tips)
    vals = np.array([tips[s] for s in names])
    mat = np.empty((len(names), permutations + 1), dtype=np.int8)
    mat[:, 0] = vals
    for b in range(permutations):
        mat[:, b + 1] = rng.permutation(vals)
    changes = sankoff(tree, {s: mat[i] for i, s in enumerate(names)})
    obs, losses, gains, root = reconstruct(tree, tips)
    assert obs == changes[0]
    perm = changes[1:]
    return {"tips": len(names), "kept": int(vals.sum()), "lost": int((vals == 0).sum()), "min_changes": obs,
            "losses": losses, "gains": gains, "root_state": "kept" if root == 1 else "lost",
            "perm_mean_changes": float(perm.mean()), "p_clustered": float((1 + (perm <= obs).sum()) / (1 + permutations))}


def main(states_file: Path, tree_file: Path, out: Path, permutations: int = PERMUTATIONS, seed: int = SEED) -> pl.DataFrame:
    states = pl.read_csv(states_file, separator="\t", infer_schema_length=None)
    tree = Tree(tree_file)
    rows = []
    for (pep, reading), g in states.group_by(["peptide", "reading"], maintain_order=True):
        g = g.filter(pl.col("state") != "not assessed").filter(pl.col("slug").is_in(list(tree.tip)))
        groups = {"all": g}
        for (c,), gc in g.group_by(["class"]):
            groups[c] = gc
        for name in sorted(groups):
            gg = groups[name]
            base = {"peptide": pep, "reading": reading, "group": name}
            if gg.height < frames.SPECIES_MIN:
                rows.append({**base, "tips": gg.height, "status": f"not assessed: under {frames.SPECIES_MIN} assessed tips"})
                continue
            tips = dict(zip(gg["slug"], (gg["state"] == "kept").cast(int)))
            if len(set(tips.values())) < 2:
                r = {"tips": len(tips), "kept": sum(tips.values()), "lost": len(tips) - sum(tips.values()),
                     "min_changes": 0, "losses": 0, "gains": 0,
                     "root_state": "kept" if next(iter(tips.values())) == 1 else "lost"}
                rows.append({**base, **r, "status": "no variation: clustering not assessed"})
                continue
            r = test_group(tree, tips, permutations, seed_for(seed, pep, reading, name))
            rows.append({**base, **r, "status": "tested"})
    res = pl.DataFrame(rows, infer_schema_length=None).with_columns(pl.lit(seed).alias("seed")).sort(
        ["peptide", "reading", "group"])
    out.parent.mkdir(parents=True, exist_ok=True)
    res.write_csv(out, separator="\t", float_precision=6)
    return res
