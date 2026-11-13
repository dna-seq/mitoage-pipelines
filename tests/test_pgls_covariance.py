"""The Brownian covariance of the phylogenetic regression: shared path from the root, on small trees."""
import io

import numpy as np
from Bio import Phylo

from mitoage_pipelines.longevity.pgls import V, shared_path_covariance


def covariance(newick: str, order: list[str]) -> np.ndarray:
    tree = Phylo.read(io.StringIO(newick), "newick")
    tips = {t.name: i for i, t in enumerate(tree.root.get_terminals())}
    C = shared_path_covariance(tree.root, tips)
    idx = [tips[n] for n in order]
    return C[np.ix_(idx, idx)]


def from_paths(paths: dict[str, list[tuple[str, float]]], order: list[str]) -> np.ndarray:
    """paths: tip -> the branches from the root to it, as (branch label, length). Shared = common branches."""
    return np.array([[sum(l for b, l in paths[x] if (b, l) in paths[y]) for y in order] for x in order])


def test_balanced_four_tips():
    order = ["A", "B", "C", "D"]
    paths = {"A": [("ab", 3.0), ("a", 1.0)], "B": [("ab", 3.0), ("b", 2.0)],
             "C": [("cd", 6.0), ("c", 4.0)], "D": [("cd", 6.0), ("d", 5.0)]}
    C = covariance("((A:1,B:2):3,(C:4,D:5):6);", order)
    assert np.array_equal(C, from_paths(paths, order))
    assert np.array_equal(C, np.array([[4, 3, 0, 0], [3, 5, 0, 0], [0, 0, 10, 6], [0, 0, 6, 11.0]]))


def test_ladder_and_tip_order():
    paths = {"A": [("abc", 5.0), ("ab", 3.0), ("a", 1.0)], "B": [("abc", 5.0), ("ab", 3.0), ("b", 2.0)],
             "C": [("abc", 5.0), ("c", 4.0)], "D": [("d", 6.0)]}
    for order in (["A", "B", "C", "D"], ["D", "B", "A", "C"]):
        C = covariance("(((A:1,B:2):3,C:4):5,D:6);", order)
        assert np.array_equal(C, from_paths(paths, order))
        assert np.array_equal(C, C.T)


def test_missing_branch_length_counts_as_zero():
    C = covariance("((A:1,B):3,C:2);", ["A", "B", "C"])
    assert np.array_equal(C, np.array([[4, 3, 0], [3, 3, 0], [0, 0, 2.0]]))


def test_lambda_scales_only_the_shared_part():
    C = covariance("((A:1,B:2):3,(C:4,D:5):6);", ["A", "B", "C", "D"])
    v = V(0.5, C)
    assert np.array_equal(np.diag(v), np.diag(C))
    off = ~np.eye(4, dtype=bool)
    assert np.array_equal(v[off], 0.5 * C[off])
    assert np.array_equal(V(1.0, C), C)
