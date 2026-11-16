"""Trait clades on the dated tree, one pair per origin, and matched control clades (H21's design, H21-H25).

The tree is the build's dated species tree, pruned to the species of a test (a tip whose state is unknown is pruned
first, so it can neither join nor split a clade); a node left with one child is merged with it, its branch lengths
added.

- A trait clade is a maximal clade whose tips all carry the trait (state 1). Its stem is the branch above its root.
  Its stem age is the age of the node that branch hangs from.
- One pair per origin (R2, L3): for each trait clade, the trait species nearest the clade's root (ties by name); its
  nearest state-0 relative outside the clade (the state-0 tips under the lowest ancestor that has any, the nearest of
  them; ties by name); and an outgroup, the tip nearest the root of the smallest clade holding both, outside it.
  Pairs are taken in order of increasing sister distance, and a pair whose path between the sisters shares a branch
  with a pair already taken is dropped, with its reason.
- Control clades (R1's draws, C2's null): H21 fixed them as state-0 clades "matched on species count (within a factor
  of 2) and on stem age (nearest), without replacement". Read here as: the candidates of a trait clade are the clades
  of the tree whose tips all have state 0 and whose size is within a factor of 2 of its; of those, the POOL nearest
  in stem age (ties by name) form its pool; a draw takes, for each trait clade in turn (largest first, ties by name),
  one member of its pool at random that shares no tip with a clade already taken in that draw. Draw k uses the
  generator seeded with (0, k). A trait clade whose pool runs out in a draw is left out of that draw and counted.
"""
import sys
from dataclasses import dataclass, field

import numpy as np
from Bio import Phylo

POOL = 10


@dataclass(eq=False)
class Node:
    name: str | None
    length: float
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None

    @property
    def is_tip(self) -> bool:
        return not self.children


def read(tree_file) -> Node:
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 100_000))
    def conv(c) -> Node:
        n = Node(c.name if c.is_terminal() else None, float(c.branch_length or 0.0))
        for ch in c.clades:
            k = conv(ch); k.parent = n; n.children.append(k)
        return n
    return conv(Phylo.read(tree_file, "newick").root)


def prune(node: Node, keep: set[str]) -> Node | None:
    """A copy of the tree with only the tips in `keep`; unary nodes merged, lengths added."""
    def rec(n: Node) -> Node | None:
        if n.is_tip:
            return Node(n.name, n.length) if n.name in keep else None
        kids = [k for k in (rec(c) for c in n.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            kids[0].length += n.length
            return kids[0]
        m = Node(None, n.length, kids)
        for k in kids: k.parent = m
        return m
    r = rec(node)
    if r is not None: r.parent = None
    return r


def walk(node: Node):
    stack = [node]
    while stack:
        n = stack.pop(); yield n; stack.extend(reversed(n.children))


def tips(node: Node) -> list[str]:
    return [n.name for n in walk(node) if n.is_tip]


def height(node: Node) -> float:
    """Distance from the node down to its nearest tip by the first child path (the tree is dated: all tips agree)."""
    h, n = 0.0, node
    while n.children:
        n = n.children[0]; h += n.length
    return h


def depth_to(node: Node, anc: Node) -> float:
    d, n = 0.0, node
    while n is not anc:
        d += n.length; n = n.parent
    return d


def trait_clades(root: Node, state: dict[str, int]) -> list[Node]:
    """Maximal clades whose tips (all with a known state in `state`) are all 1, in tree order."""
    pure = {}
    for n in reversed(list(walk(root))):
        pure[id(n)] = state[n.name] == 1 if n.is_tip else all(pure[id(c)] for c in n.children)
    out = []
    def rec(n, parent_pure):
        p = pure[id(n)]
        if p and not parent_pure: out.append(n)
        if not p:
            for c in n.children: rec(c, False)
    rec(root, False)
    return out


def stem_age(c: Node) -> float:
    return height(c) + c.length


def zero_clades(root: Node, state: dict[str, int]) -> list[Node]:
    """Every clade (tips included) whose tips are all 0."""
    zero = {}
    for n in reversed(list(walk(root))):
        zero[id(n)] = state[n.name] == 0 if n.is_tip else all(zero[id(c)] for c in n.children)
    return [n for n in walk(root) if zero[id(n)] and n.parent is not None]


def label(c: Node) -> str:
    t = sorted(tips(c))
    return t[0] if len(t) == 1 else f"{t[0]}..{t[-1]} ({len(t)})"


def nearest_tip(c: Node) -> str:
    return min(((depth_to(n, c), n.name) for n in walk(c) if n.is_tip))[1]


def nearest_outside(c: Node, ok) -> tuple[str, Node] | None:
    """The tip nearest `c`'s root outside `c` that satisfies `ok`, under the lowest ancestor that has one, and that
    ancestor."""
    inside = set(tips(c)); n = c
    while n.parent is not None:
        a = n.parent
        cands = [(depth_to(t, a), t.name) for t in walk(a) if t.is_tip and t.name not in inside and ok(t.name)]
        if cands:
            return min(cands)[1], a
        n = a
    return None


def path_edges(a: Node, b: Node, anc: Node) -> set[int]:
    out = set()
    for n in (a, b):
        while n is not anc:
            out.add(id(n)); n = n.parent
    return out


def origin_pairs(root: Node, state: dict[str, int]) -> list[dict]:
    """One pair per trait clade (see the module's docstring); every clade gets a row, taken or not, with its reason."""
    tipnode = {n.name: n for n in walk(root) if n.is_tip}
    rows, cand = [], []
    for c in trait_clades(root, state):
        t = nearest_tip(c)
        found = nearest_outside(c, lambda s: state[s] == 0)
        base = {"clade": label(c), "clade_size": len(tips(c)), "stem_age": stem_age(c), "sp_trait": t}
        if found is None:
            rows.append({**base, "taken": False, "reason": "no state-0 species outside the clade"}); continue
        r, anc = found
        og = nearest_outside(anc, lambda s: True)
        if og is None:
            rows.append({**base, "sp_other": r, "taken": False, "reason": "no outgroup outside the pair's clade"}); continue
        d = depth_to(tipnode[t], anc) + depth_to(tipnode[r], anc)
        cand.append(({**base, "sp_other": r, "outgroup": og[0], "sister_distance": d},
                     path_edges(tipnode[t], tipnode[r], anc)))
    used: set[int] = set()
    for row, edges in sorted(cand, key=lambda x: (x[0]["sister_distance"], x[0]["sp_trait"])):
        if edges & used:
            rows.append({**row, "taken": False, "reason": "path shares a branch with a pair already taken"})
        else:
            used |= edges; rows.append({**row, "taken": True, "reason": ""})
    return rows


def control_matches(root: Node, state: dict[str, int], clades: list[Node], n_draws: int,
                    pool: int = POOL) -> tuple[list[list[tuple[Node, Node]]], list[int]]:
    """`n_draws` draws of control clades, each as (trait clade, the control clade drawn for it) pairs (module
    docstring). -> (draws, trait clades left out per draw)."""
    zs = [(z, len(tips(z)), stem_age(z), label(z), frozenset(tips(z))) for z in zero_clades(root, state)]
    order = sorted(clades, key=lambda c: (-len(tips(c)), label(c)))
    pools = []
    for c in order:
        n, a = len(tips(c)), stem_age(c)
        cands = sorted((z for z in zs if n / 2 <= z[1] <= 2 * n), key=lambda z: (abs(z[2] - a), z[3]))
        pools.append(cands[:pool])
    draws, short = [], []
    for k in range(n_draws):
        rng = np.random.default_rng([0, k])
        taken: list[tuple[Node, Node]] = []; covered: set[str] = set(); missed = 0
        for c, p in zip(order, pools):
            free = [z for z in p if not (z[4] & covered)]
            if not free:
                missed += 1; continue
            z = free[int(rng.integers(len(free)))]
            taken.append((c, z[0])); covered |= z[4]
        draws.append(taken); short.append(missed)
    return draws, short


def control_draws(root: Node, state: dict[str, int], clades: list[Node], n_draws: int,
                  pool: int = POOL) -> tuple[list[list[Node]], list[int]]:
    """The control clades of each draw (`control_matches` without the trait clade each was drawn for)."""
    draws, short = control_matches(root, state, clades, n_draws, pool)
    return [[z for _, z in d] for d in draws], short


def newick(node: Node, test: set[int] | None = None, raw_names: bool = False) -> str:
    """Newick with tips named by slug ('-' written as '_' for HyPhy, unless `raw_names`); a node in `test` (by id) is
    labelled {Test}, every other {Reference}, when `test` is given."""
    def lab(n):
        return "" if test is None else ("{Test}" if id(n) in test else "{Reference}")
    def rec(n):
        if n.is_tip:
            return f"{n.name if raw_names else n.name.replace('-', '_')}{lab(n)}:{n.length:.6g}"
        return "(" + ",".join(rec(c) for c in n.children) + f"){lab(n)}:{n.length:.6g}"
    return rec(node) + ";\n"
