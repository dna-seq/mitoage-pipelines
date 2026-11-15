"""The two shuffles of RM16's second question.

Synonymous shuffle (a frame inside a protein gene): the host gene's codons are permuted among the positions that
encode the same amino acid under the species' own code. The protein is unchanged, and so is the count of every
codon; only which synonymous codon sits where changes, which is all an alternative frame can see beyond the
protein. Codons that are not a sense codon of the code (unknown bases, a terminal stop) stay where they are.

Dinucleotide shuffle (a frame inside an rRNA gene, which has no synonymous codons): a uniformly random sequence with
the same dinucleotide counts and the same first and last base (Altschul & Erickson 1985), drawn as a random
Eulerian path (Kandel et al. 1996): a last exit is chosen for every base but the final one, uniformly among its
occurrences, until those exits form a tree that reaches the final base; the other exits are put in random order,
the last exit last, and the path is walked from the first base.

Every draw comes from a numpy Generator seeded from (base seed, peptide, species) by SHA-256, so a rerun gives the
same bytes and no draw depends on the order of the work.
"""
import hashlib

import numpy as np

from . import codes


def rng_for(seed: int, *parts: str) -> np.random.Generator:
    h = hashlib.sha256(":".join([str(seed), *parts]).encode()).digest()
    return np.random.default_rng(int.from_bytes(h[:8], "little"))


def residues(codons: np.ndarray, table: codes.Table) -> np.ndarray:
    """Per codon, the residue's character code under `table`; -1 for an unknown codon or a stop (kept in place)."""
    out = np.full(len(codons), -1, dtype=np.int16)
    ok = codons >= 0
    aa = np.array([ord(a) for a in table.aa], dtype=np.int16)
    out[ok] = aa[codons[ok]]
    out[ok & table.stop[np.clip(codons, 0, 63)]] = -1
    return out


def synonymous(codons: np.ndarray, aa: np.ndarray, replicates: int, rng: np.random.Generator) -> np.ndarray:
    """`replicates` shuffles of `codons` (codon indices, -1 unknown) within the classes of `aa` (-1: fixed).
    -> (replicates, len(codons)). Each row is a uniform random permutation within every class."""
    out = np.tile(codons, (replicates, 1))
    for a in np.unique(aa[aa >= 0]):
        pos = np.flatnonzero(aa == a)
        if len(pos) < 2:
            continue
        perm = np.argsort(rng.random((replicates, len(pos))), axis=1, kind="stable")
        out[:, pos] = codons[pos][perm]
    return out


def codons_to_bases(codons: np.ndarray) -> np.ndarray:
    """(..., k) codon indices -> (..., 3k) bases 0..3."""
    b = np.stack([codons // 16, (codons // 4) % 4, codons % 4], axis=-1)
    return b.reshape(*codons.shape[:-1], codons.shape[-1] * 3)


def frame_codons(bases: np.ndarray, offset: int, m: int) -> np.ndarray:
    """(..., L) bases -> (..., m) codon indices of the frame starting at `offset`."""
    b = bases[..., offset:offset + 3 * m].astype(np.int16)
    b = b.reshape(*b.shape[:-1], m, 3)
    return 16 * b[..., 0] + 4 * b[..., 1] + b[..., 2]


def kept_open(alt: np.ndarray, table: codes.Table) -> tuple[np.ndarray, np.ndarray]:
    """(..., m) codon indices of a frame's first m codons -> (kept, open): open when none is a stop, kept when also
    the first is a start codon."""
    opened = ~table.stop[alt].any(-1)
    return opened & table.start[alt[..., 0]], opened


def dinucleotide(seq: str, rng: np.random.Generator) -> str:
    """One uniformly random sequence with the dinucleotide counts, first and last base of `seq`."""
    n = len(seq)
    if n <= 2:
        return seq
    exits: dict[str, list[str]] = {}
    for a, b in zip(seq, seq[1:]):
        exits.setdefault(a, []).append(b)
    first, last = seq[0], seq[-1]
    verts = sorted(exits)
    while True:
        lastexit = {}
        for v in verts:
            if v != last:
                lastexit[v] = int(rng.integers(len(exits[v])))
        ok = True
        for v in lastexit:                         # every last exit must lead to `last` without a cycle
            seen, u = set(), v
            while u != last:
                if u in seen or u not in lastexit:
                    ok = False; break
                seen.add(u); u = exits[u][lastexit[u]]
            if not ok:
                break
        if ok:
            break
    order = {}
    for v in verts:
        lst = list(exits[v])
        tail = lst.pop(lastexit[v]) if v in lastexit else None
        idx = rng.permutation(len(lst))
        lst = [lst[i] for i in idx]
        if tail is not None:
            lst.append(tail)
        order[v] = lst
    out, u, pos = [first], first, {v: 0 for v in verts}
    for _ in range(n - 1):
        nxt = order[u][pos[u]]; pos[u] += 1
        out.append(nxt); u = nxt
    return "".join(out)


def dinucleotide_many(seq: str, replicates: int, rng: np.random.Generator) -> np.ndarray:
    """(replicates, len(seq)) bases 0..3 of independent dinucleotide shuffles."""
    return np.stack([codes.to_index(dinucleotide(seq, rng)) for _ in range(replicates)])
