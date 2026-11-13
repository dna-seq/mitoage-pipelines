"""Rate covariation, preparation.

Nuclear proteins (the selections, <GENE>.faa): align with MAFFT L-INS-i, drop a sequence under 50 % identity
to human over aligned columns or covering under 70 % of human (wrong paralog or broken gene model), keep
columns with >= 50 % occupancy. Mitochondrial proteins: the build's alignments, same column rule.
Then the complete matrix: species present in every kept gene; a gene missing more than 12 % of the
species that most genes have is dropped instead. One fixed topology (the mitochondrial IQ-TREE one, pruned).

Writes MAFFT's alignments as it made them (mafft/<GENE>.afa), the filtered ones (aln/<GENE>.faa, the
mitochondrial ones as MT-<GENE>.faa), topology.nwk, genes.txt and species.txt.
"""
import io
from pathlib import Path

from Bio import Phylo

from .. import external
from ..genes import MT_GENES
from ..seqio import read_fasta


def mafft_linsi(mafft: str, faa: Path, out: Path, threads: int) -> dict[str, str]:
    """MAFFT L-INS-i of `faa`, kept in `out` as MAFFT wrote it. With more than one thread MAFFT's
    iterative refinement is not repeatable run to run (seen here: one alignment in a few differs in a row)."""
    with open(out, "w") as f:
        external.run([mafft, "--localpair", "--maxiterate", "1000", "--quiet", "--thread", str(threads), faa], stdout=f)
    return read_fasta(out)


def main(data: Path, sel: Path, tree_file: Path, outgroup_file: Path, mafft_dir: Path, aln_dir: Path,
         topology: Path, genes_file: Path, species_file: Path, threads: int = 4) -> None:
    mafft = external.find("mafft", env="MAFFT")
    for d in (mafft_dir, aln_dir, topology.parent, genes_file.parent, species_file.parent):
        d.mkdir(parents=True, exist_ok=True)
    selections = sorted(sel.glob("*.faa"))
    if not selections:
        raise FileNotFoundError(f"no selections under {sel}: run `mitoage mitoport select` first")
    alns, dropped = {}, {}
    for p in selections:
        g = p.stem
        a = mafft_linsi(mafft, p, mafft_dir / f"{g}.afa", threads)
        h = a["homo-sapiens"]; hcols = [i for i, c in enumerate(h) if c != "-"]; keep = {}
        for s, r in a.items():
            cov = [i for i in hcols if r[i] != "-"]
            ident = sum(r[i] == h[i] for i in cov) / max(1, len(cov))
            if len(cov) < 0.7 * len(hcols) or ident < 0.5: dropped.setdefault(g, []).append(f"{s}:{len(cov)/len(hcols):.2f}/{ident:.2f}")
            else: keep[s] = r
        alns[g] = keep
    for g in MT_GENES:
        alns["MT-" + g] = read_fasta(data / "alignments" / f"{g}.aa.fasta.gz")
    nuc = [g for g in alns if not g.startswith("MT-")]
    anynuc = set.union(*[set(alns[g]) for g in nuc])
    count = {s: sum(s in alns[g] for g in nuc) for s in anynuc}
    core = {s for s, c in count.items() if c >= 0.9 * len(nuc)}            # species most nuclear genes have
    genes = [g for g in alns if len(core - set(alns[g])) <= 0.12 * len(core)]
    species = sorted(set.intersection(*[set(alns[g]) for g in genes]) & core)
    print(f"nuclear genes {len(nuc)}; species with >=90% of them {len(core)}; genes kept {len(genes)} (dropped: {sorted(set(alns) - set(genes))}); species in every kept gene {len(species)}")
    for g, d in dropped.items(): print(f"  {g}: {len(d)} sequences dropped (coverage/identity):", " ".join(d[:6]), "..." if len(d) > 6 else "")
    for g in genes:
        rows = {s: alns[g][s] for s in species}; L = len(next(iter(rows.values())))
        cols = [i for i in range(L) if sum(r[i] != "-" for r in rows.values()) >= 0.5 * len(rows)]
        (aln_dir / f"{g}.faa").write_text("".join(f">{s.replace('-', '_')}\n{''.join(r[i] for i in cols)}\n" for s, r in rows.items()))
    tree = Phylo.read(tree_file, "newick")
    outg = outgroup_file.read_text().split()
    tree.root_with_outgroup(*[t for t in tree.get_terminals() if t.name in outg][:1])
    for t in tree.get_terminals():
        if t.name not in species: tree.prune(t)
    for cl in tree.find_clades(): cl.confidence = None; cl.branch_length = None
    for t in tree.get_terminals(): t.name = t.name.replace("-", "_")
    buf = io.StringIO(); Phylo.write(tree, buf, "newick", plain=True); topology.write_text(buf.getvalue())
    genes_file.write_text("\n".join(genes) + "\n"); species_file.write_text("\n".join(species) + "\n")
