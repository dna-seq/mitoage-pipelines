"""Branch lengths of every protein on the one fixed topology, with IQ-TREE.

`--keep-ident` is essential: without it IQ-TREE drops identical sequences and puts them back with
zero-length branches, on a topology of its own for that group, and the branches of the proteins
no longer correspond.
"""
from pathlib import Path

from .. import external


def model(gene: str) -> str:
    return "mtMAM+G4" if gene.startswith("MT-") else "LG+G4"


def main(aln_dir: Path, topology: Path, genes_file: Path, bl_dir: Path, threads: int = 2,
         seed: int | None = None) -> list[str]:
    """Runs `iqtree -s aln/<G>.faa -te topology.nwk -m <model> -T <threads> --keep-ident --prefix bl/<G>`
    for each gene of genes.txt; returns the genes IQ-TREE failed on (reported, not hidden).
    IQ-TREE draws a random seed unless `seed` is given: repeated runs differ in late digits of some lengths."""
    iqtree = external.find("iqtree")
    bl_dir.mkdir(parents=True, exist_ok=True)
    failed = []
    for g in genes_file.read_text().split():
        try:
            cmd = [iqtree, "-s", aln_dir / f"{g}.faa", "-te", topology, "-m", model(g), "-T", str(threads),
                   "--keep-ident", "--prefix", bl_dir / g, "-redo", "-quiet"]
            if seed is not None: cmd += ["--seed", str(seed)]
            external.run(cmd)
        except Exception as e:      # one protein failing does not stop the others; the count is printed
            failed.append(g); print(f"FAILED {g}: {type(e).__name__}", flush=True)
    done = sorted(p.name for p in bl_dir.glob("*.treefile"))
    print(f"branch lengths: {len(done)} tree files; failed: {failed or 'none'}")
    return failed
