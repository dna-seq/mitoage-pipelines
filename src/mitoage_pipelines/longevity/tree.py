"""The tree of one class from its concatenated proteins, with IQ-TREE (fast search, fixed model)."""
from pathlib import Path

from .. import external

# mitochondrial protein models: mammals, and the other vertebrates
MODELS = {"Mammalia": "mtMAM+F+G4", "Aves": "mtVer+F+G4"}


def main(tree_dir: Path, cls: str = "Mammalia", model: str | None = None, threads: str = "AUTO",
         seed: int | None = None) -> Path:
    """Runs `iqtree -s <Class>.faa -m <model> -fast -T <threads> --prefix <Class>.iq`; returns the tree file.
    IQ-TREE draws a random seed unless `seed` is given, so two runs need not give the same tree."""
    faa = tree_dir / f"{cls}.faa"
    if not faa.is_file():
        raise FileNotFoundError(f"{faa} is missing: run `mitoage longevity concat --class {cls}` first")
    if model is None:
        if cls not in MODELS:
            raise ValueError(f"no model is set for the class {cls}: name one (known: {MODELS})")
        model = MODELS[cls]
    iqtree = external.find("iqtree")
    prefix = tree_dir / f"{cls}.iq"
    cmd = [iqtree, "-s", faa, "-m", model, "-fast", "-T", threads, "--prefix", prefix, "-redo", "-quiet"]
    if seed is not None: cmd += ["--seed", str(seed)]
    external.run(cmd)
    return prefix.with_name(prefix.name + ".treefile")
