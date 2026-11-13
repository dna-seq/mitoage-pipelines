"""HyPhy RELAX per gene on the inputs `relax_prep` wrote: is selection relaxed or intensified on the
Foreground (long-lived) branches against the Reference (short-lived) ones."""
from pathlib import Path

from .. import external
from ..genes import MT_GENES


def main(relax_dir: Path, out_dir: Path, genes: list[str] | None = None, cpus: int = 1, minimal: bool = False) -> list[Path]:
    """Runs `hyphy relax --alignment <GENE>.fna --tree tree.nwk --code Vertebrate-mtDNA --test Foreground
    --reference Reference --output <out_dir>/<GENE>.relax.json` for each gene; returns the result files.
    `minimal` adds `--models Minimal` (the test alone, without the descriptive models)."""
    hyphy = external.find("hyphy")
    tree = relax_dir / "tree.nwk"
    if not tree.is_file():
        raise FileNotFoundError(f"{tree} is missing: run `mitoage longevity relax-prep` first")
    out_dir.mkdir(parents=True, exist_ok=True)
    done = []
    for g in genes or MT_GENES:
        out = out_dir / f"{g}.relax.json"
        cmd = [hyphy, "relax", f"CPU={cpus}", "--alignment", relax_dir / f"{g}.fna", "--tree", tree,
               "--code", "Vertebrate-mtDNA", "--test", "Foreground", "--reference", "Reference"]
        if minimal: cmd += ["--models", "Minimal"]
        external.run(cmd + ["--output", out])
        done.append(out)
    return done
