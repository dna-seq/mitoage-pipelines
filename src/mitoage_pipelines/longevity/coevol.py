"""Lifespan as a continuous trait along the tree: Coevol's run and its reading.

`run`: one chain of Coevol 1.6 on the inputs of `longevity coevol-prep` (mammals.ali, mammals.tree,
mammals.lht), with dS and dN/dS (-dsom) under the mitochondrial mammal code (-mtmam), saving every cycle until
`cycles` (published: 1,100 per chain, two chains). The inputs are copied into the chain's folder, since Coevol
records them by name and its reader opens them again from there. Days of CPU per chain at 106 species.

`read`: readcoevol over each chain after a burn-in (published: 200), which writes chain<N>.cov (covariances,
correlations, partial correlations and their posterior probabilities), and tracecomp across the chains
(tracecomp.contdiff: effective sample size and relative difference between chains per parameter). Read the
.cov files only once every effective sample size is above about 100 and the relative differences are small:
the published chains reached 7 for dN/dS, and their numbers are provisional for that reason.

The published run was on the cluster on 2026-11-14 and 15 with these exact options.
"""
import shutil
from pathlib import Path

from .. import external

INPUTS = ("mammals.ali", "mammals.tree", "mammals.lht")


def run(prep: Path, out: Path, chain: int, cycles: int = 1100) -> None:
    coevol = external.find("coevol", env="COEVOL")
    out.mkdir(parents=True, exist_ok=True)
    for f in INPUTS:
        if not (prep / f).is_file():
            raise FileNotFoundError(f"{prep / f} is missing: run `mitoage longevity coevol-prep` first")
        shutil.copyfile(prep / f, out / f)
    external.run([coevol, "-d", INPUTS[0], "-t", INPUTS[1], "-c", INPUTS[2], "-dsom", "-mtmam", "-x", "1", str(cycles),
                  f"chain{chain}"], cwd=out)


def parse_cov(text: str) -> dict[str, list[list[float | None]]]:
    """The square blocks of a .cov file by heading: covariances, correlation coefficients, posterior probs,
    precisions, partial correlation coefficients, posterior probs (the second one keyed 'partial posterior probs').
    A '-' on the diagonal of a posterior block is None."""
    blocks, name, rows, seen = {}, None, [], set()
    heads = ("covariances", "correlation coefficients", "posterior probs", "precisions", "partial correlation coefficients")
    for line in text.splitlines() + [""]:
        s = line.strip()
        if s in heads:
            if name and rows: blocks[name] = rows
            name = "partial posterior probs" if s == "posterior probs" and s in seen else s
            seen.add(s); rows = []
        elif name and s and not s.startswith("entries"):
            try:
                rows.append([None if v == "-" else float(v) for v in s.split()])
            except ValueError:
                pass
        elif name and not s and rows:
            blocks[name] = rows; name, rows = None, []
    return blocks


def read(out: Path, chains: list[int], burnin: int = 200) -> None:
    readcoevol, tracecomp = external.find("readcoevol", env="READCOEVOL"), external.find("tracecomp", env="TRACECOMP")
    for c in chains:
        external.run([readcoevol, "-x", str(burnin), f"chain{c}"], cwd=out)
    if len(chains) > 1:
        external.run([tracecomp, "-x", str(burnin), *[f"chain{c}" for c in chains]], cwd=out)
        print((out / "tracecomp.contdiff").read_text() if (out / "tracecomp.contdiff").is_file() else "no tracecomp.contdiff")
    for c in chains:
        b = parse_cov((out / f"chain{c}.cov").read_text())
        print(f"chain{c}: order dS, dN/dS, Longevity, Mass")
        for k in ("correlation coefficients", "posterior probs", "partial correlation coefficients", "partial posterior probs"):
            print(f"  {k}: " + "  ".join(" ".join("-" if v is None else f"{v:+.3f}" for v in r) for r in b.get(k, [])))
    print("PROVISIONAL until every effective sample size in tracecomp.contdiff is above ~100.")
