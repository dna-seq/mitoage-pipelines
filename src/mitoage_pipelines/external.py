"""External programs: found on PATH, never bundled. `environment.yml` pins the versions used."""
import os
import shutil
import subprocess
from pathlib import Path

# program -> the version the published results were made with
VERSIONS = {"mafft": "MAFFT 7.526 (v7.526: build 2's alignments are pinned to it)", "iqtree": "IQ-TREE 3.1.4", "hyphy": "HyPhy 2.5.104",
            "USalign": "US-align 20241108 (bioconda package usalign 20241201)", "foldmason": "FoldMason 4.dd3c235",
            "coevol": "Coevol 1.6 (github.com/bayesiancook/coevol at 2661e95fde03832a86e6b1b5a03981bad2a334ec, built from source)",
            "readcoevol": "readcoevol of Coevol 1.6 (built with coevol)", "tracecomp": "tracecomp of Coevol 1.6 (built with coevol)"}


class MissingProgram(RuntimeError):
    pass


def find(program: str, env: str | None = None) -> str:
    """Full path of `program` on PATH (or of what the environment variable `env` names)."""
    name = os.environ.get(env) if env else None
    found = shutil.which(name or program)
    if not found:
        via = f" (looked for {name!r}, from the environment variable {env})" if name else ""
        raise MissingProgram(
            f"{program} is not on PATH{via}. This step needs {VERSIONS.get(program, program)}; "
            f"`environment.yml` installs it (conda env create -f environment.yml).")
    return found


def run(cmd: list[str | Path], **kwargs) -> subprocess.CompletedProcess:
    cmd = [str(c) for c in cmd]
    print("+", " ".join(cmd), flush=True)
    return subprocess.run(cmd, check=True, **kwargs)
