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


def run(prep: Path, out: Path, chain: int, cycles: int = 1100, fixed_ages: bool = False) -> None:
    """With `fixed_ages`, -fixtimes: the node ages are held at the tree's (coevol-prep --fixed-ages)."""
    coevol = external.find("coevol", env="COEVOL")
    out.mkdir(parents=True, exist_ok=True)
    for f in INPUTS:
        if not (prep / f).is_file():
            raise FileNotFoundError(f"{prep / f} is missing: run `mitoage longevity coevol-prep` first")
        shutil.copyfile(prep / f, out / f)
    if fixed_ages and ":" not in (out / INPUTS[1]).read_text():
        raise SystemExit(f"{prep / INPUTS[1]} has no branch lengths: write it with `coevol-prep --fixed-ages`")
    external.run([coevol, "-d", INPUTS[0], "-t", INPUTS[1], "-c", INPUTS[2], "-dsom", "-mtmam",
                  *(["-fixtimes"] if fixed_ages else []), "-x", "1", str(cycles), f"chain{chain}"], cwd=out)


def saved_points(out: Path, chain: int) -> int:
    """Points a chain has saved: the lines of chain<N>.trace after its header."""
    trace = out / f"chain{chain}.trace"
    if not trace.is_file():
        raise FileNotFoundError(f"{trace} is missing: there is no chain {chain} in {out} to resume")
    return sum(1 for _ in trace.open()) - 1


def set_until(param: str, size: int, until: int) -> str:
    """The text of a chain's .param with its stopping point set to `until`. Coevol writes `every until size` on
    one line of the header and stops once size > until; it is the one line of three integers whose last is the
    chain's saved size. Refuses unless exactly one line matches and the size there is `size`."""
    lines = param.split("\n")
    hits = [i for i, s in enumerate(lines)
            if len(f := s.split("\t")) == 3 and all(x.lstrip("-").isdigit() for x in f) and int(f[2]) == size]
    if len(hits) != 1:
        raise ValueError(f"expected one 'every until size' line with size {size} in the .param, found {len(hits)}")
    if until <= size:
        raise ValueError(f"until ({until}) must be beyond the {size} points already saved")
    every, _, _ = lines[hits[0]].split("\t")
    lines[hits[0]] = f"{every}\t{until}\t{size}"
    return "\n".join(lines)


CHECKPOINT_FILES = ("param", "trace", "monitor")     # rewritten or appended each cycle; .chain is append-only


def param_size(param: str) -> int | None:
    """The saved size Coevol recorded in a .param: the third integer of its `every until size` line, found as the
    last header line of exactly three integers (line 50 of the published chains). None when there is none."""
    for s in param.split("\n"):
        f = s.split("\t")
        if len(f) == 3 and all(x.lstrip("-").isdigit() for x in f) and int(f[0]) >= 1:
            return int(f[2])
    return None


def _line_count(path: Path) -> int:
    with path.open("rb") as f:
        return sum(buf.count(b"\n") for buf in iter(lambda: f.read(1 << 20), b""))


def consistency(out: Path, chain: int) -> tuple[bool, str, dict]:
    """Whether a chain's files agree on how many points are saved: the .trace (one line per point after its
    header), the .param's size and the .chain (the same number of lines per point, 440 in the published run).
    Coevol appends to .chain, rewrites .param, then appends to .trace on every cycle, so a chain killed between
    two of these writes disagrees with itself; resuming it would misalign points silently."""
    n = saved_points(out, chain)
    p = param_size((out / f"chain{chain}.param").read_text())
    lines = _line_count(out / f"chain{chain}.chain")
    state = {"size": n, "param_size": p, "chain_lines": lines,
             "chain_bytes": (out / f"chain{chain}.chain").stat().st_size}
    if p != n:
        return False, f"the .param records {p} points, the .trace {n}", state
    if n <= 0 or lines % n:
        return False, f"the .chain's {lines} lines are not a whole number of points for {n} points", state
    state["lines_per_point"] = lines // n
    return True, "consistent", state


def _sha256_prefix(path: Path, nbytes: int) -> str:
    import hashlib
    h, left = hashlib.sha256(), nbytes
    with path.open("rb") as f:
        while left > 0:
            buf = f.read(min(1 << 20, left))
            if not buf:
                break
            h.update(buf); left -= len(buf)
    return h.hexdigest()


def checkpoint(out: Path, chain: int, keep: int = 3) -> Path:
    """Save a consistent chain's state under checkpoints/chain<N>/<size>/: copies of .param, .trace and .monitor,
    and the .chain's length and the SHA-256 of that prefix (the .chain is append-only, so its prefix is the
    checkpoint). Refuses an inconsistent chain. Keeps the newest `keep` checkpoints."""
    import json
    ok, why, state = consistency(out, chain)
    if not ok:
        raise RuntimeError(f"chain{chain} is not consistent ({why}); no checkpoint written")
    root = out / "checkpoints" / f"chain{chain}"
    d = root / f"{state['size']:07d}"
    d.mkdir(parents=True, exist_ok=True)
    for ext in CHECKPOINT_FILES:
        shutil.copyfile(out / f"chain{chain}.{ext}", d / f"chain{chain}.{ext}")
    state["chain_sha256"] = _sha256_prefix(out / f"chain{chain}.chain", state["chain_bytes"])
    (d / "state.json").write_text(json.dumps(state, indent=1, sort_keys=True) + "\n")
    for old in sorted(p for p in root.iterdir() if p.is_dir())[:-keep]:
        shutil.rmtree(old)
    print(f"chain{chain}: checkpoint at {state['size']} points ({d})")
    return d


def restore(out: Path, chain: int) -> Path:
    """Put a chain back to its newest checkpoint: truncate the .chain to the checkpoint's length, after checking
    that its prefix is the checkpoint's, and copy back .param, .trace and .monitor. Refuses when no checkpoint
    exists or the prefix differs, and says loudly what it did otherwise."""
    import json
    root = out / "checkpoints" / f"chain{chain}"
    cps = sorted(p for p in root.iterdir() if p.is_dir()) if root.is_dir() else []
    if not cps:
        raise RuntimeError(f"chain{chain} has no checkpoint under {root}: nothing to restore, not resuming")
    d = cps[-1]
    state = json.loads((d / "state.json").read_text())
    cf = out / f"chain{chain}.chain"
    if cf.stat().st_size < state["chain_bytes"] or _sha256_prefix(cf, state["chain_bytes"]) != state["chain_sha256"]:
        raise RuntimeError(f"chain{chain}.chain does not begin with checkpoint {d.name}'s bytes: not restoring")
    with cf.open("r+b") as f:
        f.truncate(state["chain_bytes"])
    for ext in CHECKPOINT_FILES:
        shutil.copyfile(d / f"chain{chain}.{ext}", out / f"chain{chain}.{ext}")
    print(f"RESTORED chain{chain} to checkpoint {d.name} ({state['size']} points); points after it were discarded")
    return d


def resume(out: Path, chain: int, until: int, stop_after_hours: float | None = None,
           segment_hours: float = 6.0, min_segment_hours: float = 0.05) -> None:
    """Continue chain<N> in `out` until `until` points are saved, in segments of `segment_hours`, with a
    checkpoint after each. A segment ends the documented way (0 into chain<N>.run): Coevol finishes the cycle it
    is in and exits with its files in agreement. Before the first segment the chain is checked; if its files
    disagree (a run killed mid-cycle) it is restored from its newest checkpoint, and refused if it has none.
    With `stop_after_hours` no segment starts that would end after that many hours; a Slurm job gives the run a
    time limit with room for one cycle (76.6 s on haumea) beyond it; a remainder under `min_segment_hours` is
    not started."""
    import threading
    import time
    coevol = external.find("coevol", env="COEVOL")
    for f in INPUTS + (f"chain{chain}.param",):
        if not (out / f).is_file():
            raise FileNotFoundError(f"{out / f} is missing: nothing to resume")
    ok, why, _ = consistency(out, chain)
    if not ok:
        print(f"chain{chain} is inconsistent: {why}")
        restore(out, chain)
        ok, why, _ = consistency(out, chain)
        if not ok:
            raise RuntimeError(f"chain{chain} is still inconsistent after the restore ({why})")
    checkpoint(out, chain)
    start = time.monotonic()
    while True:
        size = saved_points(out, chain)
        if size > until:
            print(f"chain{chain}: {size} points saved, past until ({until}); done")
            return
        hours = segment_hours
        if stop_after_hours is not None:
            hours = min(hours, stop_after_hours - (time.monotonic() - start) / 3600)
            if hours < min_segment_hours:
                print(f"chain{chain}: {size} points saved; this run's time is spent")
                return
        param = out / f"chain{chain}.param"
        param.write_text(set_until(param.read_text(), size, until))
        print(f"chain{chain}: {size} points saved, running for up to {hours:.2f} h towards {until}", flush=True)
        timer = threading.Timer(hours * 3600, lambda: (out / f"chain{chain}.run").write_text("0\n"))
        timer.daemon = True
        timer.start()
        try:
            external.run([coevol, f"chain{chain}"], cwd=out)
        finally:
            timer.cancel()
        checkpoint(out, chain)


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
