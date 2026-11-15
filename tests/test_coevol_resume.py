"""Resuming a Coevol chain with checkpoints (`longevity coevol-resume`).

The .param is the published chain 1's own (haumea, rm12/out/coevol-repo, 2026-11-15, 1,101 points saved).
Coevol writes three files per cycle in this order: it appends a point to .chain, rewrites .param with the new
size, then appends a line to .trace. The tests build a chain folder around the real .param with that layout and
kill it between the writes, which is what a time limit or a node failure does."""
import os
import signal
import stat
import subprocess
import sys
import time
from pathlib import Path

import pytest

from mitoage_pipelines.longevity import coevol

DATA = Path(__file__).parent / "data"
PARAM = (DATA / "coevol_chain1.param").read_text()
SIZE = 1101          # points saved in the published chain 1
LPP = 4              # lines per point in these test chains (440 in the real one)


def make_chain(d: Path, size: int = SIZE) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    for f in coevol.INPUTS:
        (d / f).write_text("x\n")
    (d / "chain1.param").write_text(PARAM if size == SIZE else coevol.set_until(PARAM, SIZE, SIZE + 1)
                                    .replace(f"\t{SIZE + 1}\t{SIZE}", f"\t{size + 1}\t{size}"))
    (d / "chain1.trace").write_text("#header\n" + "".join(f"{i}\t0.1\n" for i in range(size)))
    (d / "chain1.chain").write_text("".join(f"p{i} l{k}\n" for i in range(size) for k in range(LPP)))
    (d / "chain1.monitor").write_text(f"monitor {size}\n")
    (d / "chain1.run").write_text("0\n")
    return d


@pytest.fixture(autouse=True)
def hard_timeout():
    """Every test here fails after 60 s instead of hanging: a resume that ignores its time budget or its stop
    would otherwise loop forever (the mutation check relies on this)."""
    def fail(signum, frame):
        raise TimeoutError("test exceeded its 60 s hard timeout")
    old = signal.signal(signal.SIGALRM, fail)
    signal.alarm(60)
    yield
    signal.alarm(0)
    signal.signal(signal.SIGALRM, old)


def advance(d: Path) -> None:
    """Advance a chain by one whole point, as a finished cycle does: .chain, then .param, then .trace."""
    n = coevol.saved_points(d, 1)
    with (d / "chain1.chain").open("a") as f:
        f.write("".join(f"pN l{k}\n" for k in range(LPP)))
    p = (d / "chain1.param").read_text()
    (d / "chain1.param").write_text(coevol.set_until(p, n, n + 5).replace(f"\t{n + 5}\t{n}", f"\t{n + 5}\t{n + 1}"))
    with (d / "chain1.trace").open("a") as f:
        f.write("n\t0.1\n")


def kill_after_append(d: Path) -> None:
    """A chain killed after appending its next point to .chain, before rewriting .param."""
    with (d / "chain1.chain").open("a") as f:
        f.write("".join(f"pX l{k}\n" for k in range(LPP)))


def kill_after_param(d: Path) -> None:
    """A chain killed after rewriting .param, before appending to .trace."""
    kill_after_append(d)
    n = coevol.saved_points(d, 1)
    p = (d / "chain1.param").read_text()
    (d / "chain1.param").write_text(coevol.set_until(p, n, n + 10).replace(f"\t{n + 10}\t{n}", f"\t{n + 10}\t{n + 1}"))


class TestParam:
    def test_size_is_the_published_chains(self):
        assert coevol.param_size(PARAM) == SIZE

    def test_set_until_changes_one_line_and_only_its_middle(self):
        new = coevol.set_until(PARAM, SIZE, 15000)
        a, b = PARAM.split("\n"), new.split("\n")
        diff = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
        assert len(a) == len(b) and len(diff) == 1
        assert a[diff[0]].split("\t")[::2] == b[diff[0]].split("\t")[::2]       # every and size kept
        assert b[diff[0]].split("\t")[1] == "15000"
        assert coevol.param_size(new) == SIZE

    def test_set_until_refuses(self):
        with pytest.raises(ValueError):
            coevol.set_until(PARAM, SIZE + 1, 15000)          # wrong size: no line matches
        with pytest.raises(ValueError):
            coevol.set_until(PARAM, SIZE, SIZE)               # nothing left to run
        assert coevol.param_size("no header\n") is None


class TestConsistency:
    def test_a_clean_chain_agrees(self, tmp_path):
        ok, why, st = coevol.consistency(make_chain(tmp_path), 1)
        assert ok, why
        assert st["size"] == st["param_size"] == SIZE and st["lines_per_point"] == LPP

    @pytest.mark.parametrize("kill", [kill_after_append, kill_after_param])
    def test_a_killed_chain_does_not(self, tmp_path, kill):
        d = make_chain(tmp_path)
        kill(d)
        ok, _, _ = coevol.consistency(d, 1)
        assert not ok

    def test_a_param_out_of_step_with_the_trace_alone_does_not(self, tmp_path):
        """A .param one point behind a .trace and .chain that agree (a stale .param put back by hand): only the
        comparison of the .param's size with the .trace sees it, since the .chain is a whole number of points."""
        d = make_chain(tmp_path)
        p = (d / "chain1.param").read_text()
        (d / "chain1.param").write_text(p.replace(f"\t{SIZE}\n", f"\t{SIZE - 1}\n", 1))
        assert coevol.param_size((d / "chain1.param").read_text()) == SIZE - 1
        ok, why, _ = coevol.consistency(d, 1)
        assert not ok and ".param" in why


class TestCheckpoints:
    @pytest.mark.parametrize("kill", [kill_after_append, kill_after_param])
    def test_restore_returns_the_checkpointed_bytes(self, tmp_path, kill):
        d = make_chain(tmp_path)
        before = {f: (d / f"chain1.{f}").read_bytes() for f in ("chain", "param", "trace", "monitor")}
        coevol.checkpoint(d, 1)
        kill(d)
        coevol.restore(d, 1)
        assert {f: (d / f"chain1.{f}").read_bytes() for f in before} == before
        assert coevol.consistency(d, 1)[0]

    def test_no_checkpoint_of_an_inconsistent_chain(self, tmp_path):
        d = make_chain(tmp_path)
        kill_after_append(d)
        with pytest.raises(RuntimeError):
            coevol.checkpoint(d, 1)

    def test_restore_refuses_without_a_checkpoint_or_with_a_changed_prefix(self, tmp_path):
        d = make_chain(tmp_path)
        with pytest.raises(RuntimeError):
            coevol.restore(d, 1)
        coevol.checkpoint(d, 1)
        b = bytearray((d / "chain1.chain").read_bytes()); b[0:2] = b"q0"
        (d / "chain1.chain").write_bytes(bytes(b))
        with pytest.raises(RuntimeError):
            coevol.restore(d, 1)

    def test_restore_takes_the_newest_checkpoint(self, tmp_path):
        """With two checkpoints, a killed chain goes back to the later one: the older is also a prefix of the
        append-only .chain, so restoring it would succeed and silently drop a segment of points."""
        d = make_chain(tmp_path)
        coevol.checkpoint(d, 1)
        advance(d)
        coevol.checkpoint(d, 1)
        newest = {f: (d / f"chain1.{f}").read_bytes() for f in ("chain", "param", "trace", "monitor")}
        kill_after_param(d)
        assert coevol.restore(d, 1).name == f"{SIZE + 1:07d}"
        assert {f: (d / f"chain1.{f}").read_bytes() for f in newest} == newest
        assert coevol.saved_points(d, 1) == SIZE + 1

    def test_keeps_the_newest_three(self, tmp_path):
        d = make_chain(tmp_path)
        for _ in range(5):
            coevol.checkpoint(d, 1)
            advance(d)
        kept = sorted(p.name for p in (d / "checkpoints" / "chain1").iterdir())
        assert kept == [f"{SIZE + k:07d}" for k in (2, 3, 4)]


# A stand-in for the coevol binary that writes its three files in Coevol's order, one point per cycle, until
# .run holds 0 or size passes until: enough to exercise the segments, the stop and the checkpoints for real.
# FAKE_MODE=hang: then waits forever with its files consistent, as rank 0 of Coevol's MPI build does.
# FAKE_MODE=stuck: after a few points, appends one more to .chain and waits forever there (files inconsistent).
FAKE = f"""#!{sys.executable}
import os, sys, time
from pathlib import Path
(Path('.') / "fake.pid").write_text(str(os.getpid()))
mode = os.environ.get("FAKE_MODE", "")
sys.path.insert(0, {str(Path(coevol.__file__).parents[2])!r})
from mitoage_pipelines.longevity import coevol
name = sys.argv[1]; d = Path('.')
(d / f"{{name}}.run").write_text("1\\n")
while (d / f"{{name}}.run").read_text().strip() != "0":
    p = (d / f"{{name}}.param").read_text(); n = coevol.param_size(p)
    every, until, _ = [l for l in p.split("\\n") if l.split("\\t")[-1:] == [str(n)] and len(l.split("\\t")) == 3][-1].split("\\t")
    if n > int(until):
        break
    time.sleep(0.05)
    with (d / f"{{name}}.chain").open("a") as f:
        f.write("".join(f"pF l{{k}}\\n" for k in range({LPP})))
    if mode == "stuck" and n >= {SIZE} + 2:
        time.sleep(3600)
    (d / f"{{name}}.param").write_text(p.replace(f"\\t{{until}}\\t{{n}}", f"\\t{{until}}\\t{{n + 1}}"))
    with (d / f"{{name}}.trace").open("a") as f:
        f.write("f\\t0.1\\n")
(d / f"{{name}}.run").write_text("0\\n")
if mode == "hang":
    time.sleep(3600)
"""

FAST = dict(poll_s=0.02, grace_s=0.4)


class TestResume:
    def fake(self, tmp_path, monkeypatch):
        exe = tmp_path / "coevol"
        exe.write_text(FAKE)
        exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
        monkeypatch.setenv("COEVOL", str(exe))

    def test_runs_to_until_with_a_checkpoint_per_segment(self, tmp_path, monkeypatch):
        self.fake(tmp_path, monkeypatch)
        d = make_chain(tmp_path / "out")
        coevol.resume(d, 1, until=SIZE + 6, segment_hours=0.5 / 3600, **FAST)
        ok, why, st = coevol.consistency(d, 1)
        assert ok, why
        assert st["size"] == SIZE + 7                  # Coevol stops once size > until
        cps = sorted(int(p.name) for p in (d / "checkpoints" / "chain1").iterdir())
        assert cps[-1] == SIZE + 7 and len(cps) >= 2

    def test_restores_a_killed_chain_before_running(self, tmp_path, monkeypatch, capsys):
        self.fake(tmp_path, monkeypatch)
        d = make_chain(tmp_path / "out")
        coevol.checkpoint(d, 1)
        kill_after_param(d)
        coevol.resume(d, 1, until=SIZE + 2, segment_hours=1 / 3600, **FAST)
        assert "RESTORED chain1" in capsys.readouterr().out
        ok, why, st = coevol.consistency(d, 1)
        assert ok, why
        assert st["size"] == SIZE + 3

    def test_refuses_a_killed_chain_without_a_checkpoint(self, tmp_path, monkeypatch):
        self.fake(tmp_path, monkeypatch)
        d = make_chain(tmp_path / "out")
        kill_after_append(d)
        with pytest.raises(RuntimeError):
            coevol.resume(d, 1, until=SIZE + 2, **FAST)

    def test_stops_when_its_time_is_spent(self, tmp_path, monkeypatch):
        self.fake(tmp_path, monkeypatch)
        d = make_chain(tmp_path / "out")
        coevol.resume(d, 1, until=10**6, stop_after_hours=1.5 / 3600, segment_hours=1 / 3600,
                      min_segment_hours=0.2 / 3600, **FAST)
        ok, why, st = coevol.consistency(d, 1)
        assert ok, why
        assert SIZE < st["size"] < SIZE + 200

    def test_ends_a_run_that_stops_saving_but_does_not_exit(self, tmp_path, monkeypatch, capsys):
        """Coevol's MPI build: rank 0 saves its last point and the job never exits. Each segment is ended once its
        files are consistent and unchanged, checkpointed, and the next one starts, up to `until`."""
        self.fake(tmp_path, monkeypatch)
        monkeypatch.setenv("FAKE_MODE", "hang")
        d = make_chain(tmp_path / "out")
        coevol.resume(d, 1, until=SIZE + 4, segment_hours=0.3 / 3600, **FAST)
        assert "ENDED a coevol" in capsys.readouterr().out
        ok, why, st = coevol.consistency(d, 1)
        assert ok, why
        assert st["size"] == SIZE + 5
        assert int(sorted((d / "checkpoints" / "chain1").iterdir())[-1].name) == SIZE + 5

    def test_does_not_end_a_run_caught_between_its_writes(self, tmp_path, monkeypatch):
        """A run whose files disagree (stuck between appending to .chain and rewriting .param) is never ended by
        the guard: killing it there would leave the chain to be restored; the job's time limit does that."""
        import threading
        self.fake(tmp_path, monkeypatch)
        monkeypatch.setenv("FAKE_MODE", "stuck")
        d = make_chain(tmp_path / "out")
        (d / "chain1.param").write_text(coevol.set_until(PARAM, SIZE, SIZE + 100))
        exe = os.environ["COEVOL"]
        done, failed = [], []

        def segment():
            try:
                done.append(coevol.run_segment([exe, "chain1"], d, 1, 0.5 / 3600, SIZE + 100, grace_s=0.3, poll_s=0.02))
            except subprocess.CalledProcessError as e:            # the kill below
                failed.append(e)
        th = threading.Thread(target=segment, daemon=True)
        th.start()
        time.sleep(2.5)                                       # asked to stop at 0.5 s, grace 0.3 s
        assert th.is_alive() and not done, (done, (d / "chain1.log").read_text()[-500:])  # still waiting, not ended
        assert not coevol.consistency(d, 1)[0]
        os.kill(int((d / "fake.pid").read_text()), signal.SIGKILL)
        th.join(10)
        assert not th.is_alive() and failed and failed[0].returncode == -signal.SIGKILL
