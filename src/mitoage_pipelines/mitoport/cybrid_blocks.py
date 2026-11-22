"""H32: are the conserved sequence blocks of the control region special among the cybrids, or one more ruler of
distance? The blocks' differences between host and donor against random sets of windows of the same lengths.

The rule is in docs/HYPOTHESES.md (H32) and was written before this code: the lead counts only if fewer than 5% of
random window sets in the control region reach the blocks' area (reduced against not reduced) and fewer than 5% its
Spearman with the stated activities, on the mouse-host points.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from scipy import stats

from . import cybrid, cybrid_dna

MOUSE = "mus-musculus"
SETS = 10_000
SHARE_MAX = 0.05
LOCAL_MIN = 0.5            # a block placed by local alignment must score at least this share of its length


def difference_vector(host: str, donor: str, al: PairwiseAligner) -> np.ndarray:
    """One number per place of `host`: 1 where the donor's base differs or is missing in the global alignment, plus
    1 for each base the donor has inserted right after that place (before the first place: on the first)."""
    out = np.zeros(len(host), dtype=np.int64)
    i = -1
    for x, y in cybrid_dna.columns(host, donor, al):
        if x != "-":
            i += 1
            out[i] += x != y
        else:
            out[max(i, 0)] += 1
    return out


def local_aligner() -> PairwiseAligner:
    al = PairwiseAligner()
    al.mode = "local"
    al.match_score, al.mismatch_score = 1, -1
    al.open_gap_score, al.extend_gap_score = -3, -1
    return al


def blocks_of(features: dict, host: str) -> dict[str, tuple[int, int]]:
    """{block: (start, end)} as places of the host's control region (0-based, end excluded): the host record's own
    annotation where it has one, else the mouse record's blocks placed in the host's control region by local
    alignment. Refuses a placement that scores under LOCAL_MIN of the block's length."""
    h = features[cybrid.GENOMES[host]]
    if h["dloop"] is None:
        raise ValueError(f"{host}: the record annotates no control region")
    seq, start = h["dloop"]
    if set(h["blocks"]) == set(cybrid_dna.CSB):
        return {k: (a - start, b - start) for k, (a, b) in sorted(h["blocks"].items())}
    m = features[cybrid.GENOMES[MOUSE]]
    mseq, mstart = m["dloop"]
    al, out = local_aligner(), {}
    for k, (a, b) in sorted(m["blocks"].items()):
        block = mseq[a - mstart:b - mstart]
        best = al.align(seq, block)[0]
        if best.score < LOCAL_MIN * len(block):
            raise ValueError(f"{host}: {k} of the mouse record finds no place (score {best.score} for {len(block)} bases)")
        out[k] = (int(best.aligned[0][0][0]), int(best.aligned[0][-1][1]))
    return out


def random_sets(length: int, sizes: list[int], n: int, rng: np.random.Generator) -> np.ndarray:
    """[n, len(sizes)] starts of windows of the given sizes inside 0..length, the windows of one set not overlapping
    (a set with an overlap is drawn again)."""
    out = np.empty((n, len(sizes)), dtype=np.int64)
    k = 0
    while k < n:
        starts = [int(rng.integers(0, length - s + 1)) for s in sizes]
        spans = sorted(zip(starts, sizes))
        if all(a + s <= b for (a, s), (b, _) in zip(spans, spans[1:])):
            out[k] = starts; k += 1
    return out


def window_sums(vector: np.ndarray, starts: np.ndarray, sizes: list[int]) -> np.ndarray:
    """[n] the sum of the vector over each set's windows."""
    c = np.concatenate([[0], np.cumsum(vector)])
    return sum(c[starts[:, j] + s] - c[starts[:, j]] for j, s in enumerate(sizes))


def statistics(score_by_donor: dict[str, np.ndarray], points: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """(area for reduced against not reduced, Spearman with the stated activity), each [n], for n scores per donor
    laid on the measured points (one row per measured complex: donor, status, pct). Ties count one half in the
    area; a constant score gives an area of 0.5 and a Spearman of NaN."""
    use = points[points.status != "disagree"]
    s = np.stack([score_by_donor[d] for d in use.donor]).astype(float)              # [points, n]
    red = (use.status == "reduced").to_numpy()
    r = stats.rankdata(s, axis=0)
    n1 = int(red.sum())
    area = (r[red].sum(0) - n1 * (n1 + 1) / 2) / (n1 * (len(red) - n1))
    st = points[points.pct.notna()]
    t = np.stack([score_by_donor[d] for d in st.donor]).astype(float)
    rt = stats.rankdata(t, axis=0); rp = stats.rankdata(st.pct.to_numpy())
    rt = rt - rt.mean(0); rp = rp - rp.mean()
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = (rt * rp[:, None]).sum(0) / np.sqrt((rt ** 2).sum(0) * (rp ** 2).sum())
    return area, rho


def measured_points(meas: pd.DataFrame, host: str) -> pd.DataFrame:
    """One row per measured (donor, complex) of the host: status (H16's) and the mean stated activity."""
    m = meas[(meas.outcome == "measured") & (meas.host == host)]
    return (m.groupby(["donor", "complex"], sort=False)
            .apply(lambda g: pd.Series({"status": cybrid.status(g), "pct": g.pct.mean()}), include_groups=False).reset_index())


def tracks(features: dict, host: str, donor: str, al: PairwiseAligner) -> dict[str, np.ndarray]:
    """The pair's difference vectors: the control region, the two rRNA genes end to end, the 13 protein genes end to
    end (each gene aligned on its own)."""
    h, d = features[cybrid.GENOMES[host]], features[cybrid.GENOMES[donor]]
    return {"control region": difference_vector(h["dloop"][0], d["dloop"][0], al),
            "rRNA genes": np.concatenate([difference_vector(h["rrna"][k][0], d["rrna"][k][0], al) for k in ("12S", "16S")]),
            "protein genes": np.concatenate([difference_vector(cybrid_dna.cds_of(features, host, g),
                                                               cybrid_dna.cds_of(features, donor, g), al)
                                             for g in cybrid.ALL_GENES])}


def main(gb: Path, out: Path, measurements: Path = cybrid.MEASUREMENTS, sets: int = SETS, seed: int = 0) -> pd.DataFrame:
    """Writes h32_blocks.tsv (the blocks' places and each donor's differences in them) and h32_summary.tsv."""
    out.mkdir(parents=True, exist_ok=True)
    features = cybrid_dna.read_features(gb)
    meas = cybrid.read_measurements(measurements)
    al = cybrid_dna.aligner()
    rng = np.random.default_rng(seed)
    rows, per_donor = [], []
    for host in sorted(set(meas.host)):
        blocks = blocks_of(features, host)
        sizes = [b - a for a, b in blocks.values()]
        starts0 = np.array([[a for a, _ in blocks.values()]])
        donors = sorted(set(meas.donor[meas.host == host]))
        skipped = [d for d in donors if features[cybrid.GENOMES[d]]["dloop"] is None]
        donors = [d for d in donors if d not in skipped]
        tr = {d: tracks(features, host, d, al) for d in donors}
        own = {d: window_sums(tr[d]["control region"], starts0, sizes) for d in donors}
        for d in donors:
            per_donor.append({"host": host, "donor": d, "blocks": int(own[d][0]),
                              **{k: int(tr[d]["control region"][a:b].sum()) for k, (a, b) in blocks.items()},
                              "control_region": int(tr[d]["control region"].sum()),
                              "places": ";".join(f"{k}:{a + 1}-{b}" for k, (a, b) in blocks.items())})
        for d in skipped:
            rows.append({"host": host, "check": "not assessed", "windows": "", "n": 0, "blocks": np.nan, "share": np.nan,
                         "reads": f"{d}: the record annotates no control region"})
        null = {}
        for kind in ("control region", "rRNA genes", "protein genes"):
            starts = random_sets(len(next(iter(tr.values()))[kind]), sizes, sets, rng)
            null[kind] = {d: window_sums(tr[d][kind], starts, sizes) for d in donors}
        points = measured_points(meas, host)
        points = points[points.donor.isin(donors)]
        if (points.status == "reduced").any() and (points.status == "not_reduced").any() and points.donor.nunique() >= 4:
            a0, r0 = (x[0] for x in statistics(own, points))
            passed = []
            for kind, sc in null.items():
                a, r = statistics(sc, points)
                sa, sr = float((a >= a0).mean()), float((r <= r0).mean())
                if kind == "control region":
                    passed = [sa < SHARE_MAX, sr < SHARE_MAX]
                rows.append({"host": host, "check": "1: area, reduced against not", "windows": kind, "n": sets, "blocks": a0, "share": sa,
                             "reads": f"share of random sets with an area of at least the blocks'; {int((points.status == 'reduced').sum())} "
                                      f"reduced, {int((points.status == 'not_reduced').sum())} not, {points.donor.nunique()} donors"})
                rows.append({"host": host, "check": "1: Spearman with activity", "windows": kind, "n": sets, "blocks": r0, "share": sr,
                             "reads": f"share of random sets with a Spearman at most the blocks'; {int(points.pct.notna().sum())} stated activities"})
            rows.append({"host": host, "check": "1: rule", "windows": "control region", "n": sets, "blocks": np.nan, "share": np.nan,
                         "reads": ("counts" if all(passed) else "does not count") + f": both shares under {SHARE_MAX}"})
        else:
            rows.append({"host": host, "check": "1: rule", "windows": "control region", "n": 0, "blocks": np.nan, "share": np.nan,
                         "reads": f"not assessed: {points.donor.nunique()} measured donors, each reduced in the same complexes or none; "
                                  "a per-donor score cannot separate them"})
        failed = sorted(set(meas.donor[(meas.host == host) & (meas.outcome == "no_respiring_cybrid")]) & set(donors))
        respired = sorted(set(points.donor))
        for f in failed:
            above = lambda sc: sc[f] > np.max([sc[d] for d in respired], axis=0)
            for kind, sc in null.items():
                rows.append({"host": host, "check": f"2: {f} above every donor that respired", "windows": kind, "n": sets,
                             "blocks": float(above(own)[0]), "share": float(above(sc).mean()),
                             "reads": "blocks: 1 if so; share of random sets where it is so; reported, claimed as nothing"})
    blocks_table = pd.DataFrame(per_donor)
    summary = pd.DataFrame(rows)
    blocks_table.to_csv(out / "h32_blocks.tsv", sep="\t", index=False)
    summary.to_csv(out / "h32_summary.tsv", sep="\t", index=False, float_format="%.4f")
    with pd.option_context("display.width", 230, "display.max_colwidth", 120, "display.max_rows", 200):
        print(blocks_table.to_string(index=False)); print(summary.to_string(index=False))
    print("wrote", out / "h32_blocks.tsv", out / "h32_summary.tsv")
    return summary
