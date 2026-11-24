"""Which nuclear subunits and residues face the donor's mitochondrial differences in the cybrids that lost activity
(exploratory, 2026-11-24; the maintainer: the aim is to find the nuclear genes to edit, since only the nucleus can be
edited).

Inside one host the nucleus is the same in every cybrid, so a nuclear residue cannot be ranked by its own variation.
What does vary is the donor's mitochondrial protein, and with it which nuclear residues have a changed partner. For
each (host, donor, complex) with a measured outcome, a nuclear residue "faces a difference" when a mitochondrial
residue that touches it (the contact table of the structure) differs between donor and host
(`where.cybrid_sites`, which reads the donor through the host).

nuclear_subunits.tsv  one row per (complex, nuclear gene): pairs reduced and not, the mean number of residues facing a
                      difference in each, the area under the ROC curve of that number for telling them apart, the
                      same area when each donor in turn is left out (range), the mean of that area when the differing
                      sites are drawn at random from the same protein (the interface's own expectation), the share of
                      draws at or above the observed area, and the share of draws in which the best subunit of the
                      complex reaches the observed best (the family-wise p); the subunits with identical columns in
                      the table of pairs share a `profile_group`.
nuclear_exposure.tsv  the table of pairs by subunits behind it.
nuclear_residues.tsv  one row per nuclear residue that faces a difference in at least one pair: in how many reduced
                      and how many not reduced pairs, and the difference of the two shares (descriptive only).

A subunit is named only when its complex's family-wise p is below 0.05 (docs/HYPOTHESES.md, H36, the change of
2026-11-24). There are 9 reduced and 20 not reduced pairs in all, and a few per complex.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import cybrid
from .complexes import COMPLEXES


def labelled(meas: pd.DataFrame, disagree: str = "leave_out") -> pd.DataFrame:
    """One row per (host, donor, complex) with `status` (reduced / not_reduced) and the mean stated activity `pct`
    (NaN where none is stated). Pairs the papers disagree on are left out (`disagree="leave_out"`) or coded reduced
    (`"reduced"`)."""
    rows = []
    for (h, d, c), m in meas.groupby(["host", "donor", "complex"], sort=True):
        s = cybrid.status(m)
        if s == "disagree":
            if disagree == "leave_out":
                continue
            s = "reduced"
        rows.append({"host": h, "donor": d, "complex": c, "status": s, "pct": m.pct.mean()})
    return pd.DataFrame(rows)


def facing(diff: pd.DataFrame, contacts: pd.DataFrame, chains: pd.DataFrame) -> pd.DataFrame:
    """The nuclear residues that face a difference: one row per (host, donor, complex, nuclear gene, residue number,
    residue), from `where.cybrid_sites` rows and the structure's contact table of the complex (nuclear partners only)."""
    gene = dict(zip(chains.chain.astype(str), chains.gene))
    c = contacts[contacts.partner_is_mt == 0].copy()
    c["nuclear_gene"] = c.partner_chain.astype(str).map(gene)
    c = c.rename(columns={"mt_gene": "gene"})[["gene", "human_pos", "nuclear_gene", "partner_resnum", "partner_res"]].drop_duplicates()
    out = diff[["host", "donor", "complex", "gene", "human_pos"]].merge(c, on=["gene", "human_pos"], how="inner")
    return out.drop_duplicates(["host", "donor", "complex", "nuclear_gene", "partner_resnum"])


Pair = tuple[str, str, str]
Site = tuple[str, int]


def site_masks(contacts: pd.DataFrame, chains: pd.DataFrame) -> dict[Site, dict[str, int]]:
    """{(mitochondrial gene, human position): {nuclear gene: bit mask of that subunit's residues touching it}}; bit k
    is the subunit's k-th contacting residue number in order."""
    gene = dict(zip(chains.chain.astype(str), chains.gene))
    c = contacts[contacts.partner_is_mt == 0].copy()
    c["nuc"] = c.partner_chain.astype(str).map(gene)
    c = c.dropna(subset=["nuc"])
    bit = {n: {r: k for k, r in enumerate(sorted(set(g.partner_resnum)))} for n, g in c.groupby("nuc")}
    out: dict[Site, dict[str, int]] = {}
    for r in c.itertuples():
        d = out.setdefault((r.mt_gene, int(r.human_pos)), {})
        d[r.nuc] = d.get(r.nuc, 0) | (1 << bit[r.nuc][r.partner_resnum])
    return out


def pair_counts(sites, masks: dict[Site, dict[str, int]]) -> dict[str, int]:
    """{nuclear gene: how many of its residues touch at least one of the sites}."""
    acc: dict[str, int] = {}
    for site in sites:
        for u, m in masks.get(site, {}).items():
            acc[u] = acc.get(u, 0) | m
    return {u: m.bit_count() for u, m in acc.items()}


def auc_columns(counts: np.ndarray, red: np.ndarray) -> np.ndarray:
    """The area under the ROC curve of each column of `counts` (pairs x subunits) for telling the reduced rows from the
    others (ties count one half); NaN for every column when either group is empty. Equal to `cybrid.auc`."""
    if red.sum() == 0 or (~red).sum() == 0:
        return np.full(counts.shape[1], np.nan)
    a, b = counts[red][:, None, :], counts[~red][None, :, :]
    return ((a > b) + 0.5 * (a == b)).mean(axis=(0, 1))


def draw_sites(sites: dict[Pair, list[Site]], universe: dict[str, np.ndarray], rng: np.random.Generator) -> dict[Pair, list[Site]]:
    """Each pair's differing sites replaced by as many positions of the same protein, drawn without replacement from
    its resolved positions."""
    out = {}
    for k, ss in sites.items():
        n: dict[str, int] = {}
        for g, _ in ss:
            n[g] = n.get(g, 0) + 1
        out[k] = [(g, int(p)) for g, m in n.items() for p in rng.choice(universe[g], size=m, replace=False)]
    return out


def subunits(lab: pd.DataFrame, diff: pd.DataFrame, masks_by_complex: dict[str, dict[Site, dict[str, int]]],
             universe_by_complex: dict[str, dict[str, np.ndarray]], draws: int = 1000, seed: int = 0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(one row per (complex, nuclear gene), the table of pairs by subunits). The differing sites of a pair are those
    of `diff` that lie in the complex's resolved positions."""
    key = ["host", "donor", "complex"]
    rows, expo = [], []
    for cx, grp in lab.groupby("complex"):
        masks, universe = masks_by_complex[cx], universe_by_complex[cx]
        inside = {g: set(map(int, v)) for g, v in universe.items()}
        pairs = [tuple(r) for r in grp[key].itertuples(index=False)]
        sites = {k: [(g, int(p)) for g, p in zip(x.gene, x.human_pos) if int(p) in inside.get(g, ())]
                 for k, x in ((k, diff[(diff.host == k[0]) & (diff.donor == k[1]) & (diff["complex"] == k[2])]) for k in pairs)}
        subs = sorted({u for m in masks.values() for u in m})
        red = (grp.status == "reduced").to_numpy()
        donors = grp.donor.to_numpy()

        def matrix(ss):
            return np.array([[pair_counts(ss[k], masks).get(u, 0) for u in subs] for k in pairs], dtype=float)
        obs = matrix(sites)
        auc = auc_columns(obs, red)
        lodo = np.array([auc_columns(obs[donors != d], red[donors != d]) for d in sorted(set(donors))])
        rng = np.random.default_rng(seed)
        null = np.empty((draws, len(subs)))
        for i in range(draws):
            null[i] = auc_columns(matrix(draw_sites(sites, universe, rng)), red)
        best = np.nanmax(auc) if not np.isnan(auc).all() else np.nan
        fw = (1 + np.sum(np.nanmax(null, axis=1) >= best)) / (1 + draws) if not np.isnan(best) else np.nan
        same = {}
        for k, u in enumerate(subs):
            same.setdefault(tuple(obs[:, k]), []).append(u)
        group = {u: (i, len(v)) for i, v in enumerate(same.values()) for u in v}
        for k, u in enumerate(subs):
            rows.append({"complex": cx, "nuclear_gene": u, "pairs_reduced": int(red.sum()), "pairs_not_reduced": int((~red).sum()),
                         "mean_facing_reduced": obs[red, k].mean() if red.any() else np.nan,
                         "mean_facing_not_reduced": obs[~red, k].mean() if (~red).any() else np.nan,
                         "auc": auc[k], "auc_leave_donor_min": np.nanmin(lodo[:, k]) if not np.isnan(lodo[:, k]).all() else np.nan,
                         "auc_leave_donor_max": np.nanmax(lodo[:, k]) if not np.isnan(lodo[:, k]).all() else np.nan,
                         "auc_null_mean": np.nanmean(null[:, k]),
                         "p_subunit": (1 + np.sum(null[:, k] >= auc[k])) / (1 + draws) if not np.isnan(auc[k]) else np.nan,
                         "p_family": fw, "named": bool(fw < 0.05 and auc[k] == best) if not np.isnan(fw) else False,
                         "profile_group": f"{cx}{group[u][0]}", "profile_group_size": group[u][1]})
        for r, k in zip(range(len(pairs)), pairs):
            expo.append({"host": k[0], "donor": k[1], "complex": k[2], "status": grp.status.iloc[r],
                         **{u: int(obs[r, j]) for j, u in enumerate(subs)}})
    return (pd.DataFrame(rows).sort_values(["complex", "auc"], ascending=[True, False]).reset_index(drop=True), pd.DataFrame(expo))


def residues(lab: pd.DataFrame, face: pd.DataFrame) -> pd.DataFrame:
    key = ["host", "donor", "complex"]
    f = face.merge(lab[key + ["status"]], on=key, how="inner")
    n_red = lab[lab.status == "reduced"].groupby("complex").size()
    n_not = lab[lab.status == "not_reduced"].groupby("complex").size()
    g = (f.groupby(["complex", "nuclear_gene", "partner_resnum", "partner_res", "status"]).size().unstack("status", fill_value=0)
         .reindex(columns=["reduced", "not_reduced"], fill_value=0).reset_index())
    g["share_reduced"] = [r.reduced / n_red.get(r.complex, np.nan) for r in g.itertuples()]
    g["share_not_reduced"] = [r.not_reduced / n_not.get(r.complex, np.nan) for r in g.itertuples()]
    g["difference"] = g.share_reduced - g.share_not_reduced
    return g.sort_values(["difference", "reduced"], ascending=False).reset_index(drop=True)


def main(sites: Path, positions: Path, out: Path, measurements: Path = cybrid.MEASUREMENTS, draws: int = 1000) -> None:
    diff = pd.read_csv(sites, sep="\t", keep_default_na=False, na_values=[""])
    pos = pd.read_csv(positions, sep="\t", keep_default_na=False, na_values=[""])
    meas = cybrid.read_measurements(measurements)
    lab = labelled(meas)
    masks, universe, faces = {}, {}, []
    for cx in sorted(set(meas["complex"])):
        w = COMPLEXES[cx].work
        contacts = pd.read_csv(w / "contacts.tsv", sep="\t", dtype={"chain": str, "partner_chain": str})
        chains = pd.read_csv(w / "chains.tsv", sep="\t", dtype={"chain": str})
        masks[cx] = site_masks(contacts, chains)
        r = pos[(pos["complex"] == cx) & (pos.resolved == 1)]
        universe[cx] = {g: np.sort(x.human_pos.to_numpy()) for g, x in r.groupby("gene")}
        faces.append(facing(diff[diff["complex"] == cx], contacts, chains))
    face = pd.concat(faces, ignore_index=True)
    s, expo = subunits(lab, diff, masks, universe, draws)
    s2, _ = subunits(labelled(meas, "reduced"), diff, masks, universe, draws)
    s = s.merge(s2[["complex", "nuclear_gene", "named", "auc", "p_family"]].rename(
        columns={"named": "named_if_disagree_reduced", "auc": "auc_if_disagree_reduced", "p_family": "p_family_if_disagree_reduced"}),
        on=["complex", "nuclear_gene"], how="left")
    s["stable"] = s.named == s.named_if_disagree_reduced
    r = residues(lab, face)
    s.to_csv(out / "nuclear_subunits.tsv", sep="\t", index=False, float_format="%.3f")
    expo.to_csv(out / "nuclear_exposure.tsv", sep="\t", index=False)
    r.to_csv(out / "nuclear_residues.tsv", sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.width", 250, "display.max_rows", 120, "display.max_columns", 30):
        print(s.drop(columns=["pairs_reduced", "pairs_not_reduced"]).to_string(index=False)); print(r.head(25).to_string(index=False))
    print("wrote", out / "nuclear_subunits.tsv", out / "nuclear_exposure.tsv", out / "nuclear_residues.tsv")
