"""Which nuclear subunits and residues face the donor's mitochondrial differences in the cybrids that lost activity
(exploratory, 2026-11-24; the maintainer: the aim is to find the nuclear genes to edit, since only the nucleus can be
edited).

Inside one host the nucleus is the same in every cybrid, so a nuclear residue cannot be ranked by its own variation.
What does vary is the donor's mitochondrial protein, and with it which nuclear residues have a changed partner. For
each (host, donor, complex) with a measured outcome, a nuclear residue "faces a difference" when a mitochondrial
residue that touches it (the contact table of the structure) differs between donor and host
(`where.cybrid_sites`, which reads the donor through the host).

nuclear_subunits.tsv  one row per (complex, nuclear gene): how many pairs reduced and not, the mean number of
                      residues facing a difference in each, the area under the ROC curve of that number for telling
                      reduced from not reduced, the Spearman correlation with the stated activity, and the same two
                      numbers for the complex's whole count of differences, to say whether the subunit adds anything.
nuclear_residues.tsv  one row per nuclear residue that faces a difference in at least one pair: in how many reduced
                      and how many not reduced pairs, and the difference of the two shares.

No rule, no claim: there are 9 reduced and 20 not reduced pairs in all, and a few per complex.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import cybrid
from .complexes import COMPLEXES


def labelled(meas: pd.DataFrame) -> pd.DataFrame:
    """One row per (host, donor, complex) with `status` (reduced / not_reduced / disagree) and the mean stated activity
    `pct` (NaN where none is stated); pairs the papers disagree on are left out."""
    rows = []
    for (h, d, c), m in meas.groupby(["host", "donor", "complex"], sort=True):
        s = cybrid.status(m)
        if s != "disagree":
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


def subunits(lab: pd.DataFrame, diff: pd.DataFrame, face: pd.DataFrame) -> pd.DataFrame:
    rows = []
    key = ["host", "donor", "complex"]
    total = diff.groupby(key).size().rename("all").reset_index()
    for cx, grp in lab.groupby("complex"):
        g = grp.merge(total, on=key, how="left").fillna({"all": 0})
        red = (g.status == "reduced").to_numpy()
        for nuc in sorted(set(face[face["complex"] == cx].nuclear_gene.dropna())):
            n = face[(face["complex"] == cx) & (face.nuclear_gene == nuc)].groupby(key).size().rename("n").reset_index()
            x = g.merge(n, on=key, how="left").fillna({"n": 0})
            sc = pd.Series(x.n.to_numpy(float)); ref = pd.Series(x["all"].to_numpy(float))
            ok = x.pct.notna().to_numpy()
            rho = lambda s: spearmanr(s[ok], x.pct[ok])[0] if ok.sum() >= 5 and s[ok].nunique() > 1 else np.nan
            rows.append({"complex": cx, "nuclear_gene": nuc, "pairs_reduced": int(red.sum()), "pairs_not_reduced": int((~red).sum()),
                         "mean_facing_reduced": sc[red].mean(), "mean_facing_not_reduced": sc[~red].mean(),
                         "auc": cybrid.auc(sc, pd.Series(red)), "spearman_pct": rho(sc), "pairs_with_pct": int(ok.sum()),
                         "auc_complex_count": cybrid.auc(ref, pd.Series(red)), "spearman_pct_complex_count": rho(ref)})
    return pd.DataFrame(rows).sort_values(["complex", "auc"], ascending=[True, False]).reset_index(drop=True)


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


def main(sites: Path, out: Path, measurements: Path = cybrid.MEASUREMENTS) -> None:
    diff = pd.read_csv(sites, sep="\t", keep_default_na=False, na_values=[""])
    lab = labelled(cybrid.read_measurements(measurements))
    faces = []
    for cx in sorted(set(lab["complex"])):
        w = COMPLEXES[cx].work
        faces.append(facing(diff[diff["complex"] == cx], pd.read_csv(w / "contacts.tsv", sep="\t", dtype={"chain": str, "partner_chain": str}),
                            pd.read_csv(w / "chains.tsv", sep="\t", dtype={"chain": str})))
    face = pd.concat(faces, ignore_index=True)
    s, r = subunits(lab, diff, face), residues(lab, face)
    s.to_csv(out / "nuclear_subunits.tsv", sep="\t", index=False, float_format="%.3f")
    r.to_csv(out / "nuclear_residues.tsv", sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.width", 220, "display.max_rows", 80):
        print(s.to_string(index=False)); print(r.head(25).to_string(index=False))
    print("wrote", out / "nuclear_subunits.tsv", out / "nuclear_residues.tsv")
