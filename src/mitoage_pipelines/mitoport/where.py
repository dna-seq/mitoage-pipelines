"""Where the differences sit (exploratory, 2026-11-24; the maintainer: "issue is not only which compatibility method
is the most accurate but what is exactly happening").

Three tables, no test of a hypothesis and no rule:

where_variants.tsv       harmful against harmless human variants (H29's set) by what their site is: how many amino
                         acids the build's mammals hold there, buried or not, touching a nuclear subunit, touching
                         another mitochondrial subunit.
where_cybrid_sites.tsv   one row per residue at which a cybrid's donor differs from its host, with the same
                         description of the site and the BLOSUM62 score of the change.
where_cybrid_summary.tsv one row per (host, donor, complex): how many differences, how many of them at conserved
                         sites, and the measured activity beside.

A site is called conserved when the build's mammals hold at most CONSERVED different amino acids there, and buried
when its relative accessible surface in the complex is below BURIED. Both thresholds were chosen before the tables
were looked at and are not tuned.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices

from . import cybrid, variant_stab
from .complexes import COMPLEXES

CONSERVED = 2
BURIED = 0.2
FLAGS = ("conserved", "buried", "contact_nuclear", "contact_mt")


def site_table(positions: pd.DataFrame, data: Path, sites_dir=lambda c: COMPLEXES[c].work) -> pd.DataFrame:
    """One row per human position of the 13 proteins: `positions.tsv` of `mitoport cybrid` with `n_aa` (the mammals'
    amino acids at the site) and `rasa_complex` (from the complex's sites.tsv; empty where the structure has no
    residue). `conserved` and `buried` are unknown (NaN) where their number is."""
    var = {g: variant_stab.site_variability(data, g) for g in sorted(set(positions.gene))}
    t = positions.copy()
    t["n_aa"] = [var[g].get(p, np.nan) for g, p in zip(t.gene, t.human_pos)]
    rasa = {}
    for c in sorted(set(t["complex"])):
        s = pd.read_csv(sites_dir(c) / "sites.tsv", sep="\t", dtype={"chain": str})
        rasa.update({(r.gene, r.human_pos): r.rasa_complex for r in s.itertuples()})
    t["rasa_complex"] = [rasa.get((g, p), np.nan) for g, p in zip(t.gene, t.human_pos)]
    t["conserved"] = np.where(t.n_aa.isna(), np.nan, (t.n_aa <= CONSERVED).astype(float))
    t["buried"] = np.where(t.rasa_complex.isna(), np.nan, (t.rasa_complex < BURIED).astype(float))
    return t


def variants_where(variants: pd.DataFrame, sites: pd.DataFrame) -> pd.DataFrame:
    """For each flag and each label: how many variants have the flag known, and how many of them have it; then the
    variants by the number of amino acids at their site."""
    v = variants.merge(sites[["gene", "human_pos", "rasa_complex", *FLAGS]].rename(columns={"human_pos": "protein_pos"}),
                       on=["gene", "protein_pos"], how="left", validate="many_to_one")
    rows = []
    for flag in FLAGS:
        for label in ("harmful", "harmless"):
            d = v[(v.label == label) & v[flag].notna()]
            rows.append({"what": flag, "label": label, "known": len(d), "with": int(d[flag].sum()),
                         "share": d[flag].mean() if len(d) else np.nan})
    for lo, hi, name in ((1, 1, "1"), (2, 2, "2"), (3, 4, "3 to 4"), (5, 20, "5 or more")):
        for label in ("harmful", "harmless"):
            d = v[(v.label == label) & v.n_aa.notna()]
            k = int(((d.n_aa >= lo) & (d.n_aa <= hi)).sum())
            rows.append({"what": f"amino acids at the site: {name}", "label": label, "known": len(d), "with": k,
                         "share": k / len(d) if len(d) else np.nan})
    return pd.DataFrame(rows)


def cybrid_sites(sites: pd.DataFrame, human: dict[str, str], prot: dict, pairs: list[tuple[str, str]]) -> pd.DataFrame:
    """One row per (host, donor, gene, human position) at which both residues are known and differ, as
    `cybrid.pair_scores` compares them (the donor read through the host)."""
    blosum = substitution_matrices.load("BLOSUM62")
    al = cybrid.aligner()
    rows = []
    for h, d in pairs:
        for g, grp in sites.groupby("gene", sort=False):
            hp, hpart, _ = prot[h][g]; dp, dpart, _ = prot[d][g]
            hs = cybrid.on_human(human[g], hp, hpart, al); ds = cybrid.donor_through_host(human[g], hp, dp, dpart, al)
            for r in grp.itertuples():
                a, b = hs[r.human_pos - 1], ds[r.human_pos - 1]
                if "?" in (a, b) or a == b:
                    continue
                rows.append({"host": h, "donor": d, "complex": r.complex, "gene": g, "human_pos": r.human_pos,
                             "human_aa": r.human_aa, "host_aa": a, "donor_aa": b,
                             "blosum": np.nan if "-" in (a, b) else float(blosum[a][b]), "n_aa": r.n_aa,
                             "rasa_complex": r.rasa_complex, **{f: getattr(r, f) for f in FLAGS},
                             "nuclear_chains": r.nuclear_chains})
    return pd.DataFrame(rows)


def cybrid_summary(diff: pd.DataFrame, meas: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (h, d, c), m in meas.groupby(["host", "donor", "complex"], sort=True):
        x = diff[(diff.host == h) & (diff.donor == d) & (diff.complex == c)]
        cons = x[x.conserved == 1]
        rows.append({"host": h, "donor": d, "complex": c, "outcome": m.outcome.iloc[0], "pct": m.pct.mean(),
                     "significant": m.significant.iloc[0], "differences": len(x), "at_conserved": len(cons),
                     "at_conserved_buried": int((cons.buried == 1).sum()),
                     "at_conserved_contact_nuclear": int((cons.contact_nuclear == 1).sum()),
                     "at_conserved_radical": int((cons.blosum.isna() | (cons.blosum < 0)).sum()),
                     "at_contact_nuclear": int((x.contact_nuclear == 1).sum())})
    return pd.DataFrame(rows)


def main(variants: Path, positions: Path, data: Path, gb: Path, out: Path, measurements: Path = cybrid.MEASUREMENTS) -> None:
    out.mkdir(parents=True, exist_ok=True)
    sites = site_table(pd.read_csv(positions, sep="\t", keep_default_na=False, na_values=[""]), data)
    v = pd.read_csv(variants, sep="\t", keep_default_na=False, na_values=[""])
    vw = variants_where(v[["id", "gene", "protein_pos", "label", "n_aa"]], sites)
    meas = cybrid.read_measurements(measurements)
    meas = meas[meas["complex"].isin(list(COMPLEXES))]          # whole-cybrid rows (`all`) have no structure
    records = cybrid.read_proteins(gb)
    species = sorted(set(meas.host) | set(meas.donor) | {cybrid.HUMAN})
    prot = {s: cybrid.species_proteins(records, s) for s in species}
    human = {g: p for g, (p, _, _) in prot[cybrid.HUMAN].items()}
    diff = cybrid_sites(sites, human, prot, sorted(set(zip(meas.host, meas.donor))))
    summ = cybrid_summary(diff, meas)
    vw.to_csv(out / "where_variants.tsv", sep="\t", index=False, float_format="%.3f")
    diff.to_csv(out / "where_cybrid_sites.tsv", sep="\t", index=False, float_format="%.4f")
    summ.to_csv(out / "where_cybrid_summary.tsv", sep="\t", index=False, float_format="%.1f")
    with pd.option_context("display.width", 220, "display.max_rows", 200):
        print(vw.to_string(index=False)); print(summ.to_string(index=False))
    print("wrote", out / "where_variants.tsv", out / "where_cybrid_sites.tsv", out / "where_cybrid_summary.tsv")
