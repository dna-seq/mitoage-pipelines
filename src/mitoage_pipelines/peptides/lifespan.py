"""RM16's third question, second part (secondary): does keeping a frame follow maximum lifespan?

H12's model, unchanged (`aa_pgls.select`, `aa_pgls.class_covariance`, `aa_pgls.fit`): per class, H1's default
dataset; log10 maximum lifespan on the frame's state (1 kept, 0 lost; a species not assessed is left out, never
counted as lost), with log10 adult body mass and the fourfold composition (A4, C4, T4), Pagel's lambda by maximum
likelihood. The state is not standardised: its slope is the difference in log10 years between species that keep
the frame and species that lost it, given the covariates. The reading is each peptide's published one
(`own` for MTALTND4, MTALTCO1 and MENTSH, `standard` for the rest). A class needs MIN_SPECIES species with at least
MIN_EACH kept and MIN_EACH lost, or it is written as not assessed. q: Benjamini-Hochberg over the peptides, within
class.
Writes pgls.tsv.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from ..longevity import aa_pgls
from .catalogue import PEPTIDES

CLASSES = ("Mammalia", "Aves", "Actinopterygii", "Reptilia")
MIN_SPECIES = 25
MIN_EACH = 5


def primary_reading(published: str) -> str:
    return "standard" if published == "standard" else "own"


def main(data: Path, states_file: Path, composition: Path, fourfold: Path, tree_file: Path, out: Path,
         classes=CLASSES) -> pd.DataFrame:
    states = pd.read_csv(states_file, sep="\t")
    rows = []
    for cls in classes:
        sp, _ = aa_pgls.select(data, composition, fourfold, tree_file, cls)
        for pep in PEPTIDES:
            reading = primary_reading(pep.reading)
            st = states[(states.peptide == pep.name) & (states.reading == reading) & (states.state != "not assessed")]
            d = sp.merge(st[["slug", "state"]], on="slug", how="inner")
            d["kept"] = (d.state == "kept").astype(float)
            base = {"class": cls, "peptide": pep.name, "reading": reading, "n": len(d),
                    "n_kept": int(d.kept.sum()), "n_lost": int((d.kept == 0).sum())}
            if len(d) < MIN_SPECIES or base["n_kept"] < MIN_EACH or base["n_lost"] < MIN_EACH:
                rows.append({**base, "status": f"not assessed: under {MIN_SPECIES} species or {MIN_EACH} kept or lost"})
                continue
            C, order, _ = aa_pgls.class_covariance(tree_file, list(d.slug))
            d = d.set_index("slug").loc[order].reset_index()
            X = np.column_stack([d.kept.to_numpy(), d.lm.to_numpy(), d[["A4", "C4", "T4"]].to_numpy()])
            r = aa_pgls.fit(C, d.ll.to_numpy(), X)
            h = stats.t.ppf(0.975, r["df"]) * r["se"]
            rows.append({**base, "lambda": r["lambda"], "slope": r["slope"], "se": r["se"], "lo95": r["slope"] - h,
                         "hi95": r["slope"] + h, "t": r["t"], "p": r["p"], "status": "tested"})
    res = pd.DataFrame(rows)
    res["q"] = np.nan
    for cls, g in res[res.status == "tested"].groupby("class"):
        res.loc[g.index, "q"] = multipletests(g.p, method="fdr_bh")[1]
    res = res.sort_values(["class", "peptide"], kind="stable")
    out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, sep="\t", index=False, float_format="%.6g")
    return res
