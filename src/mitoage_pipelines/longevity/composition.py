"""Do the amino-acid effects in mammals survive base composition at fourfold sites?

Reads composition.tsv (from `controls`) and writes controls_composition.tsv.
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from .fourfold import fourfold_counts


def main(data: Path, composition: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    cnt = fourfold_counts(data)
    sp = pd.read_csv(composition, sep="\t")
    for b in "ACT":
        sp[b + "4"] = sp.slug.map(lambda s: cnt[s][b] / max(1, sum(cnt[s].values())))
    m = sp[(sp["class"] == "Mammalia")].dropna(subset=["lm"])
    print("mammals", len(m), "; fourfold sites per species, median", int(np.median([sum(cnt[s].values()) for s in m.slug])))
    print("lifespan against fourfold composition (mixed, with mass):")
    for b in "ACT":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = smf.mixedlm(f"{b}4 ~ ll + lm", m, groups="order", vc_formula={"family": "0 + C(family)"}).fit()
        print(f"  {b}4 slope {r.params['ll']:+.4f} p {r.pvalues['ll']:.2g}")
    rows = []
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        row = {"aa": aa}
        for label, extra in (("mass", ""), ("mass+composition", " + A4 + C4 + T4")):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                r = smf.mixedlm(f"aa_{aa} ~ ll + lm{extra}", m, groups="order", vc_formula={"family": "0 + C(family)"}).fit()
            row[label + "_slope"], row[label + "_p"] = r.params["ll"], r.pvalues["ll"]
        rows.append(row)
    res = pd.DataFrame(rows).sort_values("mass+composition_p")
    res.to_csv(out / "controls_composition.tsv", sep="\t", index=False, float_format="%.4g")
    print(res.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
