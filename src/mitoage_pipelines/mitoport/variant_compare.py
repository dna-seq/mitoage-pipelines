"""H34, part A: the predictors that already exist, on H29's variants (reported, claimed as nothing).

MitImpact 4 (Castellana et al., Nucleic Acids Research 2026; https://mitimpact.mcb2lab.org) holds, for every
single-base change of the human mitochondrial genome, the scores of published predictors. This module joins five of
them to H29's variants by position, bases and gene, and gives each one's area under the ROC curve for harmful
against harmless beside the area of the site's variability among the build's mammals on the same variants, with
the bootstrap interval of the difference (`variant_stab.compare`). KronGenome's drop (H31) and H29's own score are
read from their tables and treated the same way.

The direction of each score is fixed here from its definition, never from the areas: AlphaMissense and APOGEE 2 are
probabilities of harm; ESM1b is a log-likelihood ratio, lower for a worse change, so its sign is turned; PhyloP is
higher at conserved sites; gnomAD's local constraint (MLC) is higher where variation is missing.

What each predictor knew of labels like ours is said in the table, since an area on labels a predictor was trained
on is not a test of it: APOGEE 2 was trained on disease and population variants of these genes; AlphaMissense was
calibrated with population frequency, which our harmless variants are chosen by; ESM1b, PhyloP and KronGenome saw
no labels; MLC is built from population variation itself.

Writes h34_existing.tsv."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

from . import cybrid, variant_stab

URL = "https://mitimpact.mcb2lab.org/static/downloads/MitImpact4_V059.zip"
ZIP_SHA256 = "f52467ec3a3d5ddff4c7dec8e2a97b120a1b84f1a5bb5589d46e105a88f2edb6"
CSV = "MitImpact4_V059.csv"
CSV_SHA256 = "5045e2a156130861ea40b468c1c53a0743d04a05b0212456a80ac4aa22bc5750"
# name here -> (MitImpact's column, sign that makes higher the prediction of harm, what it knew of labels like ours)
PREDICTORS = {
    "alphamissense": ("AlphaMissense_score", 1, "calibrated with population frequency, which the harmless variants here are chosen by"),
    "esm1b": ("ESM1b_score", -1, "no labels"),
    "apogee2": ("APOGEE2_score", 1, "trained on disease and population variants of these genes: not a test"),
    "phylop_100v": ("PhyloP_100V", 1, "no labels; conservation among 100 vertebrates"),
    "mlc": ("MLC_score", 1, "built from population variation, which the harmless variants here are chosen by"),
}
OURS = {"kron": "no labels (H31)", "stab_fold": "no labels (H29)"}
BASELINE = "base_n_aa"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def download(folder: Path) -> Path:
    """Fetch MitImpact 4's archive, refuse another file than the pinned one, and unpack its table."""
    import io
    import urllib.request
    import zipfile
    folder.mkdir(parents=True, exist_ok=True)
    blob = urllib.request.urlopen(urllib.request.Request(URL, headers={"User-Agent": "mitoage-pipelines"}), timeout=300).read()
    if hashlib.sha256(blob).hexdigest() != ZIP_SHA256:
        raise SystemExit(f"{URL} is not the pinned archive")
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        (folder / CSV).write_bytes(z.read(CSV))
        (folder / "LICENSES.md").write_bytes(z.read("LICENSES.md"))
    read_mitimpact(folder / CSV)
    return folder / CSV


def read_mitimpact(csv: Path, digest: str | None = CSV_SHA256) -> pd.DataFrame:
    """One row per missense change and gene: pos, ref, alt, gene (the 13 names of `cybrid`), and one column per
    predictor with the sign turned so that higher is the prediction of harm. Refuses a file that is not the pinned
    one, and a change given twice for one gene."""
    if digest is not None and sha256(csv) != digest:
        raise SystemExit(f"{csv} is not the pinned MitImpact table: `mitoage download mitimpact`")
    d = pd.read_csv(csv, low_memory=False)
    d = d[d.functional_effect_general == "missense_variant"]
    out = pd.DataFrame({"pos": d.Start.astype(int), "ref": d.Ref, "alt": d.Alt,
                        "gene": [cybrid.gene_name(g, None) for g in d.Gene_symbol]})
    if out.gene.isna().any():
        raise ValueError(f"a gene MitImpact names is not one of the 13: {sorted(set(d.Gene_symbol[out.gene.isna().to_numpy()]))}")
    for name, (col, sign, _) in PREDICTORS.items():
        out[name] = sign * pd.to_numeric(d[col], errors="coerce").to_numpy()
    if out.duplicated(["pos", "ref", "alt", "gene"]).any():
        raise ValueError("MitImpact gives a change twice for one gene")
    return out.reset_index(drop=True)


def joined(variants: pd.DataFrame, mitimpact: pd.DataFrame, kron: pd.DataFrame) -> pd.DataFrame:
    """H29's variants (one row per variant and gene) with the predictors' columns, KronGenome's drop and the
    baseline; a variant a table does not hold keeps an empty cell there, never a zero."""
    v = variants.merge(mitimpact, on=["pos", "ref", "alt", "gene"], how="left", validate="one_to_one")
    k = kron[["pos", "ref", "alt", "kron"]].drop_duplicates()
    if k.duplicated(["pos", "ref", "alt"]).any():
        raise ValueError("KronGenome's table gives a change twice")
    v = v.merge(k, on=["pos", "ref", "alt"], how="left", validate="many_to_one")
    v["harmful"] = v.label == "harmful"
    v[BASELINE] = -v.n_aa
    return v


def table(v: pd.DataFrame) -> pd.DataFrame:
    """One row per score: how many variants it covers, its area, the baseline's on the same variants, the difference
    with its 95% bootstrap interval, and what the score knew of labels like ours."""
    rows = []
    notes = {**{k: n for k, (_, _, n) in PREDICTORS.items()}, **OURS}
    for s, note in notes.items():
        d = v[v[s].notna() & v[BASELINE].notna()]
        h = d.harmful.to_numpy()
        row = {"score": s, "n_harmful": int(h.sum()), "n_harmless": int((~h).sum()), "not_covered": int(len(v) - len(d)),
               "area": np.nan, "baseline_area": np.nan, "difference": np.nan, "low": np.nan, "high": np.nan, "labels": note}
        if h.any() and not h.all():
            a, b, lo, hi = variant_stab.compare(h, d[s].to_numpy(float), d[BASELINE].to_numpy(float))
            row.update(area=a, baseline_area=b, difference=a - b, low=lo, high=hi)
        rows.append(row)
    return pd.DataFrame(rows)


def main(variant_scores: Path, kron_variants: Path, mitimpact: Path, out: Path) -> pd.DataFrame:
    for f in (variant_scores, kron_variants):
        if not f.is_file():
            raise SystemExit(f"{f} is missing")
    v = pd.read_csv(variant_scores, sep="\t")
    v = joined(v, read_mitimpact(mitimpact), pd.read_csv(kron_variants, sep="\t"))
    t = table(v)
    out.mkdir(parents=True, exist_ok=True)
    t.to_csv(out / "h34_existing.tsv", sep="\t", index=False, float_format="%.4f")
    print(t.to_string(index=False)); print("wrote", out / "h34_existing.tsv")
    return t
