"""H20: temperature and metabolic rate against the mitochondrial substitution rate and composition.

Metabolic theory predicts that substitution rates rise with temperature and fall with body size, through their
effect on metabolism (Gillooly, Allen, West & Brown 2005, PNAS 102:140, doi:10.1073/pnas.0407735101):

    rate  ~  M^(-1/4) * exp(-E / kT)      so    log10 rate = a - (E / ln 10) * (1/kT) + b * log10 size

with k Boltzmann's constant in eV/K and T the body temperature in kelvin. A slope on 1/kT is read as
E = -ln(10) * slope, in eV: positive when warmer animals evolve faster. The body-size side is the published
positive control (mtDNA rates fall with body mass across mammals: Nabholz, Glemin & Galtier 2008, MBE 25:120).

Temperatures. Mammals: AnAge's body temperature (K), kept between 270 and 320 K. Ray-finned and cartilaginous
fishes, whose body is at the water's temperature: FishBase's preferred temperature, the median over a species'
stocks of `stocks.TempPref50` (degrees C, kept between -3 and 40), and as a sensitivity the median over stocks of
the midpoint of `TempMin` and `TempMax`. Birds have almost no body temperatures in AnAge (3 species) and are read
for metabolic rate only.
Size: log10 adult mass (AnAge) for mammals and birds; for fishes, whose mass is mostly missing, log10 of FishBase's
`species.Length` (maximum length, cm).
Metabolic rate: log10 AnAge metabolic rate (W), always with log10 mass in the model, so that its slope is that of
the mass-specific rate.

FishBase names. A build species is matched to a FishBase SpecCode by its scientific name in `species` (Genus
Species), or else by `synonyms` rows of status accepted, provisionally accepted or synonym. A name that points to
two codes is ambiguous and stays unmatched. Where the build's own lifespan came from a FishBase `species` record,
that record names a SpecCode; a name match that disagrees with it stays unmatched too. Every outcome is counted.

1. Rates (`pairs.tsv`, `pair_fits.tsv`). For each class and model, the species that have every trait of the
   model, a tip in the dated tree and all 13 protein alignments, and as many pairs of them as the tree allows with
   no branch in two pairs (pairwise comparisons, Maddison 2000, J Theor Biol 202:195; `disjoint_pairs`). Per pair: dS and dN between the two sisters over the concatenated genes (Nei-Gojobori with the
   Jukes-Cantor correction, in their mitochondrial code, as H7 counts them), divided by twice the age of their
   common ancestor to give substitutions per site per million years along each branch; each trait is the mean of
   the two sisters'. Pairs with dS above 1 (saturated) or with dS or dN of 0 are left out and counted. log10 dS
   per Myr, log10 dN per Myr and log10 dN/dS are regressed on the traits by GLS over the pairs, each pair
   standing at its first species' tip of the dated tree, with Pagel's lambda estimated by maximum likelihood.
   Beside it, every fit is repeated on the plain cherries of the whole tree whose two sisters both have the
   model's traits (`pairing` "cherries"): fewer pairs, younger ones.
2. Composition (`composition.tsv`). H11's 43 fixed features (GC content and skews, fourfold bases, the C share at
   third positions, the 32 pyrimidine-ending codons), each standardised and regressed on the traits per species on
   the dated tree, Pagel's lambda estimated. Exploratory: counted at q < 0.05 (Benjamini-Hochberg over the
   features within class, model and term).
3. Lehmann et al. retested (`lehmann.tsv`). Lehmann, Segal, Muradian & Fraifeld 2008 (Rejuvenation Res 11:409,
   doi:10.1089/rej.2008.0676) found that mtDNA GC content and resting metabolic rate together explain 77% of
   the variance of mammalian maximum lifespan over 140 species, in an analysis without the tree; Lehmann,
   Muradian & Fraifeld 2013 (Front Genet 4:111) took body temperature in as well. log10 lifespan on the GC
   content of the whole record and log10 metabolic rate, then with log10 mass, then with 1/kT: ordinary least
   squares (R^2, as published) beside the phylogenetic regression.

The pre-registered tests and their reading are in docs/HYPOTHESES.md (H20); `readout.tsv` applies them.
"""
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Data import CodonTable
from scipy import optimize, stats
from statsmodels.stats.multitest import multipletests

from ..genes import MT_GENES
from ..seqio import read_fasta
from .aa_pgls import class_covariance, default_dataset
from .generation_time import features
from .pgls import gls, nll
from .rates import SATURATED_DS, cherries, distance

K_EV = 8.617333262e-5                 # Boltzmann's constant, eV/K
BODY_K = (270.0, 320.0)               # AnAge body temperatures kept, K
WATER_C = (-3.0, 40.0)                # FishBase temperatures kept, degrees C
SYNONYM_STATUS = ("accepted name", "provisionally accepted name", "synonym")
FISHES = ("Actinopterygii", "Chondrichthyes")
CLASSES = ("Mammalia", "Aves", *FISHES)
MIN_PAIRS = 20
MIN_SPECIES = 25
GC = "gc_percent:total_mtDNA"
PAIRINGS = ("disjoint", "cherries")   # the primary pairing, and the plain cherries of the whole tree (reported beside)
RESPONSES = {"log10_dS_per_myr": "dS_per_myr", "log10_dN_per_myr": "dN_per_myr", "log10_dNdS": "dNdS"}
# (model, terms) per class; "size" is log10 mass, or log10 length in fishes
PAIR_MODELS = {
    "Mammalia": [("size", ["size"]), ("temperature", ["inv_kt"]), ("temperature+size", ["inv_kt", "size"]),
                 ("metabolic+size", ["lmr", "size"]), ("temperature+metabolic+size", ["inv_kt", "lmr", "size"])],
    "Aves": [("size", ["size"]), ("metabolic+size", ["lmr", "size"])],
    "Actinopterygii": [("size", ["size"]), ("temperature", ["inv_kt"]), ("temperature+size", ["inv_kt", "size"]),
                       ("temperature_range+size", ["inv_kt_range", "size"])],
    "Chondrichthyes": [("size", ["size"]), ("temperature", ["inv_kt"]), ("temperature+size", ["inv_kt", "size"])],
}
COMPOSITION_MODELS = {
    "Mammalia": [("temperature", ["inv_kt"]), ("temperature+size", ["inv_kt", "size"]),
                 ("metabolic+size", ["lmr", "size"]), ("lifespan+size", ["ll", "size"]),
                 ("lifespan+temperature+metabolic+size", ["ll", "inv_kt", "lmr", "size"])],
    "Aves": [("metabolic+size", ["lmr", "size"])],
    "Actinopterygii": [("temperature", ["inv_kt"]), ("temperature+size", ["inv_kt", "size"])],
    "Chondrichthyes": [("temperature", ["inv_kt"]), ("temperature+size", ["inv_kt", "size"])],
}
# The pre-registered tests (docs/HYPOTHESES.md, H20): (id, class, model, term, predicted sign of the slope)
PRIMARY = [("T1", "Actinopterygii", "temperature+size", "inv_kt", -1),
           ("T2", "Chondrichthyes", "temperature+size", "inv_kt", -1),
           ("T3", "Mammalia", "temperature+size", "inv_kt", -1),
           ("M1", "Mammalia", "metabolic+size", "lmr", +1),
           ("M2", "Aves", "metabolic+size", "lmr", +1)]
CONTROLS = [("P1", "Mammalia", "size", "size", -1), ("P2", "Actinopterygii", "size", "size", -1)]
LEHMANN_MODELS = [("gc", [GC]), ("gc+metabolic", [GC, "lmr"]), ("gc+metabolic+mass", [GC, "lmr", "size"]),
                  ("gc+metabolic+mass+temperature", [GC, "lmr", "size", "inv_kt"])]
BUILD_FB_ID = re.compile(r"^species\.\w+:(\d+)$")


# --- traits ---------------------------------------------------------------------------------------
def inv_kt(kelvin: pd.Series) -> pd.Series:
    return 1.0 / (K_EV * kelvin)


def fishbase_temperatures(fb: Path) -> pd.DataFrame:
    """Per FishBase SpecCode: tpref (median over stocks of TempPref50, C), trange (median over stocks of the
    midpoint of TempMin and TempMax, C), length_cm (species.Length). Values outside WATER_C are not used."""
    st = pd.read_parquet(fb / "stocks.parquet", columns=["SpecCode", "TempMin", "TempMax", "TempPref50"])
    ok = lambda s: s.where(s.between(*WATER_C))
    lo, hi = ok(st.TempMin), ok(st.TempMax)
    st = st.assign(tpref=ok(st.TempPref50), trange=((lo + hi) / 2).where(lo <= hi))
    temps = st.groupby("SpecCode")[["tpref", "trange"]].median()
    sp = pd.read_parquet(fb / "species.parquet", columns=["SpecCode", "Length"]).set_index("SpecCode")
    out = sp.rename(columns={"Length": "length_cm"}).join(temps, how="outer")
    out["length_cm"] = out.length_cm.where(out.length_cm > 0)
    return out


def fishbase_codes(fb: Path, species: pd.DataFrame) -> pd.DataFrame:
    """slug -> SpecCode (nullable) and how it was or was not found (see the module's docstring)."""
    sp = pd.read_parquet(fb / "species.parquet", columns=["SpecCode", "Genus", "Species"])
    syn = pd.read_parquet(fb / "synonyms.parquet", columns=["SpecCode", "SynGenus", "SynSpecies", "Status"])
    syn = syn[syn.Status.isin(SYNONYM_STATUS)]
    by_name = lambda keys, codes: pd.DataFrame({"k": keys.str.strip(), "c": codes}).dropna().groupby("k").c.agg(lambda c: sorted(set(c)))
    exact = by_name(sp.Genus + " " + sp.Species, sp.SpecCode)
    synonym = by_name(syn.SynGenus + " " + syn.SynSpecies, syn.SpecCode)
    rows = []
    for slug, name, src, rec in zip(species.slug, species.scientific_name, species.lifespan_source, species.lifespan_source_record_id):
        name = str(name).strip()
        codes, how = exact.get(name), "exact name"
        if codes is None:
            codes, how = synonym.get(name), "synonym"
        if codes is None:
            rows.append((slug, pd.NA, "no match"))
            continue
        if len(codes) > 1:
            rows.append((slug, pd.NA, f"ambiguous {how}"))
            continue
        code = int(codes[0])
        m = BUILD_FB_ID.match(str(rec)) if src == "FishBase" else None
        if m and int(m.group(1)) != code:
            rows.append((slug, pd.NA, "conflicts with the build's FishBase record"))
            continue
        rows.append((slug, code, how))
    return pd.DataFrame(rows, columns=["slug", "SpecCode", "how"]).astype({"SpecCode": "Int64"})


def traits(species: pd.DataFrame, fb: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One row per species: class, ll, size, lmr, temp_k, inv_kt, inv_kt_range, and where each came from.
    -> (traits, the FishBase matching of the fishes)."""
    log = lambda s: np.log10(s.where(s > 0))
    fish = species["class"].isin(FISHES)
    t = pd.DataFrame({"slug": species.slug, "class": species["class"], "ll": log(species.max_longevity_years),
                      "lm": log(species.adult_weight_g), "lmr": log(species.metabolic_rate_w)})
    body = species.body_temperature_k.where(species.body_temperature_k.between(*BODY_K))
    t["temp_k"] = body.where(~fish)
    t["temp_source"] = np.where(~fish & body.notna(), "AnAge body temperature", "")
    match = fishbase_codes(fb, species[fish])
    fbt = fishbase_temperatures(fb)
    m = match.merge(fbt, left_on="SpecCode", right_index=True, how="left").set_index("slug")
    f = t.slug.where(fish)
    t["fishbase_speccode"] = f.map(m.SpecCode).astype("Int64")
    tpref, trange = f.map(m.tpref), f.map(m.trange)
    t.loc[fish, "temp_k"] = tpref[fish] + 273.15
    t.loc[fish & tpref.notna(), "temp_source"] = "FishBase TempPref50"
    t["inv_kt"] = inv_kt(t.temp_k)
    t["inv_kt_range"] = inv_kt(trange + 273.15)
    t["llen"] = log(f.map(m.length_cm))
    t["size"] = t.lm.where(~fish, t.llen)
    return t, match


# --- 1. rates over sister pairs ----------------------------------------------------------------------
def disjoint_pairs(tree, members: set[str]) -> list[tuple[str, str, float]]:
    """Pairs of `members` whose paths on the dated tree share no branch (pairwise comparisons, Maddison 2000),
    as many as the tree allows: walking up from the tips, the members left unpaired below a node, one per child
    clade at most, are paired there two by two (in order of name; a polytomy can leave one over, which goes up).
    Each pair comes with the age of the node that joins it (the mean distance from it to the two tips)."""
    out = []
    def walk(clade) -> list:
        if clade.is_terminal():
            return [(clade.name, 0.0)] if clade.name in members else []
        up = []
        for ch in clade.clades:
            up += [(name, d + (ch.branch_length or 0.0)) for name, d in walk(ch)]
        up.sort()
        while len(up) >= 2:
            (x, dx), (y, dy) = up.pop(0), up.pop(0)
            out.append((x, y, (dx + dy) / 2))
        return up
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 10_000))
    walk(tree.root)
    return sorted(out)


class Distances:
    """dS and dN between two species over the 13 concatenated genes, computed once per pair."""
    def __init__(self, data: Path, species: pd.DataFrame):
        self.aln = {g: read_fasta(data / "alignments" / f"{g}.codon.fasta.gz") for g in MT_GENES}
        self.aligned = set.intersection(*[set(a) for a in self.aln.values()])
        self.table = species.set_index("slug").transl_table
        self.cache: dict[tuple[str, str], dict] = {}

    def code(self, a: str, b: str) -> int:
        tables = {int(self.table[s]) for s in (a, b) if not pd.isna(self.table.get(s))}
        t = tables.pop() if len(tables) == 1 else 2
        return t if t in CodonTable.unambiguous_dna_by_id else 2

    def __call__(self, a: str, b: str) -> dict:
        if (a, b) not in self.cache:
            sa, sb = ("".join(self.aln[g][s] for g in MT_GENES) for s in (a, b))
            self.cache[(a, b)] = distance(sa, sb, self.code(a, b))
        return self.cache[(a, b)]


def pair_table(pairs: list[tuple[str, str, float]], dist: Distances, tr: pd.DataFrame, terms: list[str]) -> pd.DataFrame:
    """One row per pair: dS, dN, the rates per million years along each branch, and the mean of the two sisters'
    value of each trait in `terms`."""
    t = tr.set_index("slug")
    rows = []
    for a, b, age in pairs:
        d = dist(a, b)
        rows.append({"sp_a": a, "sp_b": b, "age_my": age, "codons": d["codons"], "dS": d["dS"], "dN": d["dN"],
                     **{c: (t.at[a, c] + t.at[b, c]) / 2 for c in terms}})
    p = pd.DataFrame(rows, columns=["sp_a", "sp_b", "age_my", "codons", "dS", "dN", *terms])
    per_branch = 2 * p.age_my.where(p.age_my > 0)
    p["dS_per_myr"] = p.dS / per_branch
    p["dN_per_myr"] = p.dN / per_branch
    p["dNdS"] = p.dN / p.dS.where(p.dS > 0)
    p["saturated"] = p.dS > SATURATED_DS
    p["usable"] = ~p.saturated & (p.dS > 0) & (p.dN > 0) & (p.age_my > 0)
    return p


def fit(C: np.ndarray, y: np.ndarray, X: np.ndarray, names: list[str], lam: float | None = None) -> dict:
    """GLS of y on an intercept and X under Pagel's lambda (None: maximum likelihood in [0, 1])."""
    X1 = np.column_stack([np.ones(len(y)), X])
    if lam is None:
        lam = float(optimize.minimize_scalar(nll, bounds=(0.0, 1.0), args=(X1, y, C), method="bounded").x)
    b, se, df = gls(X1, y, C, lam)
    out = {"lambda": lam, "df": df}
    for j, n in enumerate(names, start=1):
        t = b[j] / se[j] if se[j] > 0 else np.nan
        out[n] = (b[j], se[j], t, 2 * stats.t.sf(abs(t), df) if np.isfinite(t) else np.nan)
    return out


def fit_pairs(species: pd.DataFrame, tr: pd.DataFrame, dist: Distances, tree_file: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per class and model: the pairs of the tree pruned to the class's species that have every term, a tip and
    the 13 alignments; GLS of each response on the terms. -> (fits, the pairs of every model)."""
    tree = Phylo.read(tree_file, "newick")
    tips = {t.name for t in tree.get_terminals()}
    rows, tables = [], []
    for cls, models in PAIR_MODELS.items():
        cl = tr[tr["class"].eq(cls) & tr.slug.isin(tips & dist.aligned)]
        for model, terms, pairing in [(m, t, k) for m, t in models for k in PAIRINGS]:
            members = set(cl.dropna(subset=terms).slug)
            found = (disjoint_pairs(tree, members) if pairing == "disjoint" else
                     [(a, b, float(np.mean([k.branch_length or 0.0 for k in node.clades]))) for node, a, b in cherries(tree, members)])
            p = pair_table(found, dist, tr, terms)
            tables.append(p.assign(**{"class": cls, "model": model, "pairing": pairing}))
            d = p[p.usable].reset_index(drop=True)
            base = {"class": cls, "model": model, "pairing": pairing, "n_species": len(members), "n_pairs_all": len(p),
                    "n_saturated": int(p.saturated.sum()), "n_pairs": len(d)}
            if len(d) < MIN_PAIRS:
                rows += [{**base, "response": r, "term": t, "status": f"not assessed: {len(d)} pairs, fewer than {MIN_PAIRS}"}
                         for r in RESPONSES for t in terms]
                continue
            C, order, _ = class_covariance(tree_file, list(d.sp_a))
            d = d.set_index("sp_a").loc[order].reset_index()
            for resp, col in RESPONSES.items():
                y = np.log10(d[col].to_numpy(dtype=float))
                r = fit(C, y, d[terms].to_numpy(dtype=float), terms)
                for t in terms:
                    b, se, tt, pv = r[t]
                    row = {**base, "response": resp, "term": t, "status": "assessed", "lambda": r["lambda"],
                           "slope": b, "se": se, "t": tt, "df": r["df"], "p": pv}
                    if t.startswith("inv_kt"):
                        q = stats.t.ppf(0.975, r["df"])
                        row |= {"E_eV": -math.log(10) * b, "E_low": -math.log(10) * (b + q * se),
                                "E_high": -math.log(10) * (b - q * se)}
                    rows.append(row)
    res = pd.DataFrame(rows)
    ok = res.status.eq("assessed")
    res.loc[ok, "q"] = res[ok].groupby(["class", "pairing", "response"]).p.transform(
        lambda p: pd.Series(multipletests(p, method="fdr_bh")[1], index=p.index))
    return res, pd.concat(tables, ignore_index=True)


# --- 2. composition ---------------------------------------------------------------------------------
def fit_composition(tr: pd.DataFrame, feats: pd.DataFrame, tree_file: Path) -> pd.DataFrame:
    rows = []
    present = {t.name for t in Phylo.read(tree_file, "newick").get_terminals()}
    names = [c for c in feats.columns if c != "slug"]
    for cls, models in COMPOSITION_MODELS.items():
        sp = tr[tr["class"].eq(cls) & tr.slug.isin(present)].merge(feats, on="slug")
        for model, terms in models:
            d = sp.dropna(subset=terms).reset_index(drop=True)
            if len(d) < MIN_SPECIES:
                rows.append({"class": cls, "model": model, "n": len(d), "status": f"not assessed: fewer than {MIN_SPECIES} species"})
                continue
            C, order, _ = class_covariance(tree_file, list(d.slug))
            d = d.set_index("slug").loc[order].reset_index()
            X = d[terms].to_numpy(dtype=float)
            for f in names:
                x = d[f].to_numpy(dtype=float)
                ok = np.isfinite(x)
                if ok.sum() < MIN_SPECIES or np.std(x[ok]) == 0:
                    continue
                y = (x[ok] - x[ok].mean()) / x[ok].std()
                r = fit(C[np.ix_(ok, ok)], y, X[ok], terms)
                for t in terms:
                    b, se, tt, p = r[t]
                    rows.append({"class": cls, "model": model, "term": t, "feature": f, "n": int(ok.sum()),
                                 "status": "assessed", "lambda": r["lambda"], "slope_sd": b, "se": se, "t": tt, "p": p})
    res = pd.DataFrame(rows)
    ok = res.status.eq("assessed")
    res.loc[ok, "q"] = res[ok].groupby(["class", "model", "term"]).p.transform(
        lambda p: pd.Series(multipletests(p, method="fdr_bh")[1], index=p.index))
    return res


# --- 3. Lehmann et al. -----------------------------------------------------------------------------
def lehmann(species: pd.DataFrame, tr: pd.DataFrame, feats: pd.DataFrame, tree_file: Path) -> pd.DataFrame:
    present = {t.name for t in Phylo.read(tree_file, "newick").get_terminals()}
    mam = species[species["class"].eq("Mammalia") & default_dataset(species)].slug
    sp = tr[tr.slug.isin(set(mam) & present)].merge(feats[["slug", GC]], on="slug")
    rows = []
    for model, terms in LEHMANN_MODELS:
        d = sp.dropna(subset=["ll", *terms]).reset_index(drop=True)
        if len(d) < MIN_SPECIES:
            rows.append({"model": model, "n": len(d), "status": f"not assessed: fewer than {MIN_SPECIES} species"})
            continue
        C, order, _ = class_covariance(tree_file, list(d.slug))
        d = d.set_index("slug").loc[order].reset_index()
        X, y = d[terms].to_numpy(dtype=float), d.ll.to_numpy(dtype=float)
        X1 = np.column_stack([np.ones(len(y)), X])
        b_ols, *_ = np.linalg.lstsq(X1, y, rcond=None)
        r2 = 1 - ((y - X1 @ b_ols) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        pg = fit(C, y, X, terms)
        for t in terms:
            b, se, tt, p = pg[t]
            rows.append({"model": model, "n": len(d), "status": "assessed", "term": t, "ols_r2": r2,
                         "ols_slope": b_ols[1 + terms.index(t)], "lambda": pg["lambda"], "slope": b, "se": se,
                         "t": tt, "p": p})
    return pd.DataFrame(rows)


def _cherries(r) -> dict:
    """The same test on the plain cherries of the whole tree, reported beside the primary pairing."""
    if r is None or r.status != "assessed":
        return {"cherries_status": "not assessed" if r is None else r.status}
    return {"cherries_status": "assessed", "cherries_n_pairs": r.n_pairs, "cherries_slope": r.slope, "cherries_p": r.p}


# --- the pre-registered readout -----------------------------------------------------------------------
def readout(fits: pd.DataFrame, leh: pd.DataFrame) -> pd.DataFrame:
    """The tests of docs/HYPOTHESES.md (H20) on log10 dS per Myr, with Holm over the primary tests that could be
    assessed, and the words each result is read with."""
    def find(cls, model, term, pairing="disjoint"):
        r = fits[fits["class"].eq(cls) & fits.model.eq(model) & fits.term.eq(term) & fits.response.eq("log10_dS_per_myr")
                 & fits.pairing.eq(pairing)]
        return r.iloc[0] if len(r) else None
    rows = []
    for kind, tests in (("control", CONTROLS), ("primary", PRIMARY)):
        for tid, cls, model, term, sign in tests:
            r = find(cls, model, term)
            assessed = r is not None and r.status == "assessed"
            rows.append({"test": tid, "kind": kind, "class": cls, "model": model, "term": term, "predicted_sign": sign,
                         "status": r.status if r is not None else "not assessed: no fit",
                         **({k: r[k] for k in ("n_pairs", "lambda", "slope", "se", "p")} if assessed else {}),
                         **({k: r[k] for k in ("E_eV", "E_low", "E_high")} if assessed and term.startswith("inv_kt") else {}),
                         **_cherries(find(cls, model, term, "cherries"))})
    out = pd.DataFrame(rows)
    prim = out.kind.eq("primary") & out.status.eq("assessed")
    out.loc[prim, "p_holm"] = multipletests(out.loc[prim, "p"], method="holm")[1]
    out.loc[out.kind.eq("control"), "p_holm"] = out.loc[out.kind.eq("control"), "p"]
    failed = set(out[out.kind.eq("control") & ~(out.status.eq("assessed") & (out.p < 0.05)
                                                & (np.sign(out.slope) == out.predicted_sign))]["class"])

    def words(r):
        if r.status != "assessed":
            return "not assessed"
        if r.kind == "primary" and r["class"] in failed:
            return "not read: the size control of this class did not pass"
        if r.p_holm < 0.05:
            return "as predicted" if np.sign(r.slope) == r.predicted_sign else "opposite to the prediction"
        return "not detected"
    out["reading"] = out.apply(words, axis=1)
    for tid, model in (("L1", "gc+metabolic"), ("L2", "gc+metabolic+mass+temperature")):
        r = leh[leh.model.eq(model) & leh.term.eq(GC)] if "term" in leh.columns else leh.iloc[0:0]
        if len(r):
            r = r.iloc[0]
            rows_l = {"test": tid, "kind": "retest", "class": "Mammalia", "model": model, "term": GC, "predicted_sign": +1,
                      "status": r.status, "n_species": r.n, "lambda": r["lambda"], "slope": r.slope,
                      "se": r.se, "p": r.p, "p_holm": r.p, "ols_r2": r.ols_r2}
            rows_l["reading"] = ("as published" if r.p < 0.05 and r.slope > 0 else
                                 "opposite to the published" if r.p < 0.05 else "not reproduced on the tree")
        else:
            rows_l = {"test": tid, "kind": "retest", "class": "Mammalia", "model": model, "term": GC,
                      "status": "not assessed", "reading": "not assessed"}
        out = pd.concat([out, pd.DataFrame([rows_l])], ignore_index=True)
    return out


def main(data: Path, fourfold: Path, tree_file: Path, fishbase: Path, out: Path) -> pd.DataFrame:
    out.mkdir(parents=True, exist_ok=True)
    species = pd.read_parquet(data / "species.parquet")
    species = species[species.record_status.ne("retired")].reset_index(drop=True)
    tr, match = traits(species, fishbase)
    tr.to_csv(out / "traits.tsv", sep="\t", index=False, float_format="%.6g")
    match.how.value_counts().rename_axis("how").reset_index().to_csv(out / "fishbase_match.tsv", sep="\t", index=False)
    print("FishBase names of the fishes:", dict(match.how.value_counts()))
    for cls in CLASSES:
        c = tr[tr["class"].eq(cls)]
        print(f"{cls}: {len(c)} species; temperature {c.inv_kt.notna().sum()}, size {c['size'].notna().sum()}, "
              f"metabolic rate {c.lmr.notna().sum()}")
    fits, pairs = fit_pairs(species, tr, Distances(data, species), tree_file)
    pairs.to_csv(out / "pairs.tsv", sep="\t", index=False, float_format="%.6g")
    fits.to_csv(out / "pair_fits.tsv", sep="\t", index=False, float_format="%.5g")
    feats = features(data, fourfold)
    fit_composition(tr, feats, tree_file).to_csv(out / "composition.tsv", sep="\t", index=False, float_format="%.5g")
    leh = lehmann(species, tr, feats, tree_file)
    leh.to_csv(out / "lehmann.tsv", sep="\t", index=False, float_format="%.5g")
    ro = readout(fits, leh)
    ro.to_csv(out / "readout.tsv", sep="\t", index=False, float_format="%.5g")
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print(ro.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return ro
