"""H3's two-step cross-check of Coevol: per-branch dS and dN/dS by maximum likelihood, then a phylogenetic
regression of the terminal branches' rates on lifespan and body mass.

Step 1 (`fit`): HyPhy fits MG94xREV with a synonymous rate alpha and a non-synonymous rate beta of its own on
every branch (hyphy/local_mg94.bf), on Coevol's inputs: the same species, alignment and topology
(`coevol-prep`). HyPhy needs the tree unrooted, so the two branches at the root become one.

Step 2 (`table`): a branch's length under that model is alpha * S + beta * N, with S and N two constants of the
fitted model (the synonymous and non-synonymous substitutions per site of a branch whose alpha and beta are
1). They are recovered exactly by least squares over the branches, and

    dS = alpha * (S + N)     substitutions per synonymous site     (HyPhy's codon branch lengths are per
    dN = beta  * (S + N)     substitutions per non-synonymous site  nucleotide site)
    dN/dS = beta / alpha

Each branch is matched to the dated tree by the species below it, which gives its duration in millions of
years, and dS per Myr.

Step 3 (`pgls`): on the terminal branches, whose species' lifespan and mass are observed, log10 dN/dS and
log10 dS per Myr are regressed on log10 lifespan and log10 mass, alone and together, with the dated tree's
shared-path covariance scaled by Pagel's lambda (estimated per response on the model with both). The reading
was fixed in docs/HYPOTHESES.md (H3) before the run.
"""
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy import optimize, stats
from statsmodels.stats.multitest import multipletests

from .. import external
from .pgls import gls, nll, shared_path_covariance

HYPHY_BF = Path(__file__).parent / "hyphy" / "local_mg94.bf"
SATURATED_DS = 2.0          # terminal branches above this dS are left out in the sensitivity rows
RESPONSES = {"log10_dNdS": "dNdS", "log10_dS_per_myr": "dS_per_myr"}
MODELS = {"lifespan": ["ll"], "mass": ["lm"], "lifespan+mass": ["ll", "lm"]}


def write_inputs(prep: Path, out: Path) -> tuple[Path, Path]:
    """mammals.fna (Coevol's codons, its missing codons ??? as NNN) and mammals.nwk (its topology)."""
    out.mkdir(parents=True, exist_ok=True)
    lines = (prep / "mammals.ali").read_text().split("\n")
    n, length = map(int, lines[0].split())
    rows = [l.split() for l in lines[1:] if l.strip()]
    if len(rows) != n or any(len(s) != length for _, s in rows):
        raise ValueError(f"{prep / 'mammals.ali'} does not hold {n} sequences of {length} sites")
    fna, nwk = out / "mammals.fna", out / "mammals.nwk"
    fna.write_text("".join(f">{name}\n{seq.replace('?', 'N')}\n" for name, seq in rows))
    nwk.write_text((prep / "mammals.tree").read_text())
    return fna, nwk


def fit(prep: Path, out: Path, cpus: int) -> Path:
    """Run HyPhy on Coevol's inputs; writes out/fit.json. Hours on 16 cores for 106 species."""
    hyphy = external.find("hyphy")
    fna, nwk = write_inputs(prep, out)
    external.run([hyphy, f"CPU={cpus}", HYPHY_BF, "--alignment", fna.name, "--tree", nwk.name,
                  "--code", "Vertebrate-mtDNA", "--output", "fit.json"], cwd=out)
    return out / "fit.json"


def _tips(clade) -> frozenset:
    return frozenset(t.name for t in clade.get_terminals())


def durations(dated_tree: Path, species: list[str], fitted_sets: set[frozenset]) -> dict[frozenset, float]:
    """Millions of years of each branch of the dated tree pruned to `species` (named as in Coevol's files, with
    _ for the dated tree's -), keyed by the species below it. A branch at the root whose sibling is not one of
    `fitted_sets` (HyPhy merged the two root branches into one) gets the sum of the two."""
    tree = Phylo.read(dated_tree, "newick")
    for t in tree.get_terminals():
        t.name = t.name.replace("-", "_")
    keep = set(species)
    missing = keep - {t.name for t in tree.get_terminals()}
    if missing:
        raise ValueError(f"{len(missing)} species are not in the dated tree, e.g. {sorted(missing)[:3]}")
    for t in tree.get_terminals():
        if t.name not in keep:
            tree.prune(t)
    out = {}
    root_children = tree.root.clades
    for cl in tree.find_clades():
        if cl is tree.root:
            continue
        out[_tips(cl)] = cl.branch_length or 0.0
    if len(root_children) == 2:
        a, b = (_tips(c) for c in root_children)
        merged = (root_children[0].branch_length or 0.0) + (root_children[1].branch_length or 0.0)
        for s, other in ((a, b), (b, a)):
            if s in fitted_sets and other not in fitted_sets:
                out[s] = merged
    return out


def table(fit_json: Path, dated_tree: Path) -> tuple[pd.DataFrame, dict]:
    """One row per branch of the fitted tree: alpha, beta, length, dS, dN, dN/dS, duration and dS per Myr."""
    d = json.loads(fit_json.read_text())
    tree = Phylo.read(io.StringIO(d["input"]["trees"]["0"]), "newick")
    bl = d["branch lengths"]
    rows = []
    for cl in tree.find_clades():
        if cl is tree.root:
            continue
        b = bl[cl.name]
        rows.append(dict(branch=cl.name, terminal=cl.is_terminal(), tips=_tips(cl), length=b["MLE"],
                         alpha=b["synonymous rate"]["MLE"], beta=b["non-synonymous rate"]["MLE"]))
    df = pd.DataFrame(rows)
    X = df[["alpha", "beta"]].to_numpy()
    (S, N), *_ = np.linalg.lstsq(X, df.length.to_numpy(), rcond=None)
    worst = float(np.abs(X @ [S, N] - df.length).max())
    if worst > 1e-6 * max(1.0, df.length.max()):
        raise ValueError(f"branch lengths are not alpha*S + beta*N (worst residual {worst:.3g}): not this model's fit")
    df["dS"] = df.alpha * (S + N)
    df["dN"] = df.beta * (S + N)
    df["dNdS"] = np.where(df.alpha > 0, df.beta / df.alpha.where(df.alpha > 0, np.nan), np.nan)
    species = [t.name for t in tree.get_terminals()]
    dur = durations(dated_tree, species, set(df.tips))
    df["myr"] = [dur.get(s, np.nan) for s in df.tips]
    if df.myr.isna().any():
        raise ValueError(f"{int(df.myr.isna().sum())} fitted branches have no branch in the dated tree")
    df["dS_per_myr"] = np.where(df.myr > 0, df.dS / df.myr.where(df.myr > 0, np.nan), np.nan)
    df["species"] = [next(iter(s)) if t else "" for s, t in zip(df.tips, df.terminal)]
    df["n_tips"] = df.tips.map(len)
    info = {"S": float(S), "N": float(N), "synonymous_fraction": float(S / (S + N)),
            "log_likelihood": d["log likelihood"], "parameters": d["parameters"], "branches": len(df)}
    return df.drop(columns="tips"), info


def read_traits(lht: Path) -> pd.DataFrame:
    """Coevol's trait file: species, lifespan (years) and mass (g), as log10 (ll, lm)."""
    lines = [l for l in lht.read_text().split("\n")[2:] if l.strip()]
    df = pd.DataFrame([l.split("\t") for l in lines], columns=["species", "longevity", "mass"])
    df["ll"] = np.log10(df.longevity.astype(float))
    df["lm"] = np.log10(df.mass.astype(float))
    return df[["species", "ll", "lm"]]


def pgls(branches: pd.DataFrame, traits: pd.DataFrame, dated_tree: Path) -> pd.DataFrame:
    """The regressions of step 3 on the terminal branches: one row per response, model, branch set and term."""
    tips = branches[branches.terminal].merge(traits, on="species", how="left", validate="one_to_one")
    if tips[["ll", "lm"]].isna().any().any():
        raise ValueError("terminal branches without lifespan or mass")
    tree = Phylo.read(dated_tree, "newick")
    for t in tree.get_terminals():
        t.name = t.name.replace("-", "_")
    keep = set(tips.species)
    for t in tree.get_terminals():
        if t.name not in keep:
            tree.prune(t)
    order = {t.name: i for i, t in enumerate(tree.get_terminals())}
    C_all = shared_path_covariance(tree.root, order)
    rows = []
    for resp, col in RESPONSES.items():
        for subset in ("all", f"dS<={SATURATED_DS:g}"):
            d = tips.copy()
            d[resp] = np.log10(d[col].where(d[col] > 0))
            dropped = int(d[resp].isna().sum())
            d = d.dropna(subset=[resp])
            if subset != "all":
                d = d[d.dS <= SATURATED_DS]
            idx = np.array([order[s] for s in d.species])
            C = C_all[np.ix_(idx, idx)]
            y = d[resp].to_numpy()
            full = np.column_stack([np.ones(len(d)), d[["ll", "lm"]].to_numpy()])
            lam = float(optimize.minimize_scalar(nll, bounds=(0.0, 1.0), args=(full, y, C), method="bounded").x)
            for model, terms in MODELS.items():
                X = np.column_stack([np.ones(len(d)), d[terms].to_numpy()])
                b, se, df = gls(X, y, C, lam)
                for k, term in enumerate(terms, start=1):
                    t = b[k] / se[k]
                    rows.append(dict(response=resp, branches=subset, model=model, term=term, slope=b[k], se=se[k],
                                     t=t, p=2 * stats.t.sf(abs(t), df), n=len(d), no_log=dropped, lam=lam))
    res = pd.DataFrame(rows)
    q = (res.branches == "all") & (res.model == "lifespan+mass") & (res.term == "ll")
    res["p_holm"] = np.nan
    res.loc[q, "p_holm"] = multipletests(res.loc[q, "p"], method="holm")[1]
    return res


def reading(res: pd.DataFrame) -> list[str]:
    """The reading fixed in H3, applied to the table: the control, then the question."""
    def row(resp, model, term, branches="all"):
        return res[(res.response == resp) & (res.model == model) & (res.term == term) & (res.branches == branches)].iloc[0]
    out = []
    c = row("log10_dNdS", "mass", "lm")
    control = c.slope > 0 and c.p < 0.05
    out.append(f"positive control, dN/dS on mass: slope {c.slope:+.3f} (p {c.p:.2g}) -> {'passes' if control else 'FAILS: nothing below is read'}")
    s = row("log10_dS_per_myr", "mass", "lm")
    out.append(f"secondary control, dS/Myr on mass: slope {s.slope:+.3f} (p {s.p:.2g}), expected negative")
    for resp, sign in (("log10_dNdS", 1), ("log10_dS_per_myr", -1)):
        r = row(resp, "lifespan+mass", "ll")
        if not control:
            verdict = "not read (control failed)"
        elif r.p_holm < 0.05:
            verdict = "agrees with Coevol's sign" if np.sign(r.slope) == sign else "DISAGREES with Coevol's sign"
        else:
            verdict = "no lifespan signal beyond body mass"
        out.append(f"{resp} on lifespan given mass: slope {r.slope:+.3f} (p {r.p:.2g}, Holm {r.p_holm:.2g}) -> {verdict}")
    return out


def main_rates(prep: Path, dated_tree: Path, out: Path, cpus: int, skip_fit: bool = False) -> pd.DataFrame:
    fit_json = out / "fit.json" if skip_fit else fit(prep, out, cpus)
    df, info = table(fit_json, dated_tree)
    df.to_csv(out / "branches.tsv", sep="\t", index=False, float_format="%.6g")
    (out / "model.json").write_text(json.dumps(info, indent=1, sort_keys=True) + "\n")
    print(f"{len(df)} branches ({int(df.terminal.sum())} terminal); synonymous fraction {info['synonymous_fraction']:.3f}; "
          f"lnL {info['log_likelihood']:.2f}; alpha 0 on {int((df.alpha == 0).sum())}, beta 0 on {int((df.beta == 0).sum())}")
    return df


def main_pgls(prep: Path, dated_tree: Path, out: Path) -> pd.DataFrame:
    branches = pd.read_csv(out / "branches.tsv", sep="\t", keep_default_na=False, na_values=[""])
    res = pgls(branches, read_traits(prep / "mammals.lht"), dated_tree)
    res.to_csv(out / "pgls.tsv", sep="\t", index=False, float_format="%.4g")
    with pd.option_context("display.width", 200):
        print(res.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    for line in reading(res):
        print(line)
    return res
