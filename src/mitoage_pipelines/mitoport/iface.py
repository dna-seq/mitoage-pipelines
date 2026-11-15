"""Do the contact residues of complex IV change with their nuclear partner, and does that differ with lifespan?

Two sides, both taken from the solved human enzyme. Mitochondrial side: the contact residues of COX1-3
against the burial-matched non-contact control of `mitoport burial`. Nuclear side: for each nuclear subunit,
the residues it uses to touch a mitochondrial one, against its remaining residues.

Test 1 (coupling): over species pairs, the distance on one side against the distance on the other, with the
patristic distance of the dated tree and its logarithm partialled out (Spearman on residuals), for the four
combinations of the two sets per side. A nuclear protein's mitochondrial-side excess is the correlation of
its interface with the mitochondrial contacts minus the same with the burial-matched control; its
nuclear-side excess is that correlation minus the one with its other columns. The cytosolic control
proteins, which touch nothing mitochondrial, give the null (random column sets of the usual interface size,
seed 0, stand in for their interface). Two distances: the fraction of differing columns (`seq_`), and the
cosine distance of ESM-2 650M per-residue embeddings averaged over the position set (`esm_`). Species pairs
are not independent, so groups of proteins are compared with the controls (Mann-Whitney over proteins).

Test 2 (lifespan): per species, the fraction of contact columns that differ from the column's modal residue
among the species of the set, the same over the matched control columns, and their difference, each
standardised; phylogenetic regression of log10 maximum lifespan on it with log10 adult mass, Pagel's lambda
estimated on the model without it. The slope comes with its standard error and 95% interval.

Writes coupling.tsv and lifespan.tsv. Ported from the cluster script `mitoport/iface/iface.py` of 2026-11-14
without change of method; the interval columns of lifespan.tsv are added here.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1
from scipy import linalg, optimize, stats

from .. import external
from ..genes import COMPLEX_II, COMPLEX_III_NUCLEAR, COMPLEX_IV_NUCLEAR, MATRIX
from ..seqio import read_fasta
from .burial import column_of_position

MTG = ["COX1", "COX2", "COX3"]
ESM_COLUMNS = ["contact_iface", "control_iface", "contact_other", "control_other", "excess_mt", "excess_nuc",
               "own_contacts_iface"]


# --- the nuclear alignments ---------------------------------------------------------------------
def align(sel: Path, out: Path, threads: int = 2, jobs: int = 8) -> None:
    """Every selection <GENE>.faa with MAFFT L-INS-i (--localpair --maxiterate 1000 --anysymbol), unfiltered,
    into <out>/<GENE>.afa, `jobs` genes at a time. With more than one thread MAFFT's refinement is not
    repeatable run to run; the published run used 2 threads and 8 jobs."""
    mafft = external.find("mafft", env="MAFFT")
    out.mkdir(parents=True, exist_ok=True)
    faas = sorted(sel.glob("*.faa"))
    if not faas:
        raise FileNotFoundError(f"no selections under {sel}: run `mitoage mitoport select` first")

    def one(p: Path) -> None:
        with open(out / f"{p.stem}.afa", "w") as f:
            external.run([mafft, "--localpair", "--maxiterate", "1000", "--anysymbol", "--quiet",
                          "--thread", str(threads), p], stdout=f)
    with ThreadPoolExecutor(jobs) as ex:
        list(ex.map(one, faas))
    for p in faas:
        a = read_fasta(out / f"{p.stem}.afa")
        print(f"{p.stem}\t{len(a)} sequences\t{len(next(iter(a.values())))} columns")


# --- distances ----------------------------------------------------------------------------------
def hamd(aln: dict[str, str], sp: list[str], cols) -> np.ndarray:
    """Share of differing residues over the columns both species fill; NaN where fewer than 5 are shared."""
    M = np.array([[aln[s][c] for c in cols] for s in sp])
    ok = M != "-"
    D = np.full((len(sp), len(sp)), np.nan)
    for i in range(len(sp)):
        both = ok[i] & ok
        diff = (M[i] != M) & both
        n = both.sum(axis=1)
        D[i] = np.where(n >= 5, diff.sum(axis=1) / np.maximum(n, 1), np.nan)
    return D


def cosd(A: np.ndarray) -> np.ndarray:
    ok = ~np.isnan(A[:, 0]); B = A.copy(); B[~ok] = 0
    norm = np.linalg.norm(B, axis=1, keepdims=True)
    B = B / np.where(norm == 0, 1, norm)
    D = 1 - B @ B.T; D[~ok, :] = np.nan; D[:, ~ok] = np.nan
    return D


def resid(y: np.ndarray, Z: list[np.ndarray]) -> np.ndarray:
    Z = np.column_stack([np.ones(len(y))] + Z); b, *_ = np.linalg.lstsq(Z, y, rcond=None)
    return y - Z @ b


def pcorr(x: np.ndarray, y: np.ndarray, T: np.ndarray) -> float:
    """Spearman correlation of x and y with T and log(1 + T) partialled out; NaN under 200 usable pairs."""
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 200:
        return np.nan
    cov = [T[ok], np.log1p(T[ok])]
    rx, ry = resid(stats.rankdata(x[ok]), cov), resid(stats.rankdata(y[ok]), cov)
    return float(np.corrcoef(rx, ry)[0, 1])


def mismatch(aln: dict[str, str], sp: list[str], cols) -> tuple[np.ndarray, np.ndarray]:
    """(share of filled columns where the species differs from the column's modal residue among `sp`,
    share of the columns it fills). Gaps are left out of both the mode and the count."""
    M = np.array([[aln[s][c] for c in cols] for s in sp]); ok = M != "-"

    def modal(col):
        v, n = np.unique(col, return_counts=True); return v[n.argmax()]
    mode = np.array([modal(M[ok[:, j], j]) if ok[:, j].any() else "-" for j in range(M.shape[1])])
    diff = (M != mode) & ok
    return diff.sum(axis=1) / np.maximum(ok.sum(axis=1), 1), ok.mean(axis=1)


# --- the phylogenetic regression -----------------------------------------------------------------
def V(lam: float, C: np.ndarray) -> np.ndarray:
    v = lam * C; np.fill_diagonal(v, np.diag(C)); return v


def nll(lam: float, X: np.ndarray, y: np.ndarray, C: np.ndarray) -> float:
    L = linalg.cholesky(V(lam, C), lower=True)
    Xw, yw = linalg.solve_triangular(L, X, lower=True), linalg.solve_triangular(L, y, lower=True)
    b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); r = yw - Xw @ b
    return 0.5 * len(y) * np.log(r @ r / len(y)) + np.log(np.diag(L)).sum()


def gls(X: np.ndarray, y: np.ndarray, C: np.ndarray, lam: float):
    L = linalg.cholesky(V(lam, C), lower=True)
    Xw, yw = linalg.solve_triangular(L, X, lower=True), linalg.solve_triangular(L, y, lower=True)
    b, *_ = np.linalg.lstsq(Xw, yw, rcond=None); r = yw - Xw @ b; df = len(y) - X.shape[1]
    cov = np.linalg.pinv(Xw.T @ Xw) * (r @ r / df)
    return b, np.sqrt(np.diag(cov)), df


# --- the run ------------------------------------------------------------------------------------
def esm_load(folder: Path, prefix: str):
    idx = pd.read_parquet(folder / f"{prefix}-esm2_650m.index.parquet")
    X = np.load(folder / f"{prefix}-esm2_650m.residues.npy", mmap_mode="r")
    return {r.id: (r.offset, r.length, r.truncated) for r in idx.itertuples()}, X


def esm_missing(folder: Path | None) -> list[str]:
    if folder is None:
        return ["no folder given"]
    return [str(folder / f"{p}-esm2_650m.{e}") for p in ("mt13", "nuclear") for e in ("index.parquet", "residues.npy")
            if not (folder / f"{p}-esm2_650m.{e}").is_file()]


def main(work: Path, nucaln: Path, contacts: Path, cif: Path, build: Path, esm: Path | None, tree_file: Path,
         out: Path, skip_esm: bool = False) -> None:
    out.mkdir(parents=True, exist_ok=True)
    missing = esm_missing(esm)
    if missing and not skip_esm:
        raise FileNotFoundError(
            "the ESM-2 per-residue embeddings are missing (" + ", ".join(missing) + "). They are not published "
            "(about 7 GB); give --esm, or --skip-esm to leave every esm_ column empty (not assessed).")
    use_esm = not missing
    if not use_esm:
        print("NOTE: no per-residue embeddings: every esm_ column of coupling.tsv is left empty, NOT ASSESSED")
    rng = np.random.default_rng(0)
    CIV, CIII, CII = COMPLEX_IV_NUCLEAR, COMPLEX_III_NUCLEAR, COMPLEX_II

    sites = pd.read_csv(work / "sites.tsv", sep="\t")
    chains = pd.read_csv(work / "chains.tsv", sep="\t")
    nuc_chain = {r.chain: r.gene for r in chains.itertuples() if not r.is_mt and isinstance(r.gene, str)}
    ct = pd.read_csv(contacts, sep="\t"); ct = ct[ct.partner_is_mt == 0]
    model = MMCIFParser(QUIET=True).get_structure("x", str(cif))[0]
    al = PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)

    nuc_aln, iface_cols, nuc_cols = {}, {}, {}
    for cid, g in nuc_chain.items():
        p = nucaln / f"{g}.afa"
        if not p.exists():
            print(f"  no alignment for {g}, skipped"); continue
        aln = read_fasta(p); nuc_aln[g] = aln
        hrow = aln["homo-sapiens"]; hseq = hrow.replace("-", "")
        col_of_pos = column_of_position(hrow)
        res = [r for r in model[cid] if r.id[0] == " "]
        seq = "".join(protein_letters_3to1.get(r.get_resname().upper(), "X") for r in res)
        a = al.align(seq, hseq)[0]
        m = {s0 + k: h0 + k + 1 for (s0, s1), (h0, h1) in zip(*a.aligned) for k in range(s1 - s0)}
        ident = sum(seq[i] == hseq[j - 1] for i, j in m.items()) / len(seq)
        pos_of_resnum = {r.id[1]: i for i, r in enumerate(res)}
        touching = sorted({pos_of_resnum[n] for n in ct[ct.partner_chain == cid].partner_resnum.unique() if n in pos_of_resnum})
        cols = sorted({col_of_pos[m[i]] for i in touching if i in m and m[i] in col_of_pos})
        iface_cols[g] = np.array(cols)
        nuc_cols[g] = np.array(sorted(set(col_of_pos.values()) - set(cols)))
        print(f"{g} chain {cid}: {len(seq)} residues, identity to the ortholog set's human protein {ident:.3f}; "
              f"interface columns {len(cols)}, other {len(nuc_cols[g])}")
    for g in CIII + CII + MATRIX + [x.stem for x in sorted(nucaln.glob("*.afa"))]:
        if g in nuc_aln or not (nucaln / f"{g}.afa").exists():
            continue
        nuc_aln[g] = read_fasta(nucaln / f"{g}.afa")

    mt_aln = {g: read_fasta(build / "alignments" / f"{g}.aa.fasta.gz") for g in MTG}
    tree = Phylo.read(tree_file, "newick")
    tips = {t.name: i for i, t in enumerate(tree.get_terminals())}
    H = np.zeros((len(tips), len(tips)))           # shared path length from the root

    def walk(cl, depth):
        for ch in cl.clades:
            d = depth + (ch.branch_length or 0.0)
            ids = [tips[t.name] for t in ch.get_terminals()]
            H[np.ix_(ids, ids)] = d
            if not ch.is_terminal():
                walk(ch, d)
    walk(tree.root, 0.0)

    def patristic(sp):
        i = np.array([tips[s] for s in sp]); S = H[np.ix_(i, i)]
        h = np.diag(S); return h[:, None] + h[None, :] - 2 * S

    if use_esm:
        mt_idx, mt_X = esm_load(esm, "mt13"); nu_idx, nu_X = esm_load(esm, "nuclear")

    def mean_embed(idx, X, gene, sp, cols, aln):
        """mean per-residue embedding over the alignment columns, per species (NaN row when unusable)"""
        res = np.full((len(sp), X.shape[1]), np.nan, dtype=np.float32)
        for k, s in enumerate(sp):
            key = f"{gene}|{s}"
            if key not in idx:
                continue
            off, ln, trunc = idx[key]
            row = aln[s]; ung = np.cumsum([ch != "-" for ch in row]) - 1
            if trunc or (len(row) - row.count("-")) != ln:
                continue
            take = [ung[c] for c in cols if c < len(row) and row[c] != "-"]
            if len(take) < 3:
                continue
            res[k] = np.asarray(X[off + np.array(take)], dtype=np.float32).mean(axis=0)
        got = int((~np.isnan(res[:, 0])).sum())
        if got < 0.5 * len(sp):
            print(f"    embeddings for {gene}: {got} of {len(sp)} species usable", flush=True)
        return res

    mt_contact = {g: sites[(sites.gene == g) & (sites.group == "contact")].aln_col.to_numpy() - 1 for g in MTG}
    mt_control = {g: sites[(sites.gene == g) & (sites.group == "control")].aln_col.to_numpy() - 1 for g in MTG}
    mt_by_partner = {}
    for g in MTG:
        s = sites[(sites.gene == g) & (sites.group == "contact")]
        for p in CIV:
            mt_by_partner.setdefault(p, {})[g] = s[[p in str(x).split(",") for x in s.partners]].aln_col.to_numpy() - 1

    def mt_side(sp, colsets, metric):
        """distance over species pairs, the three COX genes pooled (mean of the per-gene distances, each scaled to mean 1)"""
        Ds = []
        for g in MTG:
            cols = colsets[g]
            if len(cols) < 5:
                continue
            D = hamd(mt_aln[g], sp, cols) if metric == "seq" else cosd(mean_embed(mt_idx, mt_X, g, sp, cols, mt_aln[g]))
            Ds.append(D / np.nanmean(D))
        return np.nanmean(np.stack(Ds), axis=0)

    groups = {"complex IV partner": CIV, "complex III nuclear": CIII, "complex II (no mt subunit)": CII,
              "matrix, not respiratory": MATRIX}
    known = {g for v in groups.values() for g in v} | {"CYCS"}
    groups["cytochrome c"] = ["CYCS"]
    groups["cytosolic control"] = sorted(g for g in nuc_aln if g not in known)
    med_iface = int(np.median([len(v) for v in iface_cols.values()])) if iface_cols else 30

    rows_out = []
    for grp, genes in groups.items():
        for g in genes:
            if g not in nuc_aln:
                continue
            aln = nuc_aln[g]
            sp = sorted(set(aln) & set(mt_aln["COX1"]) & set(tips))
            if len(sp) < 60:
                continue
            iu = np.triu_indices(len(sp), 1)
            T = patristic(sp)[iu]
            if g in iface_cols:
                icols, rcols = iface_cols[g], nuc_cols[g]
            else:                                   # no interface: a random set of the usual size
                hrow = aln["homo-sapiens"]; have = np.array([c for c in range(len(hrow)) if hrow[c] != "-"])
                icols = rng.choice(have, size=min(med_iface, len(have)), replace=False); icols.sort()
                rcols = np.array(sorted(set(have) - set(icols)))
            r = dict(group=grp, gene=g, species=len(sp), iface_cols=len(icols), other_cols=len(rcols))
            for metric in ("seq", "esm"):
                if metric == "esm" and not use_esm:
                    for c in ESM_COLUMNS:
                        r[f"esm_{c}"] = np.nan
                    continue
                dn_i = (hamd(aln, sp, icols) if metric == "seq" else cosd(mean_embed(nu_idx, nu_X, g, sp, icols, aln)))[iu]
                dn_r = (hamd(aln, sp, rcols) if metric == "seq" else cosd(mean_embed(nu_idx, nu_X, g, sp, rcols, aln)))[iu]
                dm_c = mt_side(sp, mt_contact, metric)[iu]
                dm_k = mt_side(sp, mt_control, metric)[iu]
                r[f"{metric}_contact_iface"] = pcorr(dm_c, dn_i, T)
                r[f"{metric}_control_iface"] = pcorr(dm_k, dn_i, T)
                r[f"{metric}_contact_other"] = pcorr(dm_c, dn_r, T)
                r[f"{metric}_control_other"] = pcorr(dm_k, dn_r, T)
                r[f"{metric}_excess_mt"] = r[f"{metric}_contact_iface"] - r[f"{metric}_control_iface"]
                r[f"{metric}_excess_nuc"] = r[f"{metric}_contact_iface"] - r[f"{metric}_contact_other"]
                if g in mt_by_partner and sum(len(v) for v in mt_by_partner[g].values()) >= 10:
                    dm_s = mt_side(sp, mt_by_partner[g], metric)[iu]
                    r[f"{metric}_own_contacts_iface"] = pcorr(dm_s, dn_i, T)
            rows_out.append(r); print("  done", grp, g, flush=True)
    res = pd.DataFrame(rows_out)
    if not use_esm:          # keep the published column order: seq_ columns first, then esm_
        res = res[[c for c in res.columns if not c.startswith("esm_")] + [f"esm_{c}" for c in ESM_COLUMNS]]
    res.to_csv(out / "coupling.tsv", sep="\t", index=False, float_format="%.4f")

    ctrl = res[res.group == "cytosolic control"]
    print("\nTest 1 - median over the proteins of a group; Mann-Whitney p against the cytosolic controls")
    for grp in groups:
        d = res[res.group == grp]
        if d.empty:
            continue
        line = f"  {grp:28s} n={len(d):2d}"
        for col in ["seq_contact_iface", "seq_excess_mt", "seq_excess_nuc", "esm_contact_iface", "esm_excess_mt", "esm_excess_nuc"]:
            if d[col].isna().all():
                line += f"  {col.replace('seq_', 's:').replace('esm_', 'e:')} not assessed"; continue
            p = stats.mannwhitneyu(d[col].dropna(), ctrl[col].dropna())[1] if grp != "cytosolic control" and len(d) > 1 else float("nan")
            line += f"  {col.replace('seq_', 's:').replace('esm_', 'e:')} {d[col].median():+.3f} (p {p:.2g})"
        print(line)

    sp_tab = pd.read_parquet(build / "species.parquet")
    qok = set((work / "species.txt").read_text().split())
    lrows = []
    for label, keep in [("all species", None), ("Mammalia", "Mammalia"), ("with nuclear partners", "nuclear")]:
        sp = sorted(qok & set(mt_aln["COX1"]) & set(tips))
        tab = sp_tab[sp_tab.slug.isin(sp)].dropna(subset=["max_longevity_years", "adult_weight_g"])
        if keep == "Mammalia":
            tab = tab[tab["class"] == "Mammalia"]
        if keep == "nuclear":
            tab = tab[tab.slug.isin(set(nuc_aln["COX4I1"]))]
        sp = sorted(set(tab.slug))
        if len(sp) < 40:
            print(f"  {label}: {len(sp)} species, skipped"); continue
        tab = tab.set_index("slug").loc[sp]
        mc = np.zeros(len(sp)); mk = np.zeros(len(sp)); w = np.zeros(len(sp))
        for g in MTG:
            a, _ = mismatch(mt_aln[g], sp, mt_contact[g]); b, _ = mismatch(mt_aln[g], sp, mt_control[g])
            mc += a * len(mt_contact[g]); mk += b * len(mt_control[g]); w += len(mt_contact[g])
        mc /= w; mk /= w
        i = np.array([tips[s] for s in sp]); C = H[np.ix_(i, i)]
        Y = np.log10(tab.max_longevity_years.to_numpy()); mass = np.log10(tab.adult_weight_g.to_numpy())
        COV = np.column_stack([np.ones(len(sp)), mass])
        lam = float(optimize.minimize_scalar(nll, bounds=(0, 1), args=(COV, Y, C), method="bounded").x)
        for name, x in [("contact minus control", mc - mk), ("contact", mc), ("control", mk)]:
            xs = (x - x.mean()) / x.std()
            b, se, df = gls(np.column_stack([COV, xs]), Y, C, lam)
            t = b[-1] / se[-1]; q = stats.t.ppf(0.975, df)
            lrows.append(dict(set=label, species=len(sp), lam=lam, variable=name, slope_log10_years=b[-1],
                              t=t, p=2 * stats.t.sf(abs(t), df), se=se[-1],
                              ci95_low=b[-1] - q * se[-1], ci95_high=b[-1] + q * se[-1]))
        print(f"  {label}: {len(sp)} species, Pagel's lambda {lam:.3f}, mean contact mismatch {mc.mean():.4f}, control {mk.mean():.4f}")
    life = pd.DataFrame(lrows); life.to_csv(out / "lifespan.tsv", sep="\t", index=False, float_format="%.4g")
    print("\nTest 2 - phylogenetic regression of log10 maximum lifespan on the standardised mismatch, mass controlled")
    print(life.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("\nwrote", out / "coupling.tsv", out / "lifespan.tsv")
