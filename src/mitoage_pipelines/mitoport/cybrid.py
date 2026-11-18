"""Simulated xenomitochondrial cybrids (H16): does a mismatch score at the mito-nuclear contacts of a solved
complex order the published cybrids as their measured complex activities do?

A cybrid puts the mtDNA of a donor species into cells of a host species: the 13 mtDNA-encoded subunits come
from the donor, every nuclear-encoded subunit from the host. The simulation builds that chimera on the solved
human complexes (I 5XTD, III 5XTE, IV 5Z62, V 8H9S) and counts, for each complex, the positions of the
mtDNA-encoded subunits that touch a nuclear subunit (any two heavy atoms within 4.5 A) where the donor's
residue differs from the host's. Contacts between two mtDNA-encoded subunits are not mismatches in a cybrid
(both come from the donor) and are not counted. The same count over the resolved non-contact positions and over
every position of the proteins is the baseline: divergence alone. Positions are human positions; host and donor
proteins are placed on them by a global alignment to the human protein of the same gene.

The measured side is `cybrid_measurements.tsv` in this folder: one row per (host, donor, complex, source),
transcribed by hand from the papers named in it, with the figure or table each value is read from. A value the
paper does not state is left empty; `significant` keeps the paper's own reading.

The tests, fixed before the scores were computed:
  T1  across (host, donor, complex) with a stated activity: Spearman of activity against each score;
      support for the interface reading needs the contact score to order them better than divergence does.
  T2  across every (host, donor, complex) the papers call reduced or not reduced: how well each score separates
      the two (the area under the ROC curve; 0.5 is chance). A (host, donor, complex) whose papers disagree is
      left out and listed.
  T3  within each pair with at least one reduced complex: is the most reduced complex the one with the
      highest score?
  T4  for donors whose mtDNA gave no respiring cybrid: is their total contact mismatch above every donor that did,
      in the same host?
Every comparison has a handful of points; the output says how many, and nothing here is a significance test.
T2b (T2 within each complex) was added after the first run had been read, and is labelled exploratory.
In H17's readout (`readout`) T2b is fixed in advance like the rest; its rows keep the wording of H16's table.

Genomes: RefSeq records at fixed versions (`GENOMES`). Where a paper deposited the sequence of the cell line it
used, that sequence replaces the RefSeq one for its gene (`OWN_LINE`); such records are partial at their ends,
and the human positions they do not cover are unknown, never differences.
"""
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import protein_letters_3to1
from Bio.SeqFeature import AfterPosition, BeforePosition
from scipy import stats

from ..genes import MT_GENES
from . import complexes

CUTOFF = 4.5
MEASUREMENTS = Path(__file__).with_name("cybrid_measurements.tsv")
STRUCTURES = {k: c.pdb for k, c in complexes.COMPLEXES.items()}          # I 5XTD, III 5XTE, IV 5Z62, V 8H9S
COMPLEX_GENES = {k: tuple(c.mt) for k, c in complexes.COMPLEXES.items()}
GENOMES = {
    "homo-sapiens": "NC_012920.1", "pan-troglodytes": "NC_001643.1", "pan-paniscus": "NC_001644.1",
    "gorilla-gorilla": "NC_011120.1", "pongo-pygmaeus": "NC_001646.1", "pongo-abelii": "NC_002083.1",
    "mus-musculus": "NC_005089.1", "mus-spretus": "NC_025952.1", "mus-caroli": "NC_025268.1",
    "mus-terricolor": "NC_010650.1", "mus-pahari": "NC_036680.1", "rattus-norvegicus": "NC_001665.2",
    "otomys-irroratus": "NC_049116.1", "cricetulus-griseus": "NC_007936.1", "mesocricetus-auratus": "NC_013276.1",
}
# cytochrome b sequenced from the cell lines of McKenzie et al. 2003 (doi 10.1093/molbev/msg132, table 1)
OWN_LINE = {("mus-spretus", "CYTB"): "AY224678.1", ("mus-terricolor", "CYTB"): "AY193770.1",
            ("otomys-irroratus", "CYTB"): "AY224677.1"}
HUMAN = "homo-sapiens"
SCORES = ("d_contact", "f_contact", "d_contact_radical", "d_noncontact", "f_noncontact", "d_all", "f_all")
NAMES = {"COI": "COX1", "COII": "COX2", "COIII": "COX3", "CO1": "COX1", "CO2": "COX2", "CO3": "COX3",
         "COB": "CYTB", "CYB": "CYTB", "CYT B": "CYTB", "NAD1": "ND1", "NAD2": "ND2", "NAD3": "ND3", "NAD4": "ND4",
         "NAD4L": "ND4L", "NAD5": "ND5", "NAD6": "ND6", "ATPASE6": "ATP6", "ATPASE8": "ATP8"}
PRODUCTS = {"cytochrome b": "CYTB"}
ALL_GENES = tuple(g for gs in COMPLEX_GENES.values() for g in gs)


def accessions() -> list[str]:
    return sorted(set(GENOMES.values()) | set(OWN_LINE.values()))


def gene_name(gene: str | None, product: str | None) -> str | None:
    """One of the 13 gene names, from a CDS's /gene (any of the usual spellings) or, failing that, its /product."""
    if gene:
        g = gene.upper().removeprefix("MT-").replace("-", "")
        g = NAMES.get(g, g)
        if g in ALL_GENES:
            return g
    if product:
        return PRODUCTS.get(product.strip().lower())
    return None


def read_proteins(gb: Path) -> dict[str, dict[str, tuple[str, bool]]]:
    """{accession.version: {gene: (protein, partial)}} from the CDS translations; partial when a CDS end is fuzzy."""
    out = {}
    for rec in SeqIO.parse(str(gb), "genbank"):
        genes = {}
        for f in rec.features:
            if f.type != "CDS" or "translation" not in f.qualifiers:
                continue
            g = gene_name(f.qualifiers.get("gene", [None])[0], f.qualifiers.get("product", [None])[0])
            if g is None:
                continue
            if g in genes:
                raise ValueError(f"{rec.id}: two CDS for {g}")
            partial = isinstance(f.location.start, BeforePosition) or isinstance(f.location.end, AfterPosition)
            genes[g] = (f.qualifiers["translation"][0], partial)
        out[rec.id] = genes
    return out


def species_proteins(records: dict, slug: str) -> dict[str, tuple[str, bool, str]]:
    """{gene: (protein, partial, accession)} of one species: its RefSeq genome, with the cell line's own sequence
    in place of a gene where one is listed in OWN_LINE. Refuses when a gene is missing."""
    acc = GENOMES[slug]
    if acc not in records:
        raise KeyError(f"{acc} ({slug}) is not in the GenBank file: rerun `mitoage download cybrid`")
    out = {g: (p, part, acc) for g, (p, part) in records[acc].items()}
    for (s, g), a in OWN_LINE.items():
        if s == slug:
            if g not in records.get(a, {}):
                raise KeyError(f"{a} has no {g} CDS: rerun `mitoage download cybrid`")
            p, part = records[a][g]
            out[g] = (p, part, a)
    missing = [g for g in ALL_GENES if g not in out]
    if missing:
        raise ValueError(f"{slug} ({acc}) lacks {missing}")
    return out


def aligner() -> PairwiseAligner:
    al = PairwiseAligner()
    al.mode = "global"; al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -10; al.extend_gap_score = -0.5
    return al


def on_human(human: str, other: str, partial: bool, al: PairwiseAligner | None = None) -> list[str]:
    """`other`'s residue at each human position (index 0 is human position 1): a letter, '-' for a deletion in
    `other`, '?' for a position a partial sequence does not cover (unknown, not a difference)."""
    a = (al or aligner()).align(human, other)[0]
    out = ["-"] * len(human)
    first = last = None
    for (h0, h1), (o0, o1) in zip(*a.aligned):
        for k in range(h1 - h0):
            out[h0 + k] = other[o0 + k]
        first = h0 if first is None else first
        last = h1
    if partial and first is not None:
        for i in itertools.chain(range(0, first), range(last, len(human))):
            out[i] = "?"
    return out


# --- the contact maps -------------------------------------------------------------------------------------
def contact_map(cif: Path, complex_: str, human: dict[str, str]) -> pd.DataFrame:
    """One row per human position of the complex's mtDNA-encoded genes: resolved in the reference chain or not,
    whether it touches a nuclear-encoded chain, the number of distinct nuclear residues it touches, and whether it
    touches another mtDNA-encoded chain. Chains are named as `mitoport contacts` names them (`complexes.py`, by their
    UniProt references): each gene is read from its reference chain, and a second copy of a mitochondrial protein
    (complex III's other cytochrome b) is a mitochondrial partner, never a nuclear one."""
    cx = complexes.get(complex_)
    gene_of_chain = complexes.chain_genes(cif, cx)
    model = MMCIFParser(QUIET=True).get_structure(cif.stem, str(cif))[0]
    is_mt = {c for c, g in gene_of_chain.items() if g in MT_GENES}
    al = aligner(); mt = {}
    for g, c in cx.mt.items():
        res = [r for r in model[c] if r.id[0] == " "]
        seq = "".join(protein_letters_3to1.get(r.get_resname().upper(), "X") for r in res)
        a = al.align(seq, human[g])[0]
        mt[c] = (g, {res[s0 + k].id: h0 + k + 1 for (s0, s1), (h0, h1) in zip(*a.aligned) for k in range(s1 - s0)})
    atoms = [a for c in model for r in c if r.id[0] == " " for a in r if a.element != "H"]
    nuc_partners, mt_touch = {}, set()
    for a, b in NeighborSearch(atoms).search_all(CUTOFF):
        ra, rb = a.get_parent(), b.get_parent()
        ca, cb = ra.get_parent().id, rb.get_parent().id
        if ca == cb:
            continue
        for (r1, c1), (r2, c2) in (((ra, ca), (rb, cb)), ((rb, cb), (ra, ca))):
            if c1 not in mt or r1.id not in mt[c1][1]:
                continue
            key = (mt[c1][0], mt[c1][1][r1.id])
            if c2 in is_mt:
                mt_touch.add(key)
            else:
                nuc_partners.setdefault(key, set()).add((c2, r2.id[1]))
    resolved = {(g, p) for g, m in mt.values() for p in m.values()}
    rows = []
    for g in COMPLEX_GENES[complex_]:
        for p in range(1, len(human[g]) + 1):
            k = (g, p); partners = nuc_partners.get(k, set())
            rows.append({"complex": complex_, "pdb": cif.stem, "gene": g, "human_pos": p, "human_aa": human[g][p - 1],
                         "resolved": int(k in resolved), "contact_nuclear": int(bool(partners)),
                         "nuclear_partner_residues": len(partners),
                         "nuclear_chains": ",".join(sorted({c for c, _ in partners})),
                         "contact_mt": int(k in mt_touch)})
    return pd.DataFrame(rows)


# --- scores -----------------------------------------------------------------------------------------------
def pair_scores(positions: pd.DataFrame, human: dict[str, str], host: dict, donor: dict,
                blosum=None, al: PairwiseAligner | None = None) -> dict:
    """The counts for one (host, donor) on one complex's positions. A position is compared when both residues
    are known; a deletion against a residue is a difference."""
    blosum = blosum if blosum is not None else substitution_matrices.load("BLOSUM62")
    al = al or aligner()
    c = dict(n_contact=0, d_contact=0, d_contact_radical=0, n_noncontact=0, d_noncontact=0, n_all=0, d_all=0,
             n_unknown=0)
    for g, grp in positions.groupby("gene", sort=False):
        hp, hpart, _ = host[g]; dp, dpart, _ = donor[g]
        hs = on_human(human[g], hp, hpart, al); ds = on_human(human[g], dp, dpart, al)
        for p, resolved, contact in zip(grp.human_pos, grp.resolved, grp.contact_nuclear):
            a, b = hs[p - 1], ds[p - 1]
            if a == "?" or b == "?":
                c["n_unknown"] += 1; continue
            diff = a != b
            c["n_all"] += 1; c["d_all"] += diff
            if not resolved:
                continue
            if contact:
                c["n_contact"] += 1; c["d_contact"] += diff
                radical = diff and ("-" in (a, b) or blosum[a][b] < 0)
                c["d_contact_radical"] += radical
            else:
                c["n_noncontact"] += 1; c["d_noncontact"] += diff
    for k in ("contact", "noncontact", "all"):
        c[f"f_{k}"] = c[f"d_{k}"] / c[f"n_{k}"] if c[f"n_{k}"] else np.nan
    return c


def read_measurements(path: Path = MEASUREMENTS) -> pd.DataFrame:
    m = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for col in ("pct_low", "pct_high"):
        m[col] = pd.to_numeric(m[col].replace("", np.nan))
    m["pct"] = (m.pct_low + m.pct_high) / 2
    bad = set(m.outcome) - {"measured", "no_respiring_cybrid"}
    if bad:
        raise ValueError(f"unknown outcome in {path.name}: {bad}")
    return m


def status(rows: pd.DataFrame) -> str:
    """reduced, not_reduced or disagree for the rows of one (host, donor, complex): reduced when the paper calls
    the decrease significant, not reduced when it does not; papers that differ give `disagree`."""
    s = {"reduced" if r.significant == "yes" and not (r.pct >= 100) else "not_reduced" for r in rows.itertuples()}
    return s.pop() if len(s) == 1 else "disagree"


def auc(score: pd.Series, reduced: pd.Series) -> float:
    """Area under the ROC curve of `score` for telling reduced from not reduced (ties count one half)."""
    pos, neg = score[reduced].to_numpy(float), score[~reduced].to_numpy(float)
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    return float(sum((p > n) + 0.5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg)))


T2B_SCORES = ("d_contact", "f_contact", "d_noncontact", "d_all", "f_all")
T3_SCORES = ("f_contact", "d_contact", "f_all", "d_all")
T4_SCORES = ("d_contact", "d_all")


def complete_totals(scores: pd.DataFrame, names: tuple[str, ...],
                    complexes_by_host: dict[str, tuple[str, ...]] | None = None) -> pd.DataFrame:
    """Sum a score only when every required complex has a finite value for that donor.

    H16/H17 require I, III, IV and V. H18 explicitly supplies its host-specific structure sets.
    """
    rows = []
    for (host, donor), group in scores.groupby(["host", "donor"], sort=False):
        required = (complexes_by_host or {}).get(host, tuple(STRUCTURES))
        values = group.set_index("complex")[list(names)].reindex(required)
        values = values.where(np.isfinite(values))
        total = values.sum(min_count=len(required))
        rows.append({"host": host, "donor": donor, **total.to_dict()})
    return pd.DataFrame(rows, columns=["host", "donor", *names])


def summarise(scores: pd.DataFrame, meas: pd.DataFrame, names: tuple[str, ...] = SCORES,
              t2b: tuple[str, ...] = T2B_SCORES, t3: tuple[str, ...] = T3_SCORES,
              t4: tuple[str, ...] = T4_SCORES,
              t4_complexes: dict[str, tuple[str, ...]] | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(calibration: one row per measured (host, donor, complex) with its status and scores; summary: T1-T4).
    The defaults are H16's scores; H17 passes its own. A score that is missing for a point is left out of the
    readouts it enters, and each row's `n` says how many points it had."""
    key = ["host", "donor", "complex"]
    measured = meas[meas.outcome == "measured"]
    cal = (measured.groupby(key, sort=False)
           .apply(lambda g: pd.Series({"status": status(g), "pct": g.pct.mean(), "n_sources": g.source_doi.nunique(),
                                       "sources": ";".join(sorted(set(g.source_doi)))}), include_groups=False)
           .reset_index().merge(scores, on=key, how="left"))
    if cal[list(names)].isna().all(axis=1).any():
        raise ValueError("a measured pair has no score")
    out = []
    stated = cal[cal.pct.notna()]
    for s in names:
        st = stated[stated[s].notna()]
        if len(st) >= 3:
            rho, _ = stats.spearmanr(st.pct, st[s])
        else:
            rho = np.nan
        out.append({"test": "T1", "score": s, "n": len(st), "value": rho,
                    "reads": "Spearman of activity (% of control) against the score; negative is the prediction"})
    use = cal[cal.status != "disagree"]
    red = use.status == "reduced"
    for s in names:
        ok = use[s].notna()
        out.append({"test": "T2", "score": s, "n": f"{int((red & ok).sum())} reduced, {int((~red & ok).sum())} not",
                    "value": auc(use[s][ok], red[ok]), "reads": "area under the ROC curve for reduced against not reduced"})
    # T2b, added after the first run had been read: T2 within one complex at a time, where complex size cannot
    # carry the score. Exploratory, and labelled so in the table.
    for cx, g in use.groupby("complex", sort=False):
        r = g.status == "reduced"
        for s in t2b:
            ok = g[s].notna()
            out.append({"test": "T2b", "score": f"{s} (complex {cx})", "n": f"{int((r & ok).sum())} reduced, {int((~r & ok).sum())} not",
                        "value": auc(g[s][ok], r[ok]), "reads": "T2 within one complex; added after the first run (exploratory)"})
    for r in cal[cal.status == "disagree"].itertuples():
        out.append({"test": "T2", "score": "left out", "n": "", "value": np.nan,
                    "reads": f"{r.host} / {r.donor} complex {r.complex}: the sources disagree ({r.sources})"})
    for (h, d), g in cal.groupby(["host", "donor"], sort=False):
        g = g[g.status != "disagree"]
        if (g.status == "reduced").sum() == 0 or len(g) < 2:
            continue
        # most reduced: the lowest stated activity among the reduced complexes, else the only reduced one
        r = g[g.status == "reduced"]
        worst = r.loc[r.pct.idxmin(), "complex"] if r.pct.notna().any() else (r.complex.iloc[0] if len(r) == 1 else None)
        for s in t3:
            if g[s].isna().any():
                out.append({"test": "T3", "score": s, "n": len(g), "value": np.nan,
                            "reads": f"{h} / {d}: not assessed, the score is missing for a complex"})
                continue
            top = g.loc[g[s].idxmax(), "complex"]
            out.append({"test": "T3", "score": s, "n": len(g), "value": float(top == worst) if worst else np.nan,
                        "reads": f"{h} / {d}: most reduced {worst or 'unknown'}, highest {s} in complex {top}"})
    total = complete_totals(scores, t4, t4_complexes)
    failed = meas[meas.outcome == "no_respiring_cybrid"][["host", "donor"]].drop_duplicates()
    respiring = measured[["host", "donor"]].drop_duplicates()
    for r in failed.itertuples():
        for s in t4:
            mine = total.set_index(["host", "donor"])[s].get((r.host, r.donor), np.nan)
            others = respiring[respiring.host == r.host].merge(total, on=["host", "donor"], how="left")[s]
            known = pd.notna(mine) and len(others) > 0 and others.notna().all()
            out.append({"test": "T4", "score": s, "n": len(others), "value": float(mine > others.max()) if known else np.nan,
                        "reads": f"{r.host} / {r.donor}: total {s} {mine:g} against at most {others.max():g} "
                                 f"for a donor that respired" + ("" if known else " (not assessed: a total is missing)")})
    return cal, pd.DataFrame(out)


# --- H17: the same readouts for the co-evolution and language-model scores --------------------------------
H17_PRIMARY = ("coev_excess", "coev_excess_mean", "esm_contact", "esm_contact_mean", "esm_partner_excess")
H18_PRIMARY = ("ddE", "ddE_per_change")
H17_OTHER = ("coev_contact", "coev_control", "coev_excess_t50", "coev_excess_noorder", "esm_control", "esm_all",
             "esm_partner_effect", "esm_partner_effect_control")
H17_BASELINE = ("d_contact", "d_all", "f_all")


def readout(out: Path, measurements: Path = MEASUREMENTS) -> pd.DataFrame:
    """H16's readouts T1 to T4 for H17's scores, beside H16's baselines on the same points. Reads scores.tsv
    (`mitoport cybrid`), coev_scores.tsv (`mitoport cybrid-coev`) and esm_scores.tsv (`mitoport cybrid-esm`) in
    `out`; a table that is missing is said loudly and its scores are absent from the summary, never zero.
    Writes h17_calibration.tsv and h17_summary.tsv."""
    key = ["host", "donor", "complex"]
    if not (out / "scores.tsv").is_file():
        raise SystemExit(f"{out / 'scores.tsv'} is missing: run `mitoage mitoport cybrid` first")
    scores = pd.read_csv(out / "scores.tsv", sep="\t")[key + list(H17_BASELINE)]
    have = []
    for f in ("coev_scores.tsv", "esm_scores.tsv"):
        if not (out / f).is_file():
            print(f"NOT ASSESSED: {out / f} is missing; its scores are left out of the summary")
            continue
        t = pd.read_csv(out / f, sep="\t")
        cols = [c for c in H17_PRIMARY + H17_OTHER if c in t.columns]
        scores = scores.merge(t[key + cols], on=key, how="left"); have += cols
    if not have:
        raise SystemExit("neither coev_scores.tsv nor esm_scores.tsv is there: nothing to read out")
    primary = tuple(c for c in H17_PRIMARY if c in have)
    names = primary + tuple(c for c in H17_OTHER if c in have) + H17_BASELINE
    cal, summary = summarise(scores, read_measurements(measurements), names, t2b=primary + H17_BASELINE,
                             t3=primary + H17_BASELINE, t4=primary + ("d_contact", "d_all"))
    cal.to_csv(out / "h17_calibration.tsv", sep="\t", index=False, float_format="%.4f")
    summary.to_csv(out / "h17_summary.tsv", sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.width", 220, "display.max_colwidth", 110, "display.max_rows", 400):
        print(summary.to_string(index=False))
    print("wrote", out / "h17_calibration.tsv", out / "h17_summary.tsv")
    return summary


def main(gb: Path, pdb_dir: Path, out: Path, measurements: Path = MEASUREMENTS) -> None:
    out.mkdir(parents=True, exist_ok=True)
    records = read_proteins(gb)
    meas = read_measurements(measurements)
    species = sorted(set(meas.host) | set(meas.donor))
    prot = {s: species_proteins(records, s) for s in species + ([HUMAN] if HUMAN not in species else [])}
    human = {g: p for g, (p, _, _) in prot[HUMAN].items()}

    maps = []
    for cx, pid in STRUCTURES.items():
        cif = pdb_dir / f"{pid}.cif"
        if not cif.is_file():
            raise FileNotFoundError(f"{cif} is missing: `mitoage download pdb {' '.join(STRUCTURES.values())}`")
        m = contact_map(cif, cx, human); maps.append(m)
        res = m[m.resolved == 1]
        print(f"complex {cx} ({pid}): {len(res)} of {len(m)} positions resolved, {int(res.contact_nuclear.sum())} touch "
              f"a nuclear subunit, {int(res.contact_mt.sum())} another mtDNA-encoded one")
    positions = pd.concat(maps, ignore_index=True)
    positions.to_csv(out / "positions.tsv", sep="\t", index=False)

    al, blosum = aligner(), substitution_matrices.load("BLOSUM62")
    rows = []
    pairs = meas[["host", "donor"]].drop_duplicates().itertuples(index=False)
    for h, d in sorted(pairs):
        for cx in STRUCTURES:
            sc = pair_scores(positions[positions.complex == cx], human, prot[h], prot[d], blosum, al)
            rows.append({"host": h, "donor": d, "complex": cx, **sc})
    scores = pd.DataFrame(rows)
    sources = pd.DataFrame([{"species": s, "gene": g, "accession": a, "partial": int(part)}
                            for s in sorted(prot) for g, (_, part, a) in sorted(prot[s].items())])
    sources.to_csv(out / "sequences.tsv", sep="\t", index=False)
    scores.to_csv(out / "scores.tsv", sep="\t", index=False, float_format="%.4f")
    cal, summary = summarise(scores, meas)
    cal.to_csv(out / "calibration.tsv", sep="\t", index=False, float_format="%.4f")
    summary.to_csv(out / "summary.tsv", sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.width", 200, "display.max_colwidth", 120, "display.max_rows", 200):
        print(cal[["host", "donor", "complex", "status", "pct", "d_contact", "f_contact", "d_all", "f_all"]].to_string(index=False))
        print(summary.to_string(index=False))
    print("wrote", out)


def readout_energy(out: Path, measurements: Path = MEASUREMENTS) -> pd.DataFrame:
    """H16's readouts T1 to T4 for H18's interface energy, beside H16's baselines on the same points. Reads scores.tsv
    (`mitoport cybrid`) and energy_scores.tsv (`mitoport cybrid-energy-collect`). Only rows whose status is `computed`
    or `no substitution to model` carry a score; a point without one is left out of the readouts it would enter.
    Writes h18_calibration.tsv and h18_summary.tsv."""
    from .cybrid_energy import STRUCTURES as ENERGY_STRUCTURES
    key = ["host", "donor", "complex"]
    for f in ("scores.tsv", "energy_scores.tsv"):
        if not (out / f).is_file():
            raise SystemExit(f"{out / f} is missing")
    e = pd.read_csv(out / "energy_scores.tsv", sep="\t")
    scores = pd.read_csv(out / "scores.tsv", sep="\t")[key + list(H17_BASELINE)].merge(e[key + list(H18_PRIMARY)], on=key, how="left")
    names = H18_PRIMARY + H17_BASELINE
    cal, summary = summarise(scores, read_measurements(measurements), names, t2b=names, t3=names,
                             t4=("ddE", "d_contact", "d_all"),
                             t4_complexes={h: tuple(cs) for h, cs in ENERGY_STRUCTURES.items()})
    cal = cal.merge(e[key + [c for c in ("ddE_sd", "substituted", "interface_substituted", "unmodelled", "status") if c in e]], on=key, how="left")
    cal.to_csv(out / "h18_calibration.tsv", sep="\t", index=False, float_format="%.4f")
    summary.to_csv(out / "h18_summary.tsv", sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.width", 220, "display.max_colwidth", 110, "display.max_rows", 400):
        print(summary.to_string(index=False))
    print("wrote", out / "h18_calibration.tsv", out / "h18_summary.tsv")
    return summary
