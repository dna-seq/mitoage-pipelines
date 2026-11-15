"""H12: do the two sides of complex IV's mito-nuclear contacts change together, and does a mismatch follow lifespan?

H9 counted mismatch on the mitochondrial side only, against the consensus residue. Here each contact is a pair:
a residue of COX1-3 and a residue of a nuclear subunit within 4.5 A of it in the human enzyme (5Z62,
contacts.tsv). In a species, a pair is discordant when exactly one of its two residues differs from its column's
consensus (the modal residue among the class's species, gaps left out): one side changed, its partner did not.
The control pairs keep the nuclear residue and replace the mitochondrial one by its burial-matched non-contact
control (`mitoport burial`, which records the pairing in control_pairs.tsv), so they measure how often two residues that do not touch differ independently.

Per species: discordance over the contact pairs it fills on both sides, the same over the control pairs, and their
difference (the test). A species is kept when it fills at least half of the contact pairs. The nuclear sequences
are NCBI's vertebrate ortholog sets of the eleven subunits (`download orthologs --ortholog-taxon vertebrata`), one
protein per species (`mitoport select --class ...`), aligned by MAFFT L-INS-i (`mitoport iface-align --threads 1`).

The fit is H12's (longevity/structure_pgls.py): the class's species of H1's default dataset, log10 maximum
lifespan on the standardised variable with log10 mass and fourfold composition (and with mass only), Pagel's
lambda by maximum likelihood, the class's part of the dated tree; q within class and model.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1
from statsmodels.stats.multitest import multipletests

from ..genes import COMPLEX_IV_NUCLEAR
from ..longevity.aa_pgls import class_covariance, select
from ..longevity.structure_pgls import CLASSES, MIN_SPECIES, MODELS, fit_features
from ..seqio import read_fasta
from .burial import column_of_position

COX = ("COX1", "COX2", "COX3")
VARIABLES = ("contact minus control", "contact", "control")
MIN_SHARE = 0.5          # a species must fill both sides of at least this share of the contact pairs


def control_of(sites: pd.DataFrame, matched: pd.DataFrame) -> dict[tuple[str, int], int]:
    """{(gene, contact aln_col): its burial-matched control's aln_col}, both 1-based, from control_pairs.tsv of
    `mitoport burial`; refused unless its contacts and controls are exactly sites.tsv's."""
    for g, grp in sites.groupby("gene", sort=False):
        m = matched[matched.gene.eq(g)]
        if (set(m.contact_aln_col) != set(grp.aln_col[grp.group.eq("contact")])
                or set(m.control_aln_col) != set(grp.aln_col[grp.group.eq("control")])
                or m.control_aln_col.duplicated().any()):
            raise SystemExit(f"{g}: control_pairs.tsv does not hold sites.tsv's contacts and controls one to one; "
                             "run `mitoport burial` again to write both from one run")
    return {(r.gene, int(r.contact_aln_col)): int(r.control_aln_col) for r in matched.itertuples()}


def nuclear_columns(cif: Path, chains: pd.DataFrame, nuc_alns: dict[str, dict[str, str]]) -> dict[tuple[str, int], tuple[str, int]]:
    """{(chain, pdb residue number): (gene, 0-based column of that gene's alignment)} for every nuclear chain of the
    structure whose gene has an alignment holding the human protein; residues are placed on the human row by a
    global alignment of the chain's sequence to it (as `mitoport iface` does)."""
    model = MMCIFParser(QUIET=True).get_structure("x", str(cif))[0]
    al = PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5)
    out = {}
    for r in chains.itertuples():
        if r.is_mt or not isinstance(r.gene, str) or r.gene not in nuc_alns or "homo-sapiens" not in nuc_alns[r.gene]:
            continue
        hrow = nuc_alns[r.gene]["homo-sapiens"]; hseq = hrow.replace("-", ""); col = column_of_position(hrow)
        res = [x for x in model[r.chain] if x.id[0] == " "]
        seq = "".join(protein_letters_3to1.get(x.get_resname().upper(), "X") for x in res)
        a = al.align(seq, hseq)[0]
        m = {s0 + k: h0 + k + 1 for (s0, s1), (h0, h1) in zip(*a.aligned) for k in range(s1 - s0)}
        out.update({(r.chain, res[i].id[1]): (r.gene, col[hp]) for i, hp in m.items() if hp in col})
    return out


def pairs(contacts: pd.DataFrame, sites: pd.DataFrame, nuc_col: dict, control: dict) -> pd.DataFrame:
    """One row per distinct (mito column, nuclear column) contact pair and its control pair: mt_gene, mt_col,
    ctl_col (0-based mito columns), nuc_gene, nuc_col (0-based)."""
    ct = contacts[contacts.partner_is_mt.eq(0)]
    st = sites.set_index(["gene", "pdb_resnum"])
    rows = []
    for r in ct.itertuples():
        s = st.loc[(r.mt_gene, r.pdb_resnum)]
        n = nuc_col.get((r.partner_chain, r.partner_resnum))
        ctl = control.get((r.mt_gene, int(s.aln_col)))
        if s.group != "contact" or s.aln_col <= 0 or n is None or ctl is None:
            continue
        rows.append({"mt_gene": r.mt_gene, "mt_col": int(s.aln_col) - 1, "ctl_col": ctl - 1, "nuc_gene": n[0], "nuc_col": n[1]})
    return pd.DataFrame(rows).drop_duplicates().sort_values(["mt_gene", "mt_col", "nuc_gene", "nuc_col"]).reset_index(drop=True)


def consensus(aln: dict[str, str], slugs: list[str], col: int) -> str:
    """The modal residue of the column among `slugs` that hold the alignment, gaps out ('-' if none);
    ties go to the first in sorted order (as `iface.mismatch`)."""
    v = np.array([aln[s][col] for s in slugs if s in aln and aln[s][col] != "-"])
    if len(v) == 0:
        return "-"
    u, n = np.unique(v, return_counts=True)
    return str(u[n.argmax()])


def discordance(pr: pd.DataFrame, mt_alns: dict, nuc_alns: dict, slugs: list[str]) -> pd.DataFrame:
    """Per species of `slugs`: share of discordant pairs over the contact pairs it fills on both sides, the same over
    the control pairs, their difference, and how many contact pairs it fills. Consensus is taken among `slugs`."""
    cons = {}
    for r in pr.itertuples():
        for key, aln, col in (((r.mt_gene, r.mt_col), mt_alns[r.mt_gene], r.mt_col),
                              ((r.mt_gene, r.ctl_col), mt_alns[r.mt_gene], r.ctl_col),
                              ((r.nuc_gene, r.nuc_col), nuc_alns[r.nuc_gene], r.nuc_col)):
            if key not in cons:
                cons[key] = consensus(aln, slugs, col)
    rows = []
    for s in slugs:
        out = {"slug": s}
        for name, mcol in (("contact", "mt_col"), ("control", "ctl_col")):
            d = n = 0
            for r in pr.itertuples():
                nuc = nuc_alns[r.nuc_gene].get(s)
                if nuc is None or s not in mt_alns[r.mt_gene]:
                    continue
                c = getattr(r, mcol); a, b = mt_alns[r.mt_gene][s][c], nuc[r.nuc_col]
                if a == "-" or b == "-":
                    continue
                n += 1; d += (a != cons[(r.mt_gene, c)]) != (b != cons[(r.nuc_gene, r.nuc_col)])
            out[name] = d / n if n else np.nan
            out[f"{name}_pairs"] = n
        rows.append(out)
    res = pd.DataFrame(rows).set_index("slug")
    res.loc[res.contact_pairs < MIN_SHARE * len(pr), ["contact", "control"]] = np.nan
    res["contact minus control"] = res.contact - res.control
    return res


def main(cif: Path, data: Path, work: Path, nucaln: Path, composition: Path, fourfold: Path, tree_file: Path, out: Path,
         classes: tuple[str, ...] = CLASSES) -> pd.DataFrame:
    """Writes partners.tsv (the fits) and partners_species.tsv (every species' discordance and filled pairs)."""
    missing = [g for g in COMPLEX_IV_NUCLEAR if not (nucaln / f"{g}.afa").is_file()]
    if len(missing) == len(COMPLEX_IV_NUCLEAR) or not (work / "sites.tsv").is_file() or not (work / "control_pairs.tsv").is_file():
        raise SystemExit(f"needs the contact map ({work}/sites.tsv and control_pairs.tsv, from `mitoport burial`) and the vertebrate partner "
                         f"alignments ({nucaln}/<GENE>.afa, from `download orthologs --ortholog-taxon vertebrata`, "
                         "`mitoport select --class ...`, `mitoport iface-align`)")
    if missing:
        print(f"NOTE: no alignment for {' '.join(missing)}: their pairs are left out")
    sites = pd.read_csv(work / "sites.tsv", sep="\t")
    nuc_alns = {g: read_fasta(nucaln / f"{g}.afa") for g in COMPLEX_IV_NUCLEAR if g not in missing}
    mt_alns = {g: read_fasta(data / "alignments" / f"{g}.aa.fasta.gz") for g in COX}
    pr = pairs(pd.read_csv(work / "contacts.tsv", sep="\t"), sites,
               nuclear_columns(cif, pd.read_csv(work / "chains.tsv", sep="\t"), nuc_alns),
               control_of(sites, pd.read_csv(work / "control_pairs.tsv", sep="\t")))
    print(f"{len(pr)} contact pairs over {pr.mt_gene.nunique()} mitochondrial and {pr.nuc_gene.nunique()} nuclear genes: "
          + ", ".join(f"{g} {n}" for g, n in pr.nuc_gene.value_counts().sort_index().items()))
    out.mkdir(parents=True, exist_ok=True)
    rows, per = [], []
    for cls in classes:
        sp, _ = select(data, composition, fourfold, tree_file, cls)
        with_partner = [s for s in sp.slug if any(s in a for a in nuc_alns.values())]
        print(f"{cls}: {len(sp)} species in H12's set, {len(with_partner)} with a partner sequence")
        if len(with_partner) < MIN_SPECIES:
            rows += [{"class": cls, "family": "partners", "feature": "*", "model": m, "n": len(with_partner),
                      "status": f"not assessed: {len(with_partner)} species with a partner sequence"} for m in MODELS]
            continue
        C, order, _ = class_covariance(tree_file, list(sp.slug))
        sp = sp.set_index("slug").loc[order].reset_index()
        d = discordance(pr, mt_alns, nuc_alns, with_partner)
        print(f"{cls}: {int(d.contact.notna().sum())} species fill half the contact pairs; mean discordance contact "
              f"{d.contact.mean():.4f}, control {d.control.mean():.4f}")
        per.append(d.assign(**{"class": cls}).reset_index())
        rows += fit_features(sp, d.reindex(sp.slug)[list(VARIABLES)], C, cls, "partners")
    res = pd.DataFrame(rows)
    res["q"] = np.nan
    t = res.status.eq("tested")
    for _, idx in res[t].groupby(["class", "model"]).groups.items():
        res.loc[idx, "q"] = multipletests(res.loc[idx, "p"], method="fdr_bh")[1]
    res.to_csv(out / "partners.tsv", sep="\t", index=False, float_format="%.5g")
    if per:
        pd.concat(per).to_csv(out / "partners_species.tsv", sep="\t", index=False, float_format="%.5g")
    with pd.option_context("display.width", 220):
        cols = [c for c in ["class", "feature", "model", "n", "status", "lambda", "slope_per_sd", "ci95_low", "ci95_high", "p", "q"] if c in res]
        print(res[cols].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    return res
