"""H17, the co-evolution score: does the chimeric pairing of a cybrid lose residue pairs that vertebrates keep together?

H16 counted differing contact residues, every substitution alike. Here each contact is a pair, as in H12 and H14: a
residue of an mtDNA-encoded subunit and the nuclear residue within 4.5 A of it in the human complex. Over the
vertebrates that have both proteins (the build's mitochondrial alignments and the vertebrate partner alignments of
`mitoport iface-align`), each pair of columns gets the joint distribution of its two residues, and from it the
pointwise mutual information

    PMI(a, b) = log P(a, b) / (P(a) P(b))

of a mitochondrial residue a with a nuclear residue b: positive when vertebrates carry the two together more often
than their frequencies alone give. In a cybrid the nuclear residue is the host's and the mitochondrial one the
donor's, so a pair where donor and host differ on the mitochondrial side loses

    loss = PMI(host mt, host nuclear) - PMI(donor mt, host nuclear)

and a pair where they do not differ loses nothing. `coev_contact` is the sum of the loss over a complex's contact
pairs. `coev_control` is the same over the control pairs (the nuclear residue kept, the mitochondrial one replaced
by its burial-matched non-contact control of `mitoport burial`): what two residues that do not touch give, through
shared ancestry alone. `coev_excess` is their difference, the score; `coev_excess_mean` divides it by the number of
contact pairs scored, so that a large complex does not score high by its size.

Fixed here, before any score was computed:
- Sequences are weighted by 1 / (the number of sequences at least 80% identical to them over the complex's
  concatenated mitochondrial alignment, themselves included), as in direct-coupling analyses.
- P(a, b) = (1 - THETA) f(a, b) + THETA / 400 over the 20 amino acids, f the weighted frequency among the species
  with a residue on both sides; THETA = 0.05; the marginals are those of this P. THETA = 0.5 is run beside it as a
  sensitivity (`_t50` columns) and enters no reading.
- The host and the donor are left out of the counts of the pairs they are scored on. As a sensitivity the host's
  whole order (Primates for a human host, Rodentia for a mouse host) is left out (`_noorder` columns).
- A pair is scored when host mt, donor mt and host nuclear residues are all amino acids; the others are counted in
  `pairs_unscored` and never read as agreement. A pair with fewer than MIN_SPECIES species on both sides is unscored.

Writes coev_scores.tsv (one row per host, donor and complex) and coev_pairs.tsv (the pairs of every complex, with
the number of species and the effective number of sequences behind each).
"""
from pathlib import Path

import numpy as np
import pandas as pd

from ..seqio import read_fasta
from . import complexes, cybrid
from .burial import column_of_position
from .partners import control_of, nuclear_columns, pairs

AA = "ACDEFGHIKLMNPQRSTVWY"
IDX = {a: i for i, a in enumerate(AA)}
THETA = 0.05
THETA_SENSITIVITY = 0.5
IDENTITY = 0.8
MIN_SPECIES = 30
HOST_ORDER = {"homo-sapiens": "Primates", "mus-musculus": "Rodentia"}


def weights(rows: dict[str, str], identity: float = IDENTITY) -> dict[str, float]:
    """{name: 1 / number of rows at least `identity` identical to it, itself included}; identity is the share of
    equal letters over the columns where both rows have a residue (0 when they share none)."""
    names = list(rows)
    if not names:
        return {}
    m = np.array([np.frombuffer(rows[n].encode(), dtype=np.uint8) for n in names])
    res = m != ord("-")
    out = {}
    for i, n in enumerate(names):
        both = res & res[i]
        same = ((m == m[i]) & both).sum(axis=1)
        tot = both.sum(axis=1)
        ident = np.divide(same, tot, out=np.zeros(len(names)), where=tot > 0)
        out[n] = 1.0 / max(1, int((ident >= identity).sum()))
    return out


def joint(a_col: dict[str, str], b_col: dict[str, str], w: dict[str, float], skip: frozenset[str] = frozenset(),
          theta: float = THETA) -> tuple[np.ndarray, int, float]:
    """(P, species, effective sequences): the 20 x 20 joint distribution of the residues of two columns over the
    species that hold an amino acid in both and are not in `skip`, weighted by `w`, with the pseudocount `theta`."""
    f = np.zeros((20, 20)); n = 0
    for s, a in a_col.items():
        b = b_col.get(s)
        if s in skip or a not in IDX or b not in IDX or s not in w:
            continue
        f[IDX[a], IDX[b]] += w[s]; n += 1
    neff = float(f.sum())
    if neff > 0:
        f /= neff
    return (1 - theta) * f + theta / 400, n, neff


def pmi(p: np.ndarray, a: str, b: str) -> float:
    """log P(a, b) / (P(a) P(b)) under the joint distribution `p`."""
    i, j = IDX[a], IDX[b]
    return float(np.log(p[i, j] / (p[i].sum() * p[:, j].sum())))


def loss(p: np.ndarray, host_mt: str, donor_mt: str, host_nuc: str) -> float:
    """What the pair loses when the donor's mitochondrial residue meets the host's nuclear one; 0 when the donor's
    residue is the host's."""
    if host_mt == donor_mt:
        return 0.0
    return pmi(p, host_mt, host_nuc) - pmi(p, donor_mt, host_nuc)


def column(aln: dict[str, str], col: int) -> dict[str, str]:
    return {s: row[col].upper() for s, row in aln.items()}


def on_columns(human_row: str, protein: str, partial: bool) -> dict[int, str]:
    """{0-based alignment column: the protein's residue there}, placing `protein` on the human row of an alignment
    through `cybrid.on_human` ('-' a deletion, '?' not covered by a partial sequence)."""
    col = column_of_position(human_row)
    placed = cybrid.on_human(human_row.replace("-", "").upper(), protein, partial)
    return {col[p + 1]: aa for p, aa in enumerate(placed)}


def donor_on_columns(human_row: str, host: str, donor: str, donor_partial: bool) -> dict[int, str]:
    """As `on_columns`, for a donor beside its host: the donor is read through the host
    (`cybrid.donor_through_host`), so that residues the two share are not set apart by two separate alignments to the
    human protein (the review of 2026-11-22)."""
    col = column_of_position(human_row)
    placed = cybrid.donor_through_host(human_row.replace("-", "").upper(), host, donor, donor_partial)
    return {col[p + 1]: aa for p, aa in enumerate(placed)}


def complex_pairs(cx: complexes.Complex, nucaln: Path) -> tuple[pd.DataFrame, dict[str, dict[str, str]]]:
    """(the contact pairs of the complex with their controls, as `partners.pairs` builds them; the nuclear alignments)."""
    work = cx.work
    for f in ("contacts.tsv", "sites.tsv", "control_pairs.tsv", "chains.tsv"):
        if not (work / f).is_file():
            raise SystemExit(f"{work / f} is missing: `mitoport contacts --complex {cx.name}` and `mitoport burial --complex {cx.name}`")
    nuc_alns = {g: read_fasta(nucaln / f"{g}.afa") for g in cx.nuclear if (nucaln / f"{g}.afa").is_file()}
    if not nuc_alns:
        raise SystemExit(f"no partner alignment of complex {cx.name} in {nucaln}: `mitoport iface-align` on the vertebrate selections")
    sites = pd.read_csv(work / "sites.tsv", sep="\t", dtype={"chain": str})
    pr = pairs(pd.read_csv(work / "contacts.tsv", sep="\t", dtype={"chain": str, "partner_chain": str}), sites,
               nuclear_columns(cx.cif, pd.read_csv(work / "chains.tsv", sep="\t", dtype={"chain": str}), nuc_alns),
               control_of(sites, pd.read_csv(work / "control_pairs.tsv", sep="\t")))
    return pr, nuc_alns


def score_pair(pr: pd.DataFrame, mt_alns: dict, nuc_alns: dict, w: dict[str, float], host: str, donor: str,
               host_mt: dict[str, dict[int, str]], donor_mt: dict[str, dict[int, str]],
               skip: frozenset[str], theta: float = THETA) -> dict:
    """The sums of `loss` over the contact and the control pairs of one complex for one (host, donor)."""
    out = dict(coev_contact=0.0, coev_control=0.0, pairs_scored=0, pairs_changed=0, pairs_unscored=0)
    cache = {}
    for r in pr.itertuples():
        nuc_row = nuc_alns[r.nuc_gene].get(host)
        b = nuc_row[r.nuc_col].upper() if nuc_row is not None else "?"
        vals = {}
        for name, c in (("contact", r.mt_col), ("control", r.ctl_col)):
            ha, da = host_mt[r.mt_gene].get(c, "?"), donor_mt[r.mt_gene].get(c, "?")
            if ha not in IDX or da not in IDX or b not in IDX:
                vals = None; break
            key = (r.mt_gene, c, r.nuc_gene, r.nuc_col)
            if key not in cache:
                cache[key] = joint(column(mt_alns[r.mt_gene], c), column(nuc_alns[r.nuc_gene], r.nuc_col), w, skip, theta)
            p, n, _ = cache[key]
            if n < MIN_SPECIES:
                vals = None; break
            vals[name] = (loss(p, ha, da, b), ha != da)
        if vals is None:
            out["pairs_unscored"] += 1
            continue
        out["pairs_scored"] += 1; out["pairs_changed"] += vals["contact"][1]
        out["coev_contact"] += vals["contact"][0]; out["coev_control"] += vals["control"][0]
    out["coev_excess"] = out["coev_contact"] - out["coev_control"]
    out["coev_excess_mean"] = out["coev_excess"] / out["pairs_scored"] if out["pairs_scored"] else np.nan
    if not out["pairs_scored"]:
        out["coev_contact"] = out["coev_control"] = out["coev_excess"] = np.nan
    return out


def main(data: Path, gb: Path, nucaln: Path, out: Path, measurements: Path = cybrid.MEASUREMENTS) -> pd.DataFrame:
    out.mkdir(parents=True, exist_ok=True)
    meas = cybrid.read_measurements(measurements)
    records = cybrid.read_proteins(gb)
    species = sorted(set(meas.host) | set(meas.donor))
    prot = {s: cybrid.species_proteins(records, s) for s in species}
    sp = pd.read_parquet(data / "species.parquet", columns=["slug", "order"])
    order_of = dict(zip(sp.slug, sp.order))
    rows, pair_rows = [], []
    for name, cx in complexes.COMPLEXES.items():
        pr, nuc_alns = complex_pairs(cx, nucaln)
        mt_alns = {g: read_fasta(data / "alignments" / f"{g}.aa.fasta.gz") for g in cx.mt}
        with_partner = set().union(*[set(a) for a in nuc_alns.values()])
        common = sorted(s for s in set.intersection(*[set(a) for a in mt_alns.values()]) if s in with_partner)
        w = weights({s: "".join(mt_alns[g][s].upper() for g in cx.mt) for s in common})
        print(f"complex {name}: {len(pr)} contact pairs, {pr.nuc_gene.nunique()} nuclear genes with an alignment, "
              f"{len(common)} species with every mitochondrial protein and a partner, effective sequences {sum(w.values()):.1f}")
        for r in pr.itertuples():
            _, n, neff = joint(column(mt_alns[r.mt_gene], r.mt_col), column(nuc_alns[r.nuc_gene], r.nuc_col), w)
            pair_rows.append({"complex": name, "mt_gene": r.mt_gene, "mt_col": r.mt_col, "ctl_col": r.ctl_col,
                              "nuc_gene": r.nuc_gene, "nuc_col": r.nuc_col, "species": n, "effective": neff})
        placed = {s: {g: on_columns(mt_alns[g]["homo-sapiens"], prot[s][g][0], prot[s][g][1]) for g in cx.mt} for s in species}
        for h, d in sorted(set(zip(meas.host, meas.donor))):
            if h not in HOST_ORDER:
                raise SystemExit(f"no order is named for the host {h} (HOST_ORDER)")
            base = frozenset({h, d})
            no_order = base | frozenset(s for s in common if order_of.get(s) == HOST_ORDER[h])
            donor = {g: donor_on_columns(mt_alns[g]["homo-sapiens"], prot[h][g][0], prot[d][g][0], prot[d][g][1]) for g in cx.mt}
            row = {"host": h, "donor": d, "complex": name,
                   **score_pair(pr, mt_alns, nuc_alns, w, h, d, placed[h], donor, base)}
            for tag, skip, theta in (("_t50", base, THETA_SENSITIVITY), ("_noorder", no_order, THETA)):
                s = score_pair(pr, mt_alns, nuc_alns, w, h, d, placed[h], donor, skip, theta)
                row.update({f"coev_excess{tag}": s["coev_excess"], f"pairs_scored{tag}": s["pairs_scored"]})
            rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(out / "coev_scores.tsv", sep="\t", index=False, float_format="%.5f")
    pd.DataFrame(pair_rows).to_csv(out / "coev_pairs.tsv", sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.width", 220, "display.max_rows", 200):
        print(res[["host", "donor", "complex", "pairs_scored", "pairs_changed", "pairs_unscored", "coev_contact",
                   "coev_control", "coev_excess", "coev_excess_t50", "coev_excess_noorder"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    print("wrote", out / "coev_scores.tsv", out / "coev_pairs.tsv")
    return res
