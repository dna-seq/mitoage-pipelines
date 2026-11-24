"""H36, step 2 (the maintainer's subunit swap): does giving one nuclear subunit its donor's own ortholog restore the
binding of the donor's mitochondrial proteins to it?

For one (human host, donor, complex) and one nuclear gene u of the complex, on the host's solved structure (H18's, as
H28 uses them), with StaB-ddG (H28's weights `stabddg`, defaults of its script) read as H28 reads it:

    B0(d, u)  predicted change of binding energy, partner one the complex's mtDNA-encoded chains carrying every
              substitution of the donor (the donor read through the host, as H28 builds them), partner two the chains
              of u, human. Positive: weaker binding.
    B1(d, u)  the same, with the chains of u also carrying the donor's ortholog residues.
    restoration(d, u) = B0 - B1   (positive: the donor's own u restores binding)
    null      the same B1 with the ortholog of a random other mammal of u's alignment in place of the donor's,
              drawn among the mammals whose number of differences from human at the contact residues of u is within
              25% of the donor's (at least 3 of them, else the null is not assessed), `DRAWS` draws, seed (0, draw).
              Not run when the donor's ortholog changes no placed residue: B1 is then B0, the restoration is exactly
              0, and a restoration of 0 is never among the largest.

Fixed with the code, before any score (docs/HYPOTHESES.md, H36, the changes of 2026-11-24):
- u is a nuclear gene: partner two is every chain of that gene that the system holds (every copy; a chain is kept when
  a heavy atom lies within 4.5 A of an mtDNA-encoded chain), each carrying the same ortholog residues. The system is
  not H28's: only the mtDNA-encoded chains and the chains of u are in it (three structures: both partners together,
  partner one, partner two).
- The ortholog is read from the vertebrate alignment of u (`partners/aln/<GENE>.afa`): each residue of the human chain
  is put on the human row by a global alignment (BLOSUM62, gap open -10, extend -0.5, free end gaps), and takes the
  donor's residue in that column. A gap, an unknown residue, a residue placed against none or against another residue
  than its own, or a residue the structure does not hold as a standard amino acid leaves the human residue and is
  counted as not placed. Where the donor's residue is the structure's, nothing changes.
- A contact residue of u is one with a heavy atom within 4.5 A of a heavy atom of an mtDNA-encoded chain of the system.
  "Differences at the contact residues" counts the contact columns at which a species holds an amino acid other than
  human's (a gap or an unknown is not a difference).
- The model is read with a seed set before every term (so each term depends on its inputs and the seed alone), the
  noise and decoding order drawn per sample as `StaBddG.folding_ddG` draws them, the wild-type pass of a sample shared
  by the sequences read on that draw. Rows with no substitution are not run (they are exactly 0). Seeds 0, 1, 2; 20
  samples; the score is the mean over the three seeds. B0, B1 and the null are read by one function with the same
  settings. The partner-one term does not depend on u or on the draw and is computed once per task and seed.

Where the donor's ortholog comes from RefSeq and not from the alignments (`donor_orthologs`: the hosts and donors the
alignments lack, and the 22 nuclear genes they do not hold), it is read the other way (H36, 2026-11-24): the donor's
protein is aligned to the structure's chain itself (global, BLOSUM62, gap open -10, extend -0.5, free end gaps: the
chain may be a part of the protein, and the RefSeq protein holds a presequence the structure lacks), and each residue of
the chain takes the donor residue the alignment sets against it, as `cybrid.py` reads a protein through another. A gap
against the residue, an unknown, or a residue the structure does not hold as a standard amino acid leaves the human
residue and is counted as not placed, as before. The random-mammal null keeps using the alignments: a gene with none has
no null. The alignment reading is unchanged and is the one without a pinned source (`donor_source=None`).

The mouse host (H36, 2026-11-24, `host="mus-musculus"`, `--host`): the system is the mouse complex as solved (8OM1, 7O3H, 7O37:
`cybrid_energy.STRUCTURES`), never the human one; its mtDNA-encoded chains carry the mouse's own proteins and the donor's
substitutions are read through the mouse; the nuclear chains are named by their sequence against the mouse's RefSeq proteins
(`gene_chains_by_sequence`: the structure's chains are the mouse's own, identical to them, and a chain that is another isoform is
no chain of the gene); B0 has the mouse's own subunit on the chain (the chain as solved, which `run` checks is the host's RefSeq
protein), B1 the donor's RefSeq protein placed on it, the null the alignments' random mammals measured from the mouse's row. The
donor and the host are read from the pinned RefSeq orthologs only (never one from the alignments and the other from the table).

`swap` runs one (donor, complex, gene) and writes swap/<donor>.<complex>.<gene>.tsv (one row per kind and seed) and
.cases.tsv (one row: what was placed and why a gene is not assessed); `swap-collect` joins them into the tables and
applies the candidate rule of the entry. Needs `uv sync --group stab` and the weights (`mitoage download stabddg`).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from . import cybrid, cybrid_stab as cs
from .cybrid_energy import STRUCTURES, THREE, insertions
from .complexes import COMPLEXES

HOST = cybrid.HUMAN
MOUSE = "mus-musculus"
DONORS = ("pan-troglodytes", "pan-paniscus", "gorilla-gorilla", "pongo-pygmaeus")     # with orthologs in the alignments
SEEDS = cs.SEEDS
DRAWS = 50
WINDOW = (3, 5)             # a candidate's differences k satisfy 3 n <= 4 k <= 5 n, n the donor's (within 25%)
MIN_CANDIDATES = 3
BATCH = 4                   # sequences read in one pass over a structure (the model's memory grows with it)
KINDS = ("B0", "B1", "null")
CLASS = "Mammalia"


# --- the orthologs --------------------------------------------------------------------------------------------------


def read_alignment(path: Path) -> dict[str, str]:
    """{species slug: its row, upper case} of a nuclear subunit's alignment (headers are `slug accession`)."""
    from ..seqio import read_fasta
    return {k: v.upper() for k, v in read_fasta(path).items()}


def chain_columns(chain_seq: str, human_row: str) -> list[int | None]:
    """For each residue of the structure's chain, the 0-based column of the alignment at which the human row holds that
    same residue, where a global alignment of the chain to the ungapped human row places it; None where it is placed
    against no residue or against a different one (the structure holds a sequence the alignment's human protein does
    not, a variant, another isoform or a misplaced stretch, and its residues are not read from the alignment)."""
    from Bio.Align import PairwiseAligner, substitution_matrices
    cols = [c for c, ch in enumerate(human_row) if ch != "-"]
    al = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"), open_gap_score=-10,
                         extend_gap_score=-0.5, end_gap_score=0)       # a chain may be a part of the protein, at either end
    a = al.align(chain_seq, human_row.replace("-", ""))[0]
    out: list[int | None] = [None] * len(chain_seq)
    for (s0, s1), (h0, h1) in zip(*a.aligned):
        for k in range(s1 - s0):
            if chain_seq[s0 + k] == human_row[cols[h0 + k]]:
                out[s0 + k] = cols[h0 + k]
    return out


def ortholog_subs(chain_seq: str, cols: list[int | None], row: str) -> tuple[dict[int, str], int, int]:
    """({residue index in the chain: the species' residue}, placed, not placed) for one chain. A residue is placed
    when it has a column, the structure holds a standard amino acid there and the species holds one in the column; it
    is substituted when the species' residue is another one. A residue not placed keeps the human one."""
    subs, placed, lost = {}, 0, 0
    for i, (aa, col) in enumerate(zip(chain_seq, cols)):
        if col is None or aa not in THREE or row[col] not in THREE:
            lost += 1
            continue
        placed += 1
        if row[col] != aa:
            subs[i] = row[col]
    return subs, placed, lost


def place_protein(chain_seq: str, protein: str) -> str:
    """The residue of `protein` that a global alignment (BLOSUM62, gap open -10, extend -0.5, free end gaps) sets against each
    residue of the structure's chain, as a string of the chain's length: '-' where the chain's residue is set against a gap.
    The chain is not required to be the protein's whole (an end of either may be left out) nor to be identical to it."""
    from Bio.Align import PairwiseAligner, substitution_matrices
    al = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"), open_gap_score=-10,
                         extend_gap_score=-0.5, end_gap_score=0)
    a = al.align(chain_seq, protein.upper())[0]
    out = ["-"] * len(chain_seq)
    for (s0, s1), (p0, p1) in zip(*a.aligned):
        for k in range(s1 - s0):
            out[s0 + k] = protein[p0 + k].upper()
    return "".join(out)


def protein_subs(chain_seq: str, protein: str) -> tuple[dict[int, str], int, int, str]:
    """(`ortholog_subs` of the chain read from the protein placed on it, and the placement row): ({residue index: the donor's
    residue}, placed, not placed, the row). A residue is placed when the donor's residue set against it is an amino acid
    and the structure holds a standard amino acid there."""
    row = place_protein(chain_seq, protein)
    return (*ortholog_subs(chain_seq, list(range(len(chain_seq))), row), row)


def protein_contact_differences(seqs: dict[int, str], rows: dict[int, str], ids: dict[int, list],
                                flags: dict[int, np.ndarray]) -> int:
    """How many contact positions of the gene hold a donor amino acid other than the chain's: the positions are the
    identities `ids[copy][i]` of a copy's residues (an alignment column, or the residue's number where the gene has no
    alignment) so that a position seen in several copies counts once; a residue with identity None, a gap or an unknown
    in the donor or the chain is not a difference."""
    diff = set()
    for k, f in flags.items():
        for i, hit in enumerate(f):
            if hit and ids[k][i] is not None and rows[k][i] in THREE and seqs[k][i] in THREE and rows[k][i] != seqs[k][i]:
                diff.add(ids[k][i])
    return len(diff)


def contact_flags(chains: list[dict], mt_idx: list[int], which: list[int]) -> dict[int, np.ndarray]:
    """{chain index: bool per residue} for the chains `which`: a heavy atom within CONTACT_A of a heavy atom of an
    mtDNA-encoded chain (`mt_idx`). Needs the chains' `n_heavy` (`cybrid_stab.read_chains`)."""
    from scipy.spatial import cKDTree
    tree = cKDTree(np.concatenate([chains[i]["heavy"] for i in mt_idx]))
    out = {}
    for c in which:
        d = tree.query(chains[c]["heavy"], distance_upper_bound=cs.CONTACT_A + 1e-9)[0]
        edges = np.concatenate([[0], np.cumsum(chains[c]["n_heavy"])])
        out[c] = np.array([bool((d[a:b] <= cs.CONTACT_A).any()) for a, b in zip(edges[:-1], edges[1:])], dtype=bool)
    return out


def contact_columns(cols: dict[int, list[int | None]], flags: dict[int, np.ndarray]) -> list[int]:
    """The alignment columns of the contact residues of every copy, sorted."""
    return sorted({c for k, f in flags.items() for c, hit in zip(cols[k], f) if hit and c is not None})


def contact_differences(row: str, human_row: str, columns: list[int]) -> int:
    """At how many of `columns` the row holds an amino acid other than the human row's (a gap or an unknown in either
    is not a difference)."""
    return sum(1 for c in columns if row[c] in THREE and human_row[c] in THREE and row[c] != human_row[c])


def null_candidates(rows: dict[str, str], mammals: set[str], donor: str, columns: list[int],
                    window: tuple[int, int] = WINDOW, n: int | None = None, host: str = HOST) -> tuple[list[str], int]:
    """(the sorted slugs of the other mammals of the alignment whose differences from human at the contact columns are
    within 25% of the donor's, the donor's number; differences are from the host's row, `rows[host]`). Neither the donor nor the host is a candidate. The donor's number is
    read from its row of the alignment, or is `n` where the donor's ortholog was read from its protein (`protein_subs`)."""
    human = rows[host]
    if n is None:
        n = contact_differences(rows[donor], human, columns)
    lo, hi = window
    pool = sorted(s for s in rows if s in mammals and s not in (donor, host))
    return [s for s in pool if lo * n <= 4 * contact_differences(rows[s], human, columns) <= hi * n], n


def draw_species(candidates: list[str], draw: int, seed: int = 0) -> str:
    """The species of draw number `draw`, from the generator seeded (seed, draw)."""
    return candidates[int(np.random.default_rng((seed, draw)).integers(len(candidates)))]


# --- the model ------------------------------------------------------------------------------------------------------


def ddg_rows(model, domain: dict, seqs: list[np.ndarray], seed: int, samples: int = cs.MC_SAMPLES, batch: int = BATCH) -> np.ndarray:
    """The model's folding term of `domain` for each sequence of `seqs` (index arrays), positive for worse: the mean
    over `samples` draws of what `StaBddG.folding_ddG` computes, with the generator seeded `seed` first, one decoding
    order and one backbone noise per sample shared by the wild-type pass and every sequence, and the wild-type pass
    made once per sample. A sequence equal to the domain's own is not run: it is 0."""
    import torch
    from stabddg.mpnn_utils import featurize
    wt = cs.encode(domain["seq"], {})
    out = np.zeros(len(seqs))
    todo = [k for k, s in enumerate(seqs) if not (s == wt).all()]
    if not todo:
        return out
    torch.manual_seed(seed)
    acc = np.zeros(len(todo))
    with torch.no_grad():
        X, wt_seq, _, _, chain_M, *_ = featurize([domain], model.device)
        T = torch.from_numpy(np.stack([seqs[k] for k in todo])).to(model.device)
        for _ in range(samples):
            order = model._get_decoding_order(chain_M) if model.use_antithetic_variates else None
            noise = model._get_backbone_noise(X) if model.use_antithetic_variates else None
            wt_dG = model.folding_dG(domain, wt_seq, decoding_order=order, backbone_noise=noise)
            for b in range(0, len(todo), batch):
                mut = model.folding_dG(domain, T[b:b + batch], decoding_order=order, backbone_noise=noise)
                acc[b:b + batch] += (mut - wt_dG).reshape(-1).double().cpu().numpy()
    out[todo] = -acc / samples
    return out


def binding(whole: np.ndarray, one: float, two: np.ndarray) -> np.ndarray:
    """The predicted change of binding energy, as `cybrid_stab.predict` forms it: the complex's term minus the two
    partners' terms."""
    return whole - (one + two)


# --- one task -------------------------------------------------------------------------------------------------------


def sourced_species(species_table: Path | None = None) -> dict[str, bool]:
    """{species slug: whether the pinned table of donor orthologs has a RefSeq source for it} (`donor_species.tsv`)."""
    from . import donor_orthologs as do
    t = pd.read_csv(species_table or do.SPECIES_TABLE, sep="\t", dtype=str, keep_default_na=False)
    return {r.species: bool(r.assembly) for r in t.itertuples()}


def outcome_pairs(measurements: Path, host: str) -> list[tuple[str, str]]:
    """The sorted (complex, donor) pairs of `host` with a row for a complex of the repository's (the whole-cybrid rows have no
    structure)."""
    m = cybrid.read_measurements(measurements)
    m = m[(m.host == host) & m["complex"].isin(list(COMPLEXES))]
    return sorted({(c, d) for c, d in zip(m["complex"], m.donor)})


def known_donors(pairs: list[tuple[str, str]]) -> dict[str, bool]:
    """`sourced_species()`, refusing when a donor of `pairs` is not in the table of donor species (`donor_species.tsv`)."""
    src = sourced_species()
    unknown = sorted({d for _, d in pairs if d not in src})
    if unknown:
        raise SystemExit(f"{unknown} are not in the table of donor species (`donor_species.tsv`)")
    return src


def tasks(measurements: Path = cybrid.MEASUREMENTS, host: str = HOST) -> list[tuple[str, str, str]]:
    """Every (donor, complex, nuclear gene) the swap runs, in a fixed order (the array's tasks): the donors of the host with a
    stated outcome for the complex, and every nuclear gene of the complex's structure. For the human host the donors are the
    apes that have orthologs in the alignments (`DONORS`); for another host they are the donors with a RefSeq source in the
    pinned table (`sourced_species`); a donor without one has no task and is listed by `unassessed_donors`. A donor the
    species table does not know is refused."""
    pairs = outcome_pairs(measurements, host)
    if host == HOST:
        pairs = [(c, d) for c, d in pairs if d in DONORS]
    else:
        src = known_donors(pairs)
        pairs = [(c, d) for c, d in pairs if src[d]]
    return [(d, c, g) for c, d in pairs for g in COMPLEXES[c].nuclear]


def unassessed_donors(measurements: Path = cybrid.MEASUREMENTS, host: str = MOUSE) -> pd.DataFrame:
    """One row per (host, donor, complex) of a host other than human whose donor has no RefSeq source: not run, and said."""
    pairs = outcome_pairs(measurements, host)
    src = known_donors(pairs)
    return pd.DataFrame([{"host": host, "donor": d, "complex": c, "status": "not assessed: no ortholog source (no RefSeq annotation of the species' genome)"}
                         for c, d in pairs if not src[d]], columns=["host", "donor", "complex", "status"])


def gene_chains(chains: list[dict], mt_chain: dict, chain_gene: dict[str, str]) -> dict[str, list[int]]:
    """{nuclear gene: the indices of its chains among the system's}; a nuclear chain whose gene chains.tsv does not
    name is refused."""
    out: dict[str, list[int]] = {}
    for i, c in enumerate(chains):
        if i in mt_chain:
            continue
        g = chain_gene.get(str(c["id"]))
        if g is None:
            raise SystemExit(f"chain {c['id']} has no gene in chains.tsv")
        out.setdefault(g, []).append(i)
    return out


def gene_chains_by_sequence(chains: list[dict], mt_chain: dict, proteins: dict[str, str]) -> dict[str, list[int]]:
    """{nuclear gene: the indices of its chains among the system's} for a structure of another host than human, which has no
    `chains.tsv`: a chain is the gene's when it is a copy of the gene's protein in `proteins` (`cybrid_energy.match_protein`:
    90% identical over aligned residues that cover 80% of the chain; the best gene wins). A chain of no gene (another isoform,
    a protein of no gene of the complex) is left out, never given a gene."""
    from .cybrid_energy import match_protein
    al = cybrid.aligner()
    out: dict[str, list[int]] = {}
    for i, c in enumerate(chains):
        if i in mt_chain:
            continue
        hit = match_protein(c["seq"], proteins, al)
        if hit is not None:
            out.setdefault(hit[0], []).append(i)
    return out


def case_row(donor: str, complex_: str, gene: str, status: str, host: str = HOST, **kw) -> dict:
    return {"host": host, "donor": donor, "complex": complex_, "gene": gene, "status": status, **kw}


def run(donor: str, complex_: str, gene: str, gb: Path, pdb_dir: Path, aln_dir: Path, species: Path, weights_dir: Path,
        out: Path, device: str = "cpu", seeds: tuple[int, ...] = SEEDS, samples: int = cs.MC_SAMPLES, draws: int = DRAWS,
        batch: int = BATCH, model_for=None, donor_source=None, host: str = HOST) -> tuple[pd.DataFrame, dict]:
    """(the scores, the case row) of one (donor, complex, nuclear gene); both written under `out`/swap.
    `host` is the species whose nucleus the donor's mitochondrial proteins are put against: human (the structure, the chain's own
    sequence as the baseline, the alignments or the table for the donor) or another (the mouse): then `donor_source` is required, its
    structure is that host's own in `cybrid_energy.STRUCTURES` and no other is used, and the host's ortholog and the donor's are
    both RefSeq proteins of the table.
    `model_for(weights name)` replaces the model (for tests). `donor_source` is a `donor_orthologs.Source`: a donor it
    covers is read from its RefSeq protein (`protein_subs`), and a gene it has none for is not assessed, never read from
    the alignment instead; without it, or for a donor it does not cover, the ortholog is read from the alignment."""
    if host != HOST and donor_source is None:
        raise SystemExit(f"a host other than {HOST} reads its subunits from the pinned RefSeq orthologs: give the donor orthologs")
    pid = STRUCTURES.get(host, {}).get(complex_)
    d = out / "swap"; d.mkdir(parents=True, exist_ok=True)
    stem = d / f"{donor}.{complex_}.{gene}"

    def done(case: dict, scores: pd.DataFrame | None = None):
        pd.DataFrame([case]).to_csv(f"{stem}.cases.tsv", sep="\t", index=False, float_format="%.6f")
        scores = pd.DataFrame(columns=SCORE_COLUMNS) if scores is None else scores
        scores.to_csv(f"{stem}.tsv", sep="\t", index=False, float_format="%.6f")
        print("wrote", f"{stem}.tsv", case["status"])
        return scores, case

    if pid is None:                                                 # never the human structure for another host
        return done(case_row(donor, complex_, gene, f"not assessed: no solved structure of the {host} complex {complex_}", host=host))
    cif = pdb_dir / f"{pid}.cif"
    if not cif.is_file():
        raise SystemExit(f"{cif} is missing: `mitoage download pdb {pid}`")
    aln_file = aln_dir / f"{gene}.afa"
    host_ortholog = None
    if host != HOST:
        host_ortholog, reason = donor_source.get(host, gene)
        if host_ortholog is None:
            why = "no pinned source for the host" if reason is None else reason.replace("not assessed: ", "").replace("of the donor", "of the host")
            return done(case_row(donor, complex_, gene, f"not assessed: {why}", host=host))
    ortholog, reason = donor_source.get(donor, gene) if donor_source is not None else (None, None)
    if reason is not None:
        return done(case_row(donor, complex_, gene, reason, host=host))
    if ortholog is None and host != HOST:                           # one host and one donor are never read from different sources
        return done(case_row(donor, complex_, gene, "not assessed: the pinned table does not cover the donor", host=host))
    rows = read_alignment(aln_file) if aln_file.is_file() else None
    if rows is not None and host not in rows:
        raise SystemExit(f"{aln_file} holds no row of {host}")
    if ortholog is None:                                            # the alignment is the source of the donor's ortholog
        if rows is None:
            return done(case_row(donor, complex_, gene, "not assessed: no alignment of the gene"))
        if donor not in rows:
            return done(case_row(donor, complex_, gene, "not assessed: no ortholog of the donor in the alignment"))
    records = cybrid.read_proteins(gb)
    hp, dp = (cybrid.species_proteins(records, s) for s in (host, donor))
    genes = tuple(COMPLEXES[complex_].mt)
    mt = {g: hp[g][0] for g in genes}
    chains, mt_chain = cs.system(cs.read_chains(cif), mt, {g: p for g, (p, _, _) in hp.items() if g not in genes})
    if host == HOST:
        ct = pd.read_csv(COMPLEXES[complex_].work / "chains.tsv", sep="\t", dtype={"chain": str})
        chain_gene = dict(zip(ct.chain, ct.gene))
        by_gene = gene_chains(chains, mt_chain, chain_gene)
        no_chain = "not assessed: no chain of the gene touches the mitochondrial chains"
    else:
        nuclear = {g: o.protein for g in COMPLEXES[complex_].nuclear for o in [donor_source.get(host, g)[0]] if o is not None}
        by_gene = gene_chains_by_sequence(chains, mt_chain, nuclear)
        no_chain = (f"not assessed: no chain of {pid} touching the mitochondrial chains is the host's {gene} (the structure holds none "
                    f"that is 90% identical to its RefSeq protein: it touches none, or it holds another isoform)")
    if gene not in by_gene:
        return done(case_row(donor, complex_, gene, no_chain, host=host))
    u_idx, one_idx = by_gene[gene], sorted(mt_chain)
    extra: dict = {}
    if host != HOST:                                                # B0's subunit is the chain as solved: it must be the host's own protein
        differs = sum(len(protein_subs(chains[c]["seq"], host_ortholog.protein)[0]) for c in u_idx)
        if differs:
            return done(case_row(donor, complex_, gene, f"not assessed: the chains of {pid} differ from the host's RefSeq protein at {differs} placed residues, "
                                 "so the chain as solved is not the host's own subunit", host=host))
        extra = {"host_protein": host_ortholog.accession, "host_differs": differs}
    flags = contact_flags(chains, one_idx, u_idx)
    cols = {c: chain_columns(chains[c]["seq"], rows[host]) for c in u_idx} if rows is not None else None
    columns = contact_columns(cols, flags) if cols is not None else []
    sp = pd.read_parquet(species, columns=["slug", "class"])
    pool = set(sp.slug[sp["class"] == CLASS])

    def u_subs(row: str) -> tuple[dict[tuple[int, int], str], int, int]:
        sub, placed, lost = {}, 0, 0
        for c in u_idx:
            s, p, l = ortholog_subs(chains[c]["seq"], cols[c], row)
            sub.update({(c, i): aa for i, aa in s.items()}); placed += p; lost += l
        return sub, placed, lost

    if ortholog is None:
        donor_sub, placed, lost = u_subs(rows[donor])
        candidates, n_contact_diff = null_candidates(rows, pool, donor, columns)
    else:
        donor_sub, placed, lost, placement = {}, 0, 0, {}
        for c in u_idx:
            s_, p_, l_, placement[c] = protein_subs(chains[c]["seq"], ortholog.protein)
            donor_sub.update({(c, i): aa for i, aa in s_.items()}); placed += p_; lost += l_
        ids = cols if cols is not None else {c: list(chains[c]["num"]) for c in u_idx}
        n_contact_diff = protein_contact_differences({c: chains[c]["seq"] for c in u_idx}, placement, ids, flags)
        candidates = null_candidates(rows, pool, donor, columns, n=n_contact_diff, host=host)[0] if rows is not None else []
        if cols is None:
            columns = sorted({n for c in u_idx for n, hit in zip(chains[c]["num"], flags[c]) if hit})
        extra = {**extra, "reading": "protein", "donor_protein": ortholog.accession, "donor_gene_id": ortholog.gene_id}
    al = cybrid.aligner()
    donor_on = {g: cybrid.on_human(mt[g], dp[g][0], dp[g][1], al) for g in genes}
    sets, unmodelled = cs.plan(chains, mt_chain, mt, donor_on, {g: {} for g in genes})
    unmodelled += sum(insertions(mt[g], dp[g][0], al) for g in genes if not dp[g][1])
    mt_sub = sets["all"]
    if not donor_sub:
        null = "not run: the donor's ortholog changes no placed residue, so its restoration is 0 and cannot be a candidate's"
    elif rows is None:
        null = "not assessed: no alignment of the gene (the null draws from the alignments)"
    elif len(candidates) < MIN_CANDIDATES:
        null = f"not assessed: {len(candidates)} candidates, at least {MIN_CANDIDATES} needed"
    else:
        null = "assessed"
    case = case_row(donor, complex_, gene, "scored", host=host, copies=len(u_idx), contact_residues=len(columns), placed=placed, not_placed=lost,
                    substituted=len(donor_sub), contact_differences=n_contact_diff, candidates=len(candidates), null=null,
                    mt_substituted=len(mt_sub), unmodelled=unmodelled, pdb=pid, **extra)
    # the substitution sets read: B0, B1 and the species of the null, each once however often it is drawn
    drawn = [draw_species(candidates, k) for k in range(draws)] if null == "assessed" else []
    species_sub = {s: u_subs(rows[s])[0] for s in sorted(set(drawn))}
    parts_whole = cs.structure(f"{pid}_{gene}", chains, sorted(one_idx + u_idx))
    parts_one = cs.structure(f"{pid}_mt", chains, one_idx)
    parts_two = cs.structure(f"{pid}_{gene}_nuc", chains, u_idx)
    uniq: dict[tuple, dict] = {}
    for sub in [{}, donor_sub, *species_sub.values()]:
        uniq.setdefault(tuple(sorted(sub.items())), sub)
    keys = list(uniq)
    if model_for is None:
        w = weights_dir / cs.WEIGHTS["stabddg"][0]
        if not w.is_file() or cs.sha256(w) != cs.WEIGHTS["stabddg"][1]:
            raise SystemExit(f"{w} is missing or is not the pinned file: `mitoage download stabddg`")
        model = cs.load_model(w, device)
    else:
        model = model_for("stabddg")

    def enc(part, extra: dict) -> np.ndarray:
        dom, at = part
        return cs.encode(dom["seq"], {at[k]: aa for k, aa in {**mt_sub, **extra}.items() if k in at})

    res = []
    for seed in seeds:
        fold_one = float(ddg_rows(model, parts_one[0], [enc(parts_one, {})], seed, samples, batch)[0])
        whole = ddg_rows(model, parts_whole[0], [enc(parts_whole, uniq[k]) for k in keys], seed, samples, batch)
        two = ddg_rows(model, parts_two[0], [cs.encode(parts_two[0]["seq"], {parts_two[1][c]: aa for c, aa in uniq[k].items() if c in parts_two[1]})
                                              for k in keys], seed, samples, batch)
        value = {k: (whole[j], two[j], float(binding(whole[j], fold_one, two[j]))) for j, k in enumerate(keys)}

        def row(kind: str, sub: dict, draw: int, sp: str) -> dict:
            fw, ft, dg = value[tuple(sorted(sub.items()))]
            return {**{k: case[k] for k in ("host", "donor", "complex", "gene")}, "kind": kind, "draw": draw, "species": sp, "seed": seed,
                    "substituted": len(sub), "fold": fw, "fold_one": fold_one, "fold_two": ft, "ddG": dg}
        b0, b1 = row("B0", {}, -1, host), row("B1", donor_sub, -1, donor)
        res += [b0, b1, *(row("null", species_sub[s], k, s) for k, s in enumerate(drawn))]
        print(f"{donor} {complex_} {gene} seed {seed}: B0 {b0['ddG']:.4f}, B1 {b1['ddG']:.4f}, "
              f"{len(drawn)} null draws on {len(species_sub)} species", flush=True)
    cs.report_gpu(device)
    return done(case, pd.DataFrame(res, columns=SCORE_COLUMNS))


SCORE_COLUMNS = ["host", "donor", "complex", "gene", "kind", "draw", "species", "seed", "substituted", "fold", "fold_one", "fold_two", "ddG"]



# --- the tables and the rule ----------------------------------------------------------------------------------------


def top_size(n_genes: int) -> int:
    """How many of a complex's chains count as the largest: 3, or a quarter (rounded up) when it has more than 12."""
    return 3 if n_genes <= 12 else -(-n_genes // 4)


def restoration_table(cases: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    """One row per (donor, complex, gene) that was scored: B0, B1 (means over the seeds), the restoration B0 - B1, its
    spread over the seeds, the null's mean and 95th percentile of the restoration (a draw's value is the mean over the
    seeds), the share of draws at or above the real one (the observed counted in), whether it exceeds the 95th
    percentile, and its rank among the complex's genes of the donor."""
    key = ["donor", "complex", "gene"]
    rows = []
    for r in cases[cases.status == "scored"].itertuples():
        s = scores[(scores.donor == r.donor) & (scores["complex"] == r.complex) & (scores.gene == r.gene)]
        b0 = s[s.kind == "B0"].set_index("seed").ddG
        b1 = s[s.kind == "B1"].set_index("seed").ddG
        if b0.empty or b1.empty or list(b0.index) != list(b1.index):
            raise SystemExit(f"{r.donor} {r.complex} {r.gene}: B0 and B1 are not both there for every seed")
        rest = float(b0.mean() - b1.mean())
        nul = s[s.kind == "null"]
        null = (b0.mean() - nul.groupby("draw").ddG.mean()).to_numpy() if len(nul) else np.array([])
        rows.append({"donor": r.donor, "complex": r.complex, "gene": r.gene, "B0": b0.mean(), "B1": b1.mean(), "restoration": rest,
                     "restoration_seed_sd": float((b0 - b1).std(ddof=1)) if len(b0) > 1 else np.nan, "seeds": len(b0),
                     "substituted": int(r.substituted), "contact_differences": int(r.contact_differences),
                     "candidates": int(r.candidates), "null_draws": len(null),
                     "null_mean": null.mean() if len(null) else np.nan, "null_p95": np.percentile(null, 95) if len(null) else np.nan,
                     "null_p": (1 + np.sum(null >= rest)) / (1 + len(null)) if len(null) else np.nan,
                     "exceeds_null": bool(len(null) and rest > np.percentile(null, 95))})
    res = pd.DataFrame(rows, columns=key + ["B0", "B1", "restoration", "restoration_seed_sd", "seeds", "substituted", "contact_differences",
                                            "candidates", "null_draws", "null_mean", "null_p95", "null_p", "exceeds_null"])
    if len(res):
        grp = res.groupby(["donor", "complex"])
        res["genes_ranked"] = grp.gene.transform("size")
        res["rank"] = grp.restoration.rank(method="min", ascending=False).astype(int)        # 1 + the number strictly larger
        res["in_top"] = (res.restoration > 0) & (res["rank"] <= res.genes_ranked.map(top_size))
    else:
        res["genes_ranked"], res["rank"], res["in_top"] = [], [], []
    return res.sort_values(["complex", "donor", "rank"]).reset_index(drop=True)


def candidate_table(cases: pd.DataFrame, res: pd.DataFrame) -> pd.DataFrame:
    """One row per (complex, gene) of the cases: the donors with an ortholog (the gene scored), those in which it is among
    the largest, those in which it exceeds the null's 95th percentile, those in which it is both, and whether it is a
    candidate: both in more than half of the donors that have an ortholog for it."""
    rows = []
    for (cx, g), c in cases.groupby(["complex", "gene"], sort=False):
        scored = c[c.status == "scored"]
        r = res[(res["complex"] == cx) & (res.gene == g)]
        both = int((r.in_top & r.exceeds_null).sum())
        rows.append({"complex": cx, "gene": g, "donors": len(c), "donors_with_ortholog": len(scored), "in_top": int(r.in_top.sum()),
                     "exceeds_null": int(r.exceeds_null.sum()), "both": both,
                     "candidate": bool(len(scored) and 2 * both > len(scored)),
                     "best_restoration": r.restoration.max() if len(r) else np.nan,
                     "not_assessed": "; ".join(sorted({f"{d}: {s}" for d, s in zip(c.donor, c.status) if s != "scored"}))})
    return pd.DataFrame(rows)


def activity_table(res: pd.DataFrame, cases: pd.DataFrame, meas: pd.DataFrame, host: str = HOST) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(one row per (donor, complex): the stated activity, the outcome and the restoration summed over the genes scored for
    every donor of that complex, and the Spearman correlation of activity with that sum within each complex and pooled,
    over the cases with a stated activity)."""
    from .nuclear_where import labelled
    lab = labelled(meas[meas.host == host], complexes=tuple(COMPLEXES))
    rows = []
    for cx, c in cases.groupby("complex", sort=False):
        donors = sorted(set(c.donor))
        ok = [set(c[(c.donor == d) & (c.status == "scored")].gene) for d in donors]
        common = set.intersection(*ok) if ok else set()
        for d in donors:
            r = res[(res.donor == d) & (res["complex"] == cx) & res.gene.isin(common)]
            rows.append({"donor": d, "complex": cx, "genes_summed": len(common), "restoration_sum": r.restoration.sum() if len(common) else np.nan,
                         "genes_scored": int((c[c.donor == d].status == "scored").sum())})
    sums = pd.DataFrame(rows).merge(lab[["donor", "complex", "status", "pct"]], on=["donor", "complex"], how="left")
    out = []
    for name, g in [*((cx, g) for cx, g in sums.groupby("complex", sort=False)), ("pooled", sums)]:
        g = g[g.pct.notna() & g.restoration_sum.notna()]
        rho = p = np.nan
        if len(g) >= 3 and g.pct.nunique() > 1 and g.restoration_sum.nunique() > 1:
            rho, p = stats.spearmanr(g.pct, g.restoration_sum)
        out.append({"complex": name, "n": len(g), "spearman": rho, "p": p})
    return sums, pd.DataFrame(out)


def reduced_table(res: pd.DataFrame, cases: pd.DataFrame, meas: pd.DataFrame, host: str) -> pd.DataFrame:
    """The candidate rule read on the cases the papers state as reduced (H36, mouse host): one row per (complex, gene) of the
    cases with, among the donors whose pair is reduced for that complex and that have an ortholog for the gene, how many are
    both among the largest and above the null's 95th percentile (`in_top & exceeds_null`), whether that is more than half of
    them (`candidate_reduced`: the rule as coded, on those donors alone, and never with no reduced donor), and the same count
    for the pairs stated as not reduced, which a subunit that matters for the loss should not have. Pairs the papers disagree
    on are in neither."""
    from .nuclear_where import labelled
    lab = labelled(meas[meas.host == host], complexes=tuple(COMPLEXES))[["donor", "complex", "status"]]
    r = res.merge(lab, on=["donor", "complex"], how="left")
    r["both"] = r.in_top & r.exceeds_null
    rows = []
    for (cx, g), c in cases.groupby(["complex", "gene"], sort=False):
        x = r[(r["complex"] == cx) & (r.gene == g)]
        red, nred = x[x.status == "reduced"], x[x.status == "not_reduced"]
        rows.append({"complex": cx, "gene": g, "reduced_donors": len(red), "reduced_both": int(red.both.sum()),
                     "candidate_reduced": bool(2 * int(red.both.sum()) > len(red)),
                     "not_reduced_donors": len(nred), "not_reduced_both": int(nred.both.sum()),
                     "best_reduced_restoration": red.restoration.max() if len(red) else np.nan})
    return pd.DataFrame(rows)


def collect(out: Path, measurements: Path = cybrid.MEASUREMENTS, host: str = HOST) -> pd.DataFrame:
    """swap_scores.tsv (every task's rows), swap_cases.tsv, swap_restoration.tsv, swap_candidates.tsv, swap_sums.tsv and
    swap_spearman.tsv from the tasks' tables. A task without a table is listed as missing, never as zero. For a host other than
    human, also swap_unassessed.tsv (the donors with an outcome and no RefSeq source, which have no task) and swap_reduced.tsv (the
    candidate rule on the reduced pairs alone)."""
    want = tasks(measurements, host)
    cases, scores, missing = [], [], []
    for d, c, g in want:
        f = out / "swap" / f"{d}.{c}.{g}"
        if Path(f"{f}.cases.tsv").is_file() and Path(f"{f}.tsv").is_file():
            # "null" is a kind of row here, not a missing value: pandas would read it as NaN by default
            cases.append(pd.read_csv(f"{f}.cases.tsv", sep="\t", keep_default_na=False, na_values=[""]))
            scores.append(pd.read_csv(f"{f}.tsv", sep="\t", keep_default_na=False, na_values=[""]))
        else:
            missing.append(f"{d} / {c} / {g}")
    if missing:
        print(f"MISSING: {len(missing)} of {len(want)} tasks have no table: " + "; ".join(missing))
    if not cases:
        raise SystemExit("no table to collect")
    cases, scores = pd.concat(cases, ignore_index=True), pd.concat(scores, ignore_index=True)
    res = restoration_table(cases, scores)
    cand = candidate_table(cases, res)
    sums, rho = activity_table(res, cases, cybrid.read_measurements(measurements), host)
    tables = [("scores", scores), ("cases", cases), ("restoration", res), ("candidates", cand), ("sums", sums), ("spearman", rho)]
    if host != HOST:
        tables.append(("unassessed", unassessed_donors(measurements, host)))
        tables.append(("reduced", reduced_table(res, cases, cybrid.read_measurements(measurements), host)))
    for name, t in tables:
        t.to_csv(out / f"swap_{name}.tsv", sep="\t", index=False, float_format="%.5f")
    with pd.option_context("display.width", 250, "display.max_rows", 400, "display.max_columns", 40):
        print(res.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print(cand.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print(sums.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        print(rho.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    return res
