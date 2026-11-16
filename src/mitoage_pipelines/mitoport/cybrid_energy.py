"""H18, the interface energy: is the mito-nuclear interface of a cybrid's chimeric complex weaker than the host's own?

For one (host, donor, complex) the chimera is built on a solved structure of the host's complex: every mtDNA-encoded
chain keeps its backbone and takes the donor's residue wherever donor and host differ; the nuclear chains are the
host's as solved. The same protocol, with no substitution, gives the host's own complex as the reference. Then

    E_int = E(complex) - E(mitochondrial chains alone) - E(nuclear chains alone)

on the minimised complex (AMBER14 with the GBn2 implicit solvent, OpenMM), and the score is
`ddE = E_int(chimera) - E_int(reference)`, in kJ/mol: positive when the chimeric interface holds less.

Fixed here, before any energy was computed:
- Structures: human 5XTD (I), 5XTE (III), 5Z62 (IV), 8H9S (V); mouse 8OM1 (I), 7O3H (III), 7O37 (IV, a III2-IV
  supercomplex of which only complex IV and what touches it is kept). No mouse complex V: no murid cybrid has it measured.
- The system is the complex's mtDNA-encoded chains (every copy) and each chain with a heavy atom within 4.5 A of
  one of them. The mtDNA-encoded chains are found by sequence (at least 90% identical to the host's protein over at
  least 80% of the chain); a chain that is an mtDNA-encoded protein of another complex is dropped.
- Only protein is modelled: hemes, iron-sulfur clusters, copper, lipids, waters and ions are removed, as are
  residues without a side chain (UNK). Non-standard residues are replaced by their standard parent. Unresolved
  stretches are not built: a chain is cut at each break, and the ends are ordinary charged termini (the same in
  chimera and reference).
- Substitutions only. Where donor and host differ by a deletion or an insertion, or the residue is not resolved,
  nothing is modelled and the position is counted (`unmodelled`).
- Restraints: every heavy atom is held to its starting place (1000 kJ/mol/nm2) except the side chains of the
  residues with a heavy atom within 8 A of a substituted residue; the same residues are free in the reference.
  Minimisation: L-BFGS to 1 kJ/mol/nm or 2000 iterations. Non-bonded cutoff 1.6 nm, no periodicity.
- Noise: chimera and reference are each minimised from REPLICATES starts (free atoms displaced by Gaussian noise of
  0.005 nm, seeds 0, 1, 2). `ddE` is the difference of the means; `ddE_sd` combines the two standard deviations.
  A |ddE| below 2 `ddE_sd` is not distinguished from zero, and is still entered in the readouts as its value.
- `ddE_per_change` divides by the number of substituted residues at the interface (within 4.5 A of a nuclear chain).

`cybrid-energy` runs one (host, donor, complex) and writes energy/<host>.<donor>.<complex>.tsv; `cybrid-energy-collect`
joins them into energy_scores.tsv. Needs `uv sync --group energy`; runs on OpenMM's CPU platform.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from . import cybrid

STRUCTURES = {
    "homo-sapiens": {"I": "5XTD", "III": "5XTE", "IV": "5Z62", "V": "8H9S"},
    "mus-musculus": {"I": "8OM1", "III": "7O3H", "IV": "7O37"},
}
CONTACT_NM = 0.45
FREE_NM = 0.8
RESTRAINT_K = 1000.0
CUTOFF_NM = 1.6
MAX_ITERATIONS = 2000
TOLERANCE = 1.0
REPLICATES = 3
JITTER_NM = 0.005
MIN_IDENTITY = 0.9
MIN_COVER = 0.8
BACKBONE = {"N", "CA", "C", "O", "OXT"}
THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS",
         "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
         "Y": "TYR", "V": "VAL"}
ONE = {v: k for k, v in THREE.items()}


def triples(meas: pd.DataFrame) -> list[tuple[str, str, str]]:
    """Every (host, donor, complex) with a structure of the host's complex, in a fixed order: the array's tasks."""
    pairs = sorted(set(zip(meas.host, meas.donor)))
    return [(h, d, c) for h, d in pairs for c in STRUCTURES.get(h, {})]


def chain_sequence(chain) -> tuple[list, str]:
    res = list(chain.residues())
    return res, "".join(ONE.get(r.name, "X") for r in res)


def match_protein(seq: str, proteins: dict[str, str], al=None) -> tuple[str, dict[int, int]] | None:
    """(gene, {index in `seq`: 0-based position in the protein}) for the protein `seq` is a copy of: at least
    MIN_IDENTITY identical over aligned residues that cover at least MIN_COVER of `seq`; None when there is none."""
    al = al or cybrid.aligner(); best = None
    if len(seq) < 20:
        return None
    for g, p in proteins.items():
        a = al.align(seq, p)[0]
        pairs = [(s0 + k, p0 + k) for (s0, s1), (p0, p1) in zip(*a.aligned) for k in range(s1 - s0)]
        if not pairs or len(pairs) < MIN_COVER * len(seq):
            continue
        ident = sum(seq[i] == p[j] for i, j in pairs) / len(pairs)
        if ident >= MIN_IDENTITY and (best is None or ident > best[0]):
            best = (ident, g, dict(pairs))
    return None if best is None else (best[1], best[2])


def substitutions(struct_seq: str, to_protein: dict[int, int], host: str, donor_on_host: list[str]) -> tuple[dict[int, str], int]:
    """({index in the chain: the donor's residue}, positions not modelled) for one mtDNA-encoded chain.
    `donor_on_host[p]` is the donor's residue at host position p ('-' deleted, '?' unknown), as `cybrid.on_human`
    gives it with the host's protein in the place of the human one. A resolved position where the donor has an amino
    acid other than the host's is substituted; a deletion, an unknown residue or a residue the structure lacks is
    counted as not modelled. Insertions in the donor are not seen here and are counted by the caller."""
    subs, unmodelled = {}, 0
    resolved = set(to_protein.values())
    for i, p in to_protein.items():
        d = donor_on_host[p]
        if d == host[p]:
            continue
        if d in THREE and struct_seq[i] in THREE:
            if d != struct_seq[i]:
                subs[i] = d
        else:
            unmodelled += 1
    unmodelled += sum(1 for p, d in enumerate(donor_on_host) if p not in resolved and d != host[p])
    return subs, unmodelled


def insertions(host: str, donor: str, al=None) -> int:
    """Residues of the donor's protein that a global alignment places against no residue of the host's."""
    a = (al or cybrid.aligner()).align(host, donor)[0]
    return len(donor) - sum(o1 - o0 for _, (o0, o1) in zip(*a.aligned))


def split_at_breaks(topology, positions, max_bond_nm: float = 0.2):
    """A copy of the topology in which a chain is cut wherever the C of one residue is farther than `max_bond_nm`
    from the N of the next (or either atom is missing), each piece a chain of its own with a new id; residue names,
    ids and atom order are kept, and the standard bonds are made again. Returns (topology, {new chain index: old chain index})."""
    from openmm import app, unit
    pos = np.array(positions.value_in_unit(unit.nanometer))
    new = app.Topology(); origin = {}
    for chain in topology.chains():
        cur, prev_c = None, None
        for res in chain.residues():
            atoms = {a.name: a for a in res.atoms()}
            n = atoms.get("N")
            joined = (cur is not None and prev_c is not None and n is not None
                      and np.linalg.norm(pos[n.index] - pos[prev_c.index]) <= max_bond_nm)
            if not joined:
                cur = new.addChain(f"c{new.getNumChains()}"); origin[cur.index] = chain.index
            r = new.addResidue(res.name, cur, res.id, res.insertionCode)
            for a in res.atoms():
                new.addAtom(a.name, a.element, r)
            prev_c = atoms.get("C")
    new.createStandardBonds()
    new.createDisulfideBonds(positions)
    return new, origin


def heavy(topology) -> np.ndarray:
    return np.array([a.element is not None and a.element.symbol != "H" for a in topology.atoms()])


def near(pos: np.ndarray, a_idx: np.ndarray, b_idx: np.ndarray, cutoff: float) -> np.ndarray:
    """Which atoms of `a_idx` lie within `cutoff` of any atom of `b_idx` (a boolean array over `a_idx`)."""
    from scipy.spatial import cKDTree
    if len(a_idx) == 0 or len(b_idx) == 0:
        return np.zeros(len(a_idx), dtype=bool)
    d, _ = cKDTree(pos[b_idx]).query(pos[a_idx], distance_upper_bound=cutoff)
    return np.isfinite(d)


def prepare(cif: Path, mt_proteins: dict[str, str], other_mt: dict[str, str]):
    """The host's system from a solved entry: (topology, positions, {chain index: (gene, {residue index in the
    chain: 0-based protein position})}). Protein only, standard residues, chains cut at breaks, the complex's
    mtDNA-encoded chains and the chains that touch them; missing heavy atoms are not yet added."""
    from openmm import app, unit
    from pdbfixer import PDBFixer
    with open(cif) as f:
        fx = PDBFixer(pdbxfile=f)
    fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
    fx.removeHeterogens(keepWater=False)
    m = app.Modeller(fx.topology, fx.positions)
    m.delete([r for r in m.topology.residues() if r.name not in ONE])
    top, origin = split_at_breaks(m.topology, m.positions)
    pos = np.array(m.positions.value_in_unit(unit.nanometer))
    al = cybrid.aligner()
    # mtDNA-encoded chains are recognised on the chains as deposited, so that a short piece is not missed
    by_old = {}
    for chain in m.topology.chains():
        res, seq = chain_sequence(chain)
        hit = match_protein(seq, {**mt_proteins, **other_mt}, al)
        if hit:
            by_old[chain.index] = (hit[0], {res[i].index: p for i, p in hit[1].items()})
    found = {g for g, _ in by_old.values()}
    if not set(mt_proteins) <= found:
        raise SystemExit(f"{cif.name}: no chain of {sorted(set(mt_proteins) - found)}")
    # residue indices are the same in `top` as in the modeller's topology: residues were copied in order
    mt_chain, drop = {}, set()
    for c in top.chains():
        hit = by_old.get(origin[c.index])
        if hit is None:
            continue
        if hit[0] in mt_proteins:
            mt_chain[c.index] = (hit[0], {k: hit[1][r.index] for k, r in enumerate(c.residues()) if r.index in hit[1]})
        else:
            drop.add(c.index)
    hv = heavy(top)
    chain_of = np.array([a.residue.chain.index for a in top.atoms()])
    mt_atoms = np.flatnonzero(hv & np.isin(chain_of, list(mt_chain)))
    rest = np.flatnonzero(hv & ~np.isin(chain_of, list(mt_chain)) & ~np.isin(chain_of, list(drop)))
    touching = set(chain_of[rest[near(pos, rest, mt_atoms, CONTACT_NM)]])
    keep = set(mt_chain) | touching
    m2 = app.Modeller(top, m.positions)
    m2.delete([c for c in top.chains() if c.index not in keep])
    kept = sorted(keep)
    remap = {old: new for new, old in enumerate(kept)}
    return m2.topology, m2.positions, {remap[c]: v for c, v in mt_chain.items()}


def build(topology, positions, subs: dict[tuple[int, int], str]):
    """(topology, positions) with the substitutions {(chain index, residue index in the chain): one-letter residue}
    made, missing heavy atoms and hydrogens added; with no substitution, the reference built the same way."""
    from openmm import app
    from pdbfixer import PDBFixer
    import io
    buf = io.StringIO(); app.PDBxFile.writeFile(topology, positions, buf, keepIds=True); buf.seek(0)
    fx = PDBFixer(pdbxfile=buf)
    chains = list(fx.topology.chains())
    by_chain = {}
    for (c, k), aa in subs.items():
        res = list(chains[c].residues())[k]
        by_chain.setdefault(c, []).append(f"{res.name}-{res.id}-{THREE[aa]}")
    for c, muts in sorted(by_chain.items()):
        fx.applyMutations(muts, chains[c].id)
    fx.findMissingResidues(); fx.missingResidues = {}
    fx.findMissingAtoms(); fx.addMissingAtoms()
    m = app.Modeller(fx.topology, fx.positions)
    m.addHydrogens(forcefield(), pH=7.0)
    return m.topology, m.positions


_FF = None


def forcefield():
    global _FF
    if _FF is None:
        from openmm import app
        _FF = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    return _FF


def potential(topology, positions, restrained: np.ndarray | None = None, minimise: bool = False,
              jitter: np.ndarray | None = None, seed: int = 0, threads: int | None = None):
    """(energy in kJ/mol without the restraint term, positions). With `minimise`, the heavy atoms flagged in
    `restrained` are held to their starting places and the atoms flagged in `jitter` start displaced."""
    import openmm
    from openmm import app, unit
    system = forcefield().createSystem(topology, nonbondedMethod=app.CutoffNonPeriodic,
                                       nonbondedCutoff=CUTOFF_NM * unit.nanometer, constraints=None)
    pos = np.array(positions.value_in_unit(unit.nanometer))
    pos = pos - pos.mean(axis=0)          # single precision: keep the coordinates small
    if minimise and restrained is not None:
        f = openmm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        f.addGlobalParameter("k", RESTRAINT_K)
        for n in ("x0", "y0", "z0"):
            f.addPerParticleParameter(n)
        for i in np.flatnonzero(restrained):
            f.addParticle(int(i), pos[i].tolist())
        f.setForceGroup(1); system.addForce(f)
    props = {"Threads": str(threads)} if threads else {}
    ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("CPU"), props)
    start = pos.copy()
    if minimise and jitter is not None:
        start[jitter] += np.random.default_rng(seed).normal(0, JITTER_NM, (int(jitter.sum()), 3))
    ctx.setPositions(start)
    if minimise:
        openmm.LocalEnergyMinimizer.minimize(ctx, TOLERANCE, MAX_ITERATIONS)
    st = ctx.getState(getEnergy=True, getPositions=True, groups={0})
    return st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole), st.getPositions()


def interaction(topology, positions, mt_chains: set[int], threads: int | None = None) -> float:
    """E(complex) - E(mtDNA-encoded chains alone) - E(the other chains alone) at the given coordinates, kJ/mol."""
    from openmm import app
    total, _ = potential(topology, positions, threads=threads)
    parts = []
    for drop_mt in (False, True):
        m = app.Modeller(topology, positions)
        m.delete([c for c in topology.chains() if (c.index in mt_chains) == drop_mt])
        parts.append(potential(m.topology, m.positions, threads=threads)[0])
    return total - parts[0] - parts[1]


def masks(topology, positions, sub_residues: set[tuple[int, int]]) -> tuple[np.ndarray, np.ndarray]:
    """(restrained, free heavy atoms): free are the side-chain heavy atoms of the residues with a heavy atom within
    FREE_NM of a substituted residue (the substituted ones among them); restrained are the other heavy atoms."""
    from openmm import unit
    pos = np.array(positions.value_in_unit(unit.nanometer))
    hv = heavy(topology)
    key = np.array([(a.residue.chain.index, a.residue.index) for a in topology.atoms()])
    first = {c.index: next(c.residues()).index for c in topology.chains()}
    is_sub = np.array([(c, r - first[c]) in sub_residues for c, r in key])
    side = np.array([a.name not in BACKBONE for a in topology.atoms()])
    hidx = np.flatnonzero(hv)
    close = hidx[near(pos, hidx, np.flatnonzero(hv & is_sub), FREE_NM)]
    free_res = {tuple(k) for k in key[close]}
    in_free = np.array([tuple(k) in free_res for k in key])
    free = hv & side & in_free
    return hv & ~free, free


def run(host: str, donor: str, complex_: str, gb: Path, pdb_dir: Path, out: Path, threads: int | None = None,
        replicates: int = REPLICATES) -> dict:
    from .complexes import COMPLEXES
    pid = STRUCTURES[host][complex_]
    cif = pdb_dir / f"{pid}.cif"
    if not cif.is_file():
        raise SystemExit(f"{cif} is missing: `mitoage download pdb {pid}`")
    records = cybrid.read_proteins(gb)
    hp, dp = cybrid.species_proteins(records, host), cybrid.species_proteins(records, donor)
    genes = tuple(COMPLEXES[complex_].mt)
    mt = {g: hp[g][0] for g in genes}
    other = {g: p for g, (p, _, _) in hp.items() if g not in genes}
    top, pos, mt_chain = prepare(cif, mt, other)
    al = cybrid.aligner()
    subs, unmodelled = {}, 0
    donor_on = {g: cybrid.on_human(mt[g], dp[g][0], dp[g][1], al) for g in genes}
    for c, (g, to_protein) in sorted(mt_chain.items()):
        _, seq = chain_sequence(list(top.chains())[c])
        s, u = substitutions(seq, to_protein, mt[g], donor_on[g])
        subs.update({(c, k): aa for k, aa in s.items()}); unmodelled += u
    unmodelled += sum(insertions(mt[g], dp[g][0], al) for g in genes if not dp[g][1])
    row = {"host": host, "donor": donor, "complex": complex_, "pdb": pid, "chains": top.getNumChains(),
           "residues": top.getNumResidues(), "mt_chains": len(mt_chain), "substituted": len(subs), "unmodelled": unmodelled}
    if not subs:
        row.update(status="no substitution to model", ddE=0.0, ddE_sd=0.0, interface_substituted=0, ddE_per_change=np.nan)
        return write(row, out)
    ref_top, ref_pos = build(top, pos, {})
    chi_top, chi_pos = build(top, pos, subs)
    # interface substitutions: substituted residues with a heavy atom within CONTACT_NM of a chain that is not mtDNA-encoded
    from openmm import unit
    rp = np.array(ref_pos.value_in_unit(unit.nanometer)); hv = heavy(ref_top)
    chain_of = np.array([a.residue.chain.index for a in ref_top.atoms()])
    first = {c.index: next(c.residues()).index for c in ref_top.chains()}
    rkey = [(a.residue.chain.index, a.residue.index - first[a.residue.chain.index]) for a in ref_top.atoms()]
    sub_atoms = np.flatnonzero(hv & np.array([k in subs for k in rkey]))
    nuc_atoms = np.flatnonzero(hv & ~np.isin(chain_of, list(mt_chain)))
    at_iface = {rkey[i] for i in sub_atoms[near(rp, sub_atoms, nuc_atoms, CONTACT_NM)]}
    e = {}
    for name, (t, p) in (("reference", (ref_top, ref_pos)), ("chimera", (chi_top, chi_pos))):
        restrained, free = masks(t, p, set(subs))
        vals = []
        for seed in range(replicates):
            _, mp = potential(t, p, restrained, minimise=True, jitter=free, seed=seed, threads=threads)
            vals.append(interaction(t, mp, set(mt_chain), threads))
            print(f"{host} {donor} {complex_} {name} seed {seed}: E_int {vals[-1]:.1f} kJ/mol, "
                  f"{t.getNumAtoms()} atoms, {int(free.sum())} free heavy atoms", flush=True)
        e[name] = np.array(vals)
        row.update({f"e_int_{name}": e[name].mean(), f"e_int_{name}_sd": e[name].std(ddof=1) if replicates > 1 else np.nan,
                    f"atoms_{name}": t.getNumAtoms()})
    row["ddE"] = e["chimera"].mean() - e["reference"].mean()
    row["ddE_sd"] = float(np.hypot(row["e_int_chimera_sd"], row["e_int_reference_sd"])) if replicates > 1 else np.nan
    row["interface_substituted"] = len(at_iface)
    row["ddE_per_change"] = row["ddE"] / len(at_iface) if at_iface else np.nan
    row["status"] = "computed"
    return write(row, out)


def write(row: dict, out: Path) -> dict:
    d = out / "energy"; d.mkdir(parents=True, exist_ok=True)
    f = d / f"{row['host']}.{row['donor']}.{row['complex']}.tsv"
    pd.DataFrame([row]).to_csv(f, sep="\t", index=False, float_format="%.3f")
    print("wrote", f, row)
    return row


def collect(out: Path, measurements: Path = cybrid.MEASUREMENTS) -> pd.DataFrame:
    """energy_scores.tsv from the per-task tables; a task without a table is listed as missing, never as zero."""
    want = triples(cybrid.read_measurements(measurements))
    rows, missing = [], []
    for h, d, c in want:
        f = out / "energy" / f"{h}.{d}.{c}.tsv"
        if f.is_file():
            rows.append(pd.read_csv(f, sep="\t"))
        else:
            missing.append(f"{h} / {d} / {c}")
    if missing:
        print(f"NOT ASSESSED: {len(missing)} of {len(want)} tasks have no table: " + "; ".join(missing))
    if not rows:
        raise SystemExit("no energy table to collect")
    res = pd.concat(rows, ignore_index=True)
    res.to_csv(out / "energy_scores.tsv", sep="\t", index=False, float_format="%.3f")
    print(res.to_string(index=False)); print("wrote", out / "energy_scores.tsv")
    return res
