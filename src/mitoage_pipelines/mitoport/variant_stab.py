"""H29: does StaB-ddG tell harmful from harmless substitutions within a mitochondrial protein?

Harmful: ClinVar's records on the mitochondrial genome that are single-nucleotide, classified pathogenic or likely
pathogenic, and reviewed by its expert panel. Harmless: gnomAD v3.1's single-nucleotide variants that pass its filters
and are homoplasmic in at least AF_HOM_MIN of the people. Only missense variants of the 13 protein genes are kept, and
the amino-acid change is worked out here from the human reference (the vertebrate mitochondrial code), never taken
from either source. A variant that ClinVar calls pathogenic at any review level is left out of the harmless ones and
counted.

Each variant is scored alone on the human structures and systems of H28 (`cybrid_stab`), every copy of its chain:
`fold`, the predicted change in the stability of the whole system, and `ddG`, the change of binding of the
mtDNA-encoded chains to the rest, both positive for worse. The nuclear partner alone holds no substitution, so its
term is exactly 0 and is not computed. The weights `stabddg` with SEEDS; the two control weights with seed 0.
Variants of a task that are scored in one call share their draws of noise and decoding order, as in the model's own
script.

Baselines: `blosum`, the BLOSUM62 score of the substitution, and `n_aa`, the number of different amino acids the
build's mammals hold at the site. The rule (fixed in the entry before any score): the area under the ROC curve of
`stab_fold` for harmful against harmless is at least AUC_MIN and above each baseline's on the same variants, the
95% bootstrap interval of the difference above 0.

`download variants` fetches the two sources; `variant-stab` runs one task (a complex and a share of its variants);
`variant-stab-collect` and `variant-stab-readout` follow. Needs `uv sync --group stab` and `download stabddg`.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import cybrid
from . import cybrid_stab as cs
from .cybrid_energy import STRUCTURES

GNOMAD_URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
AF_HOM_MIN = 0.001
PANEL = "reviewed by expert panel"
CHUNKS = {"I": 12, "III": 2, "IV": 2, "V": 3}        # tasks per complex: its variants are dealt out in turn
SEEDS = {"stabddg": (0, 1, 2), "proteinmpnn": (0,), "stability_finetuned": (0,)}
BATCH = 4
AUC_MIN = 0.70
BOOTSTRAPS = 2000
MIN_PER_PROTEIN = 10
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def fetch(dest: Path, manifest: Path) -> None:
    """gnomAD's table of mitochondrial sites and ClinVar's mitochondrial records (a reduced copy: id, title, type,
    SPDI, classification, review status) into `dest`, with each file's SHA-256 in the manifest. ClinVar changes
    daily: the digest and the date say which state a result was computed on."""
    import hashlib
    import time
    from ..download import _request, append_manifest
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "gnomad.chrM.tsv").write_bytes(_request(GNOMAD_URL, accept="*/*"))
    ids = json.loads(_request(f"{EUTILS}/esearch.fcgi?db=clinvar&retmode=json&retmax=20000&term=MT%5Bchr%5D"))["esearchresult"]["idlist"]
    rows = []
    for i in range(0, len(ids), 200):
        d = json.loads(_request(f"{EUTILS}/esummary.fcgi?db=clinvar&retmode=json&id=" + ",".join(ids[i:i + 200])))["result"]
        for u in d["uids"]:
            r = d[u]; vs = (r.get("variation_set") or [{}])[0]; g = r.get("germline_classification") or {}
            rows.append({"uid": int(u), "title": r.get("title"), "type": r.get("obj_type"), "spdi": vs.get("canonical_spdi"),
                         "classification": g.get("description"), "review": g.get("review_status")})
        time.sleep(0.4)
    if len(rows) != len(ids):
        raise SystemExit(f"ClinVar listed {len(ids)} records and gave {len(rows)}: nothing written")
    (dest / "clinvar_mt.json").write_text(json.dumps(sorted(rows, key=lambda r: r["uid"]), indent=0) + "\n")
    digests = {f: hashlib.sha256((dest / f).read_bytes()).hexdigest() for f in ("gnomad.chrM.tsv", "clinvar_mt.json")}
    append_manifest(manifest, "gnomad-v3.1, clinvar-esummary", "chrM", dest, sha256=digests)
    for f, d in digests.items():
        print(f"{f}  sha256 {d}")


def reference(gb: Path) -> tuple[str, list[tuple[str, int, int, int]]]:
    """(the human mitochondrial genome, [(gene, start 0-based, end, strand)] of its 13 protein genes)."""
    from Bio import SeqIO
    acc = cybrid.GENOMES[cybrid.HUMAN]
    rec = next((r for r in SeqIO.parse(str(gb), "genbank") if r.id == acc), None)
    if rec is None:
        raise SystemExit(f"{acc} is not in {gb}: `mitoage download cybrid`")
    cds = []
    for f in rec.features:
        if f.type == "CDS":
            g = cybrid.gene_name(f.qualifiers.get("gene", [None])[0], f.qualifiers.get("product", [None])[0])
            if g:
                cds.append((g, int(f.location.start), int(f.location.end), int(f.location.strand)))
    if sorted(g for g, *_ in cds) != sorted(cybrid.ALL_GENES):
        raise SystemExit(f"{acc}: the protein genes found are {sorted(g for g, *_ in cds)}")
    return str(rec.seq).upper(), cds


def consequence(seq: str, cds: list[tuple[str, int, int, int]], pos: int, ref: str, alt: str) -> list[tuple[str, str, int, str, str]]:
    """[(kind, gene, 1-based protein position, reference amino acid, new amino acid)] for a single-nucleotide change
    at the 1-based position `pos`, one entry per protein gene that holds it (ATP8 and ATP6 overlap); kind is missense,
    synonymous or nonsense (a stop gained or lost; a codon the genome leaves incomplete counts as a stop). Empty
    outside the protein genes. Refuses a reference base that is not the genome's."""
    from Bio.Data import CodonTable
    table = CodonTable.unambiguous_dna_by_id[2].forward_table
    if seq[pos - 1] != ref:
        raise ValueError(f"position {pos} is {seq[pos - 1]} in the reference, not {ref}")
    out = []
    for g, s, e, strand in cds:
        if not s < pos <= e:
            continue
        off = pos - 1 - s if strand == 1 else e - pos
        k, i = divmod(off, 3)
        codon = seq[s + 3 * k: s + 3 * k + 3] if strand == 1 else seq[e - 3 * k - 3: e - 3 * k][::-1].translate(COMPLEMENT)
        a = alt if strand == 1 else alt.translate(COMPLEMENT)
        new = codon[:i] + a + codon[i + 1:]
        x, y = (table.get(c, "*") if len(c) == 3 else "*" for c in (codon, new))
        out.append(("synonymous" if x == y else "nonsense" if "*" in (x, y) else "missense", g, k + 1, x, y))
    return out


def is_pathogenic(classification: str | None) -> bool:
    return bool(classification) and classification.lower().startswith(("pathogenic", "likely pathogenic"))


def variants(gb: Path, clinvar: Path, gnomad: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    """(one row per missense variant and gene: id, gene, protein_pos, ref_aa, alt_aa, label, pos, ref, alt, af_hom;
    counts of what was left out and why)."""
    seq, cds = reference(gb)
    return classify(seq, cds, json.loads(clinvar.read_text()), pd.read_csv(gnomad, sep="\t"))


def classify(seq: str, cds: list[tuple[str, int, int, int]], clinvar: list[dict], g: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    """`variants` on the sources already read: ClinVar's reduced records and gnomAD's table."""
    counts = {"harmless left out: pathogenic in ClinVar": 0, "harmful without a sequence change given": 0}
    rows, any_pathogenic = [], set()
    for r in clinvar:
        if r["type"] != "single nucleotide variant" or not is_pathogenic(r["classification"]):
            continue
        sp = (r["spdi"] or "").split(":")
        if len(sp) != 4 or len(sp[2]) != 1 or len(sp[3]) != 1:
            counts["harmful without a sequence change given"] += r["review"] == PANEL
            continue
        key = (int(sp[1]) + 1, sp[2], sp[3])
        any_pathogenic.add(key)
        if r["review"] == PANEL:
            rows += [(key, "harmful", np.nan, c) for c in consequence(seq, cds, *key)]
    g = g[(g.filters == "PASS") & (g.ref.str.len() == 1) & (g.alt.str.len() == 1) & (g.AF_hom >= AF_HOM_MIN)]
    for r in g.itertuples():
        key = (int(r.position), r.ref, r.alt)
        cons = [c for c in consequence(seq, cds, *key) if c[0] == "missense"]
        if key in any_pathogenic:
            counts["harmless left out: pathogenic in ClinVar"] += bool(cons)
            continue
        rows += [(key, "harmless", float(r.AF_hom), c) for c in cons]
    out = pd.DataFrame([{"id": f"m.{k[0]}{k[1]}>{k[2]}:{c[1]}", "gene": c[1], "protein_pos": c[2], "ref_aa": c[3], "alt_aa": c[4],
                         "label": label, "pos": k[0], "ref": k[1], "alt": k[2], "af_hom": af}
                        for k, label, af, c in rows if c[0] == "missense"]).drop_duplicates("id")
    return out.sort_values(["gene", "protein_pos", "alt_aa", "pos"]).reset_index(drop=True), counts


def tasks() -> list[tuple[str, int]]:
    return [(c, k) for c in STRUCTURES[cybrid.HUMAN] for k in range(CHUNKS[c])]


def share(v: pd.DataFrame, complex_: str, chunk: int) -> pd.DataFrame:
    """The variants of the complex's genes that fall to this task: every CHUNKS[complex]-th, in the table's order."""
    from .complexes import COMPLEXES
    mine = v[v.gene.isin(list(COMPLEXES[complex_].mt))].reset_index(drop=True)
    return mine[mine.index % CHUNKS[complex_] == chunk].reset_index(drop=True)


def place(v, chains: list[dict], mt_chain: dict[int, tuple[str, dict[int, int]]]) -> tuple[dict[tuple[int, int], str], str]:
    """({(chain index, residue index): the new amino acid} in every copy of the variant's chain, status). Status is
    "scored", or why not: the position is not resolved in any copy, or the structure holds another residue there."""
    subs = {}
    for c, (g, to_protein) in sorted(mt_chain.items()):
        if g != v.gene:
            continue
        for i, p in to_protein.items():
            if p == v.protein_pos - 1:
                if chains[c]["seq"][i] != v.ref_aa:
                    return {}, f"not assessed: the structure holds {chains[c]['seq'][i]}, the reference {v.ref_aa}"
                subs[(c, i)] = v.alt_aa
    return subs, "scored" if subs else "not assessed: the position is not resolved"


def score_batch(model, whole: dict, one: dict, seqs_whole: np.ndarray, seqs_one: np.ndarray, samples: int) -> tuple[np.ndarray, np.ndarray]:
    """(fold, ddG) for each row of the two index arrays: the means over `samples` draws, positive for worse."""
    import torch
    fw, f1 = [], []
    with torch.no_grad():
        tw, t1 = (torch.from_numpy(s).to(model.device) for s in (seqs_whole, seqs_one))
        for _ in range(samples):
            fw.append(model.folding_ddG(whole, tw).reshape(-1).cpu().numpy())
            f1.append(model.folding_ddG(one, t1).reshape(-1).cpu().numpy())
    fw, f1 = -np.mean(fw, axis=0), -np.mean(f1, axis=0)
    return fw, fw - f1


def run(complex_: str, chunk: int, gb: Path, pdb_dir: Path, source: Path, weights_dir: Path, out: Path, device: str = "cpu",
        samples: int = cs.MC_SAMPLES, batch: int = BATCH, model_for=None, seeds: dict | None = None) -> pd.DataFrame:
    from .complexes import COMPLEXES
    v, _ = variants(gb, source / "clinvar_mt.json", source / "gnomad.chrM.tsv")
    v = share(v, complex_, chunk)
    pid = STRUCTURES[cybrid.HUMAN][complex_]
    hp = cybrid.species_proteins(cybrid.read_proteins(gb), cybrid.HUMAN)
    genes = tuple(COMPLEXES[complex_].mt)
    chains, mt_chain = cs.system(cs.read_chains(pdb_dir / f"{pid}.cif"), {g: hp[g][0] for g in genes},
                                 {g: p for g, (p, _, _) in hp.items() if g not in genes})
    one_idx = sorted(mt_chain)
    (whole, at), (one, at_one) = cs.structure(pid, chains, list(range(len(chains)))), cs.structure(f"{pid}_mt", chains, one_idx)
    placed = [place(r, chains, mt_chain) for r in v.itertuples()]
    todo = [k for k, (_, status) in enumerate(placed) if status == "scored"]
    rows = []
    print(f"{complex_} task {chunk}: {len(todo)} variants to score, {batch} in a call, on {device}")
    for name, (file, digest) in cs.WEIGHTS.items():
        if model_for is None:
            w = weights_dir / file
            if not w.is_file() or cs.sha256(w) != digest:
                raise SystemExit(f"{w} is missing or is not the pinned file: `mitoage download stabddg`")
            model = cs.load_model(w, device)
        else:
            model = model_for(name)
        for seed in (seeds or SEEDS)[name]:
            if model_for is None:
                import torch
                torch.manual_seed(seed)
            for b in range(0, len(todo), batch):
                part = todo[b:b + batch]
                sw = np.stack([cs.encode(whole["seq"], {at[k]: aa for k, aa in placed[i][0].items()}) for i in part])
                s1 = np.stack([cs.encode(one["seq"], {at_one[k]: aa for k, aa in placed[i][0].items()}) for i in part])
                fold, ddg = score_batch(model, whole, one, sw, s1, samples)
                for i, f, d in zip(part, fold, ddg):
                    rows.append({**v.iloc[i].to_dict(), "complex": complex_, "weights": name, "seed": seed, "copies": len(placed[i][0]),
                                 "status": "scored", "fold": float(f), "ddG": float(d)})
                print(f"{complex_} task {chunk} {name} seed {seed}: {min(b + batch, len(todo))} of {len(todo)} variants", flush=True)
    for i, (_, status) in enumerate(placed):
        if status != "scored":
            rows.append({**v.iloc[i].to_dict(), "complex": complex_, "weights": "", "seed": -1, "copies": 0, "status": status,
                         "fold": np.nan, "ddG": np.nan})
    res = pd.DataFrame(rows, columns=[*v.columns, "complex", "weights", "seed", "copies", "status", "fold", "ddG"])
    d = out / "variants"; d.mkdir(parents=True, exist_ok=True)
    f = d / f"{complex_}.{chunk}.tsv"
    res.to_csv(f, sep="\t", index=False, float_format="%.6f")
    print("wrote", f, f"({len(todo)} scored, {len(v) - len(todo)} not assessed)")
    cs.report_gpu(device)
    return res


def site_variability(data: Path, gene: str, cls: str = "Mammalia") -> dict[int, int]:
    """{1-based human protein position: how many different amino acids the build's species of the class hold there}."""
    from ..seqio import read_fasta
    from .burial import column_of_position
    aln = read_fasta(data / "alignments" / f"{gene}.aa.fasta.gz")
    sp = pd.read_parquet(data / "species.parquet", columns=["slug", "class"])
    rows = [aln[s].upper() for s in sp.slug[sp["class"] == cls] if s in aln]
    if len(rows) < 50:
        raise SystemExit(f"{gene}: {len(rows)} species of {cls} in the alignment, fewer than 50")
    return {p: len({r[c] for r in rows} & set(cs.ALPHABET[:20])) for p, c in column_of_position(aln[cybrid.HUMAN]).items()}


def collect(out: Path, data: Path, sites_dir) -> pd.DataFrame:
    """variant_scores.tsv: one row per variant with `stab_fold` and `stab_bind` (means over the seeds, and the spread
    of the first), the control weights' `fold`, the two baselines and whether the site is a mito-nuclear contact. A
    task without a table is said, never passed over; its variants are absent, not zero."""
    from Bio.Align import substitution_matrices
    from .complexes import COMPLEXES
    blosum = substitution_matrices.load("BLOSUM62")
    tabs, missing = [], []
    for c, k in tasks():
        f = out / "variants" / f"{c}.{k}.tsv"
        (tabs.append(pd.read_csv(f, sep="\t", keep_default_na=False, na_values=[""])) if f.is_file() else missing.append(f"{c} task {k}"))
    if missing:
        print(f"NOT ASSESSED: {len(missing)} of {len(tasks())} tasks have no table: " + "; ".join(missing))
    if not tabs:
        raise SystemExit("no table to collect")
    long = pd.concat(tabs, ignore_index=True)
    info = long.drop_duplicates("id").set_index("id")[["gene", "protein_pos", "ref_aa", "alt_aa", "label", "pos", "ref", "alt", "af_hom", "complex", "copies", "status"]]
    sc = long[long.status == "scored"]
    main = sc[sc.weights == "stabddg"].groupby("id")
    cols = {"stab_fold": main.fold.mean(), "stab_fold_sd": main.fold.std(ddof=1), "stab_bind": main.ddG.mean(), "seeds": main.seed.nunique()}
    for w in ("proteinmpnn", "stability_finetuned"):
        cols[f"{cs.PREFIX[w]}_fold"] = sc[sc.weights == w].groupby("id").fold.mean()
    res = info.join(pd.DataFrame(cols)).reset_index()
    res["blosum"] = [float(blosum[a][b]) for a, b in zip(res.ref_aa, res.alt_aa)]
    var = {g: site_variability(data, g) for g in sorted(set(res.gene))}
    res["n_aa"] = [var[g].get(p, np.nan) for g, p in zip(res.gene, res.protein_pos)]
    contact = set()
    for c in sorted(set(res["complex"])):
        s = pd.read_csv(sites_dir(c) / "sites.tsv", sep="\t", dtype={"chain": str})
        contact |= {(r.gene, r.human_pos) for r in s[s.group == "contact"].itertuples()}
    res["contact"] = [(g, p) in contact for g, p in zip(res.gene, res.protein_pos)]
    res = res.sort_values(["gene", "protein_pos", "alt_aa", "pos"]).reset_index(drop=True)
    res.to_csv(out / "variant_scores.tsv", sep="\t", index=False, float_format="%.5f")
    print(res.groupby(["label", "status"]).size().to_string()); print("wrote", out / "variant_scores.tsv")
    return res


def compare(harmful: np.ndarray, score: np.ndarray, base: np.ndarray, n: int = BOOTSTRAPS, seed: int = 0) -> tuple[float, float, float, float]:
    """(area of the score, area of the baseline, the 2.5% and 97.5% points of their difference over `n` resamples of
    the variants). A resample without both kinds is drawn again."""
    from scipy.stats import rankdata

    def a(s, h):          # the same area as `cybrid.auc` (ties count one half), from ranks
        n1 = int(h.sum())
        return float((rankdata(s)[h].sum() - n1 * (n1 + 1) / 2) / (n1 * (len(h) - n1)))
    if harmful.all() or not harmful.any():
        raise ValueError("the variants are of one kind only: no area, and no resample would hold both")
    rng = np.random.default_rng(seed); diffs = []
    while len(diffs) < n:
        i = rng.integers(0, len(harmful), len(harmful))
        if harmful[i].all() or not harmful[i].any():
            continue
        diffs.append(a(score[i], harmful[i]) - a(base[i], harmful[i]))
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return a(score, harmful), a(base, harmful), float(lo), float(hi)


def readout(out: Path) -> pd.DataFrame:
    """The rule on variant_scores.tsv, and beside it the same areas per protein (where it has MIN_PER_PROTEIN of each
    kind), at contacts and away from them, and for the other scores. Writes h29_summary.tsv."""
    f = out / "variant_scores.tsv"
    if not f.is_file():
        raise SystemExit(f"{f} is missing: `mitoage mitoport variant-stab-collect`")
    v = pd.read_csv(f, sep="\t")
    ok = v[(v.status == "scored") & v.stab_fold.notna() & v.n_aa.notna()].copy()
    ok["harmful"] = ok.label == "harmful"
    ok["base_blosum"], ok["base_n_aa"] = -ok.blosum, -ok.n_aa        # higher is the prediction of harm, as for the scores
    rows = [{"what": "variants", "group": "all", "score": "", "n_harmful": int((v.label == "harmful").sum()), "n_harmless": int((v.label == "harmless").sum()),
             "value": np.nan, "reads": f"{len(ok)} scored with both baselines; {len(v) - len(ok)} not assessed"}]

    def area(d, s, group, what="area"):
        rows.append({"what": what, "group": group, "score": s, "n_harmful": int(d.harmful.sum()), "n_harmless": int((~d.harmful).sum()),
                     "value": cybrid.auc(d[s], d.harmful), "reads": "area under the ROC curve, harmful against harmless"})

    for s in ("stab_fold", "stab_bind", "mpnn_fold", "stabfold_fold", "base_blosum", "base_n_aa"):
        area(ok[ok[s].notna()], s, "all")
    passed = []
    h = ok.harmful.to_numpy()
    if h.any() and not h.all():
        own = cybrid.auc(ok.stab_fold, ok.harmful)
        passed.append(bool(own >= AUC_MIN))
        for b in ("base_blosum", "base_n_aa"):
            a, ab, lo, hi = compare(h, ok.stab_fold.to_numpy(), ok[b].to_numpy())
            passed.append(bool(lo > 0))
            rows.append({"what": "difference", "group": "all", "score": f"stab_fold - {b}", "n_harmful": int(h.sum()), "n_harmless": int((~h).sum()),
                         "value": a - ab, "reads": f"95% bootstrap interval {lo:.3f} to {hi:.3f}"})
        verdict = "counts" if all(passed) else "does not count"
    else:
        verdict = "not assessed"
    rows.append({"what": "rule", "group": "all", "score": "stab_fold", "n_harmful": int(h.sum()), "n_harmless": int((~h).sum()), "value": np.nan,
                 "reads": f"{verdict}: an area of at least {AUC_MIN} and above both baselines with the interval above 0"})
    for g, d in ok.groupby("gene"):
        if d.harmful.sum() >= MIN_PER_PROTEIN and (~d.harmful).sum() >= MIN_PER_PROTEIN:
            for s in ("stab_fold", "base_blosum", "base_n_aa"):
                area(d, s, g, "area (reported, claimed as nothing)")
    for name, d in (("contact", ok[ok.contact]), ("not contact", ok[~ok.contact])):
        for s in ("stab_fold", "stab_bind"):
            area(d, s, name, "area (reported, claimed as nothing)")
    res = pd.DataFrame(rows)
    res.to_csv(out / "h29_summary.tsv", sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.width", 220, "display.max_colwidth", 110, "display.max_rows", 200):
        print(res.to_string(index=False))
    print("wrote", out / "h29_summary.tsv")
    return res
