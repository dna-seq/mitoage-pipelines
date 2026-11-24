"""H36: the nuclear OXPHOS subunits of the cybrid hosts and donors that the vertebrate alignments do not hold, taken from
NCBI RefSeq as pinned input.

The vertebrate alignments of `mitoage mitoport iface-align` (50 genes) hold none of M. spretus, M. caroli, M. pahari,
M. terricolor, Otomys irroratus and the Chinese hamster. This module chooses, for each nuclear gene of the repository's
complexes (`complexes.COMPLEXES`, 72 genes) and each species, one RefSeq protein, and writes the choice as a table
(`donor_orthologs.tsv`, next to this file) that `mitoage download donor-orthologs` fetches from by accession.version.

The choice, fixed here (docs/HYPOTHESES.md, H36):
- A species has a source when NCBI holds a RefSeq annotation of its genome: its RefSeq assembly with refseq_category
  "reference genome" (`reference_assembly`). A species with none is "no ortholog source" and nothing is substituted: a
  GenBank assembly without a RefSeq annotation (M. spretus) is none, and neither is a species with no nuclear genome at
  all (M. terricolor, Otomys irroratus).
- The ortholog of a human gene is NCBI's own: the genes of the species in the human gene's NCBI Ortholog group (NCBI
  Datasets `gene/id/<id>/orthologs`). Exactly one gene of the species in the group is the ortholog. None is "no
  ortholog". Two or more are "ambiguous (paralogs)": none is taken and the candidates are listed.
- The ortholog must be protein coding and annotated on the species' reference assembly, and have a RefSeq protein on it.
- The protein is the longest RefSeq protein of the gene on that assembly. Proteins of equal length are told apart by
  their accession, NP_ before XP_, then the number, so the choice does not depend on the order NCBI lists them in; the
  row says how many proteins the gene has and which proteins of the same length the one taken was chosen over.
"""
import hashlib
import json
import time
import urllib.parse
from pathlib import Path

import pandas as pd

from ..download import NCBI_API, NCBI_SYMBOL, EUTILS, _request, append_manifest
from .complexes import COMPLEXES

TABLE = Path(__file__).with_name("donor_orthologs.tsv")
SPECIES_TABLE = Path(__file__).with_name("donor_species.tsv")
# slug -> NCBI taxon; the hosts (mouse, Chinese hamster) and every donor of cybrid_measurements.tsv that is a rodent
SPECIES = {
    "mus-musculus": 10090, "rattus-norvegicus": 10116, "mus-spretus": 10096, "mus-caroli": 10089, "mus-pahari": 10093,
    "mus-terricolor": 254704, "otomys-irroratus": 121570, "cricetulus-griseus": 10029,
}
COLUMNS = ["species", "gene", "status", "human_gene_id", "gene_id", "symbol", "protein", "length", "sha256", "proteins", "tied_with", "tied_differ", "candidates"]
SPECIES_COLUMNS = ["species", "taxon", "assembly", "source"]
OK = "ok"
THROTTLE = 0.4          # seconds between requests: under the 3 a second NCBI allows without a key
_last = [0.0]


def genes() -> list[str]:
    """Every nuclear gene of the repository's complexes, in the order of I, III, IV, V, each once."""
    seen: dict[str, None] = {}
    for c in COMPLEXES.values():
        for g in c.nuclear:
            seen.setdefault(g)
    return list(seen)


# --- what NCBI answers: parsed, with no network -----------------------------------------------------------------------


def reference_assembly(report: dict) -> str | None:
    """The accession of the RefSeq reference genome in an NCBI Datasets `genome/taxon/<id>/dataset_report` answer made with
    `filters.assembly_source=refseq`; None when the answer holds no assembly. Refused when it holds two reference genomes
    (a choice would then be a guess), and a RefSeq assembly that is not the reference is not taken: a species whose
    only RefSeq assembly is another one (a strain's) has no pinned source."""
    ref = [r["accession"] for r in report.get("reports", [])
           if r.get("assembly_info", {}).get("refseq_category") == "reference genome" and r["accession"].startswith("GCF_")]
    if len(ref) > 1:
        raise ValueError(f"more than one RefSeq reference genome: {ref}")
    return ref[0] if ref else None


def human_gene_id(report: dict, symbol: str) -> int:
    """The id of the one human gene with exactly this symbol in a `gene/symbol/<symbol>/taxon/human` answer."""
    ids = [int(r["gene"]["gene_id"]) for r in report.get("reports", []) if r["gene"]["symbol"] == symbol]
    if len(ids) != 1:
        raise ValueError(f"{symbol}: expected one human gene with this symbol, found {len(ids)}")
    return ids[0]


def ortholog_genes(pages: list[dict]) -> dict[int, list[dict]]:
    """{taxon id: [the gene reports of that taxon in the ortholog group]} from the pages of one `orthologs` answer, each
    taxon's genes in gene id order."""
    out: dict[int, list[dict]] = {}
    for page in pages:
        for r in page.get("reports", []):
            out.setdefault(int(r["gene"]["tax_id"]), []).append(r["gene"])
    return {t: sorted(g, key=lambda x: int(x["gene_id"])) for t, g in out.items()}


def assembly_accessions(gene: dict, assembly: str) -> set[str]:
    """The genomic accessions on which a gene report says the gene is annotated on `assembly` (empty when it is not)."""
    return {loc["genomic_accession_version"] for a in gene.get("annotations", []) if a.get("assembly_accession") == assembly
            for loc in a.get("genomic_locations", [])}


def refseq_proteins(product: dict, genomic: set[str]) -> list[dict]:
    """[{accession, length, genomic}] of the RefSeq proteins of one gene's `product_report` entry whose transcript lies on
    one of the genomic accessions `genomic`, each protein once."""
    seen: dict[str, dict] = {}
    for t in product.get("transcripts", []):
        p = t.get("protein")
        if not p:
            continue
        on = {loc["genomic_accession_version"] for loc in t.get("genomic_locations", [])} & genomic
        if on:
            seen.setdefault(p["accession_version"], {"accession": p["accession_version"], "length": int(p["length"]),
                                                     "genomic": sorted(on)[0]})
    return list(seen.values())


def accession_key(acc: str) -> tuple[int, str, int, int]:
    """Sorts NP_ before XP_, then by the number and the version: NP_001234.2 -> (0, 'NP', 1234, 2)."""
    prefix, rest = acc.split("_", 1)
    num, ver = rest.split(".")
    return (0 if prefix == "NP" else 1, prefix, int(num), int(ver))


def choose_protein(proteins: list[dict]) -> tuple[dict, list[str]]:
    """(the longest protein, the accessions of the others of the same length, in the order of `accession_key`). The protein
    of equal length that comes first by `accession_key` is the one taken."""
    if not proteins:
        raise ValueError("no protein to choose from")
    top = max(p["length"] for p in proteins)
    tied = sorted((p for p in proteins if p["length"] == top), key=lambda p: accession_key(p["accession"]))
    return tied[0], [p["accession"] for p in tied[1:]]


def choose_gene(genes_of_species: list[dict], assembly: str) -> tuple[dict | None, str]:
    """(the ortholog gene, status) from the genes of one species in a human gene's ortholog group. Status: `ok`, `no
    ortholog`, `ambiguous (paralogs)`, `ortholog is not protein coding` or `ortholog not annotated on the assembly`."""
    if not genes_of_species:
        return None, "no ortholog"
    if len(genes_of_species) > 1:
        return None, "ambiguous (paralogs)"
    g = genes_of_species[0]
    if g.get("type") != "PROTEIN_CODING":
        return None, "ortholog is not protein coding"
    if not assembly_accessions(g, assembly):
        return None, "ortholog not annotated on the assembly"
    return g, OK


def sequence_digest(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()


def chosen_genes(found: dict[int, list[dict]], species: dict[str, int], assembly: dict[str, str]) -> dict[str, tuple[dict | None, str]]:
    """{species slug: (the ortholog gene or None, status)} for the species that have an assembly (`choose_gene` for each)."""
    return {slug: choose_gene(found.get(species[slug], []), acc) for slug, acc in assembly.items()}


def gene_rows(symbol: str, human_id: int, found: dict[int, list[dict]], chosen: dict, products: dict[int, dict],
              species: dict[str, int], assembly: dict[str, str]) -> list[dict]:
    """The table's rows of one human gene, one per species with an assembly: the status, the ortholog gene, and the protein
    `choose_protein` takes among the RefSeq proteins of the gene on the species' assembly (`products`: {gene id: the gene's
    product_report entry}). The sequence's SHA-256 is left empty (`fill_sequences`)."""
    rows = []
    for slug, (g, status) in chosen.items():
        row = {"species": slug, "gene": symbol, "status": status, "human_gene_id": human_id, "gene_id": "", "symbol": "", "protein": "",
               "length": "", "sha256": "", "proteins": "", "tied_with": "", "tied_differ": "", "candidates": ""}
        cand = found.get(species[slug], [])
        if len(cand) > 1:
            row["candidates"] = ";".join(f"{c['gene_id']}:{c['symbol']}:{c['type']}" for c in cand)
        if g:
            prots = refseq_proteins(products.get(int(g["gene_id"]), {}), assembly_accessions(g, assembly[slug]))
            row.update(gene_id=g["gene_id"], symbol=g["symbol"])
            if not prots:
                row["status"] = "ortholog has no RefSeq protein on the assembly"
            else:
                p, tied = choose_protein(prots)
                row.update(protein=p["accession"], length=p["length"], proteins=len(prots), tied_with=";".join(tied))
        rows.append(row)
    return rows


# --- the network ------------------------------------------------------------------------------------------------------


def _get(url: str, **kw) -> dict:
    wait = THROTTLE - (time.monotonic() - _last[0])
    if wait > 0:
        time.sleep(wait)
    out = _request(url, **kw)
    _last[0] = time.monotonic()
    return json.loads(out) if kw.get("accept", "application/json") == "application/json" else out


def _pages(url: str) -> list[dict]:
    pages, token = [], None
    while True:
        page = _get(url + (f"&page_token={urllib.parse.quote(token)}" if token else ""))
        pages.append(page)
        token = page.get("next_page_token")
        if not token:
            return pages


def select(out: Path, symbols: list[str] | None = None, species: dict[str, int] | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Choose the orthologs (the rules in the module's text) and write donor_orthologs.tsv and donor_species.tsv under `out`
    (with every answer NCBI gave kept in `out`/ncbi-answers.json). Needs the network."""
    symbols = symbols or genes()
    species = species or SPECIES
    q = urllib.parse.quote
    sp_rows, assembly = [], {}
    answers: dict = {"genome": {}, "gene": {}}
    for slug, tax in species.items():
        rep = _get(f"{NCBI_API}/genome/taxon/{tax}/dataset_report?filters.assembly_source=refseq&page_size=100")
        answers["genome"][slug] = rep
        acc = reference_assembly(rep)
        if acc:
            assembly[slug] = acc
        sp_rows.append({"species": slug, "taxon": tax, "assembly": acc or "", "source": "RefSeq" if acc else "no ortholog source: no RefSeq reference genome"})
    taxa = "&".join(f"taxon_filter={species[s]}" for s in assembly)
    rows = []
    for sym in symbols:
        ncbi = NCBI_SYMBOL.get(sym, sym)
        hid = human_gene_id(_get(f"{NCBI_API}/gene/symbol/{q(ncbi)}/taxon/human"), ncbi)
        pages = _pages(f"{NCBI_API}/gene/id/{hid}/orthologs?{taxa}&page_size=1000")
        found = ortholog_genes(pages)
        chosen = chosen_genes(found, species, assembly)
        ids = sorted({int(g["gene_id"]) for g, _ in chosen.values() if g})
        products = {}
        if ids:
            for page in _pages(f"{NCBI_API}/gene/id/{','.join(map(str, ids))}/product_report?page_size=100"):
                answers["gene"].setdefault(str(hid), []).append(page)
                for r in page.get("reports", []):
                    products[int(r["product"]["gene_id"])] = r["product"]
        new_rows = gene_rows(sym, hid, found, chosen, products, species, assembly)
        rows += new_rows
        print(f"{sym}: " + ", ".join(f"{r['species']} {r['status'] if r['status'] != OK else r['protein']}" for r in new_rows), flush=True)
    table = fill_sequences(pd.DataFrame(rows, columns=COLUMNS))
    sp = pd.DataFrame(sp_rows, columns=SPECIES_COLUMNS)
    out.mkdir(parents=True, exist_ok=True)
    (out / "ncbi-answers.json").write_text(json.dumps(answers, sort_keys=True))
    return table, sp


def fill_sequences(table: pd.DataFrame, fetch=None) -> pd.DataFrame:
    """The table with `sha256` filled from the chosen proteins' sequences, each checked against the length NCBI's product
    report gave, and `tied_differ` (how many of the proteins of equal length that were not taken have another sequence)."""
    fetch = fetch or fetch_proteins
    t = table.copy()
    ok = t.status == OK
    ties = [a for v in t.tied_with[ok] for a in v.split(";") if a]
    seqs = fetch(list(t.protein[ok]) + ties)
    for i in t.index[ok]:
        s = seqs[t.at[i, "protein"]]
        if len(s) != int(t.at[i, "length"]):
            raise ValueError(f"{t.at[i, 'protein']}: the sequence has {len(s)} residues, the product report said {t.at[i, 'length']}")
        t.at[i, "sha256"] = sequence_digest(s)
        t.at[i, "tied_differ"] = str(sum(seqs[a] != s for a in t.at[i, "tied_with"].split(";") if a))
    return t


def fetch_proteins(accessions: list[str], batch: int = 100) -> dict[str, str]:
    """{accession.version: sequence} from E-utilities efetch (db protein, FASTA), refused when an accession is missing from
    the answer or the answer holds one that was not asked for."""
    got: dict[str, str] = {}
    acc = sorted(set(accessions))
    for i in range(0, len(acc), batch):
        chunk = acc[i:i + batch]
        text = _get(f"{EUTILS}/efetch.fcgi?db=protein&id={','.join(chunk)}&rettype=fasta&retmode=text", accept="*/*").decode()
        got.update(parse_fasta(text))
    extra, missing = sorted(set(got) - set(acc)), sorted(set(acc) - set(got))
    if extra or missing:
        raise ValueError(f"efetch answered the wrong records: missing {missing[:5]}, not asked for {extra[:5]}")
    return got


def parse_fasta(text: str) -> dict[str, str]:
    """{first word of the header: the sequence, upper case, lines joined}; a header twice is refused."""
    out: dict[str, str] = {}
    name = None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in out:
                raise ValueError(f"{name} twice in the FASTA")
            out[name] = ""
        elif name is not None and line.strip():
            out[name] += line.strip().upper()
    return out


# --- the pinned input -------------------------------------------------------------------------------------------------


def read_table(path: Path = TABLE) -> pd.DataFrame:
    """The pinned choice. Every column is read as text (an empty cell is an empty string, never a missing value)."""
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_table(table: pd.DataFrame, path: Path) -> None:
    """The table in a fixed order: the species of `SPECIES`, the genes of `genes()`."""
    order_s, order_g = list(SPECIES), genes()
    t = table.assign(_s=table.species.map(order_s.index), _g=table.gene.map(order_g.index)).sort_values(["_s", "_g"])
    t.drop(columns=["_s", "_g"]).to_csv(path, sep="\t", index=False)


def download(dest: Path, manifest: Path, table_path: Path = TABLE, fetch=fetch_proteins) -> list[Path]:
    """The proteins of the pinned table, by accession.version, into `dest`/<species>.faa (headers `accession gene=<GENE>
    gene_id=<id>`), each sequence checked against the table's length and SHA-256 and the file's SHA-256 put in the
    manifest. A sequence that differs from the pin refuses the whole download: nothing is written."""
    t = read_table(table_path)
    t = t[t.status == OK]
    seqs = fetch(list(t.protein))
    for r in t.itertuples():
        s = seqs[r.protein]
        if len(s) != int(r.length) or sequence_digest(s) != r.sha256:
            raise SystemExit(f"{r.protein} ({r.species} {r.gene}) is not the pinned sequence (length {len(s)} against {r.length}, "
                             f"or another SHA-256): nothing written")
    dest.mkdir(parents=True, exist_ok=True)
    files, digests = [], {}
    for slug, g in t.groupby("species", sort=True):
        p = dest / f"{slug}.faa"
        p.write_text("".join(f">{r.protein} gene={r.gene} gene_id={r.gene_id}\n{seqs[r.protein]}\n" for r in g.itertuples()))
        digests[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
        files.append(p)
        print(f"{p.name}  {len(g)} proteins  sha256 {digests[p.name]}")
    append_manifest(manifest, "ncbi-protein-efetch", f"donor-orthologs ({len(t)} proteins of {t.species.nunique()} species)", dest,
                    sha256=digests, table_sha256=hashlib.sha256(table_path.read_bytes()).hexdigest())
    return files


class Ortholog:
    """One species' ortholog of a human gene: the protein, its accession.version and the NCBI gene id."""
    __slots__ = ("protein", "accession", "gene_id")

    def __init__(self, protein: str, accession: str, gene_id: str):
        self.protein, self.accession, self.gene_id = protein, accession, gene_id

    def __eq__(self, o):
        return isinstance(o, Ortholog) and (self.protein, self.accession, self.gene_id) == (o.protein, o.accession, o.gene_id)

    def __repr__(self):
        return f"Ortholog({self.accession}, gene {self.gene_id}, {len(self.protein)} aa)"


def read(folder: Path, table_path: Path = TABLE) -> dict[tuple[str, str], Ortholog]:
    """{(species slug, human gene symbol): Ortholog} of every `ok` row of the table, from the FASTA files of `download`,
    each sequence checked against the table's SHA-256; a file missing for a species of the table, a sequence that is
    missing or changed, or a record the table does not hold refuses. The rows that are not `ok` are not here: they are the
    gaps, and `gaps()` lists them."""
    t = read_table(table_path)
    t = t[t.status == OK]
    out: dict[tuple[str, str], Ortholog] = {}
    for slug, g in t.groupby("species", sort=True):
        f = folder / f"{slug}.faa"
        if not f.is_file():
            raise SystemExit(f"{f} is missing: `mitoage download donor-orthologs`")
        seqs = parse_fasta(f.read_text())
        if set(seqs) != set(g.protein):
            raise SystemExit(f"{f} does not hold the table's proteins of {slug} (a rerun of `mitoage download donor-orthologs`)")
        for r in g.itertuples():
            if sequence_digest(seqs[r.protein]) != r.sha256:
                raise SystemExit(f"{f}: {r.protein} is not the pinned sequence")
            out[(slug, r.gene)] = Ortholog(seqs[r.protein], r.protein, r.gene_id)
    return out


def gaps(table_path: Path = TABLE) -> pd.DataFrame:
    """The (species, gene) of the table that have no ortholog, with the reason (and the candidates of an ambiguous one)."""
    t = read_table(table_path)
    return t[t.status != OK][["species", "gene", "status", "candidates"]].reset_index(drop=True)


class Source:
    """What the pinned table says of a donor's ortholog of a gene, for the swap: `get(species, gene)` is
    (the Ortholog, None) where there is one; (None, the reason it is not assessed) where the species is covered and has none
    (no source for the species, or a gap of the table); (None, None) where the species is not covered at all, so that the
    alignments are the source. A species is never read from the table for one gene and from the alignments for another."""

    def __init__(self, orthologs: dict[tuple[str, str], Ortholog], status: dict[tuple[str, str], str], species: dict[str, str]):
        self.orthologs, self.status, self.species = orthologs, status, species

    def get(self, species: str, gene: str) -> tuple[Ortholog | None, str | None]:
        if species not in self.species:
            return None, None
        if not self.species[species]:
            return None, "not assessed: no ortholog source (no RefSeq annotation of the species' genome)"
        o = self.orthologs.get((species, gene))
        if o is not None:
            return o, None
        st = self.status.get((species, gene))
        if st is None:
            raise SystemExit(f"{species} {gene} is not in the pinned table of donor orthologs (`mitoage download donor-orthologs-select`)")
        return None, f"not assessed: no ortholog of the donor in RefSeq ({st})"


def read_source(folder: Path, table_path: Path = TABLE, species_path: Path = SPECIES_TABLE) -> Source:
    """The `Source` of the pinned table and the FASTA files of `download`."""
    sp = pd.read_csv(species_path, sep="\t", dtype=str, keep_default_na=False)
    t = read_table(table_path)
    return Source(read(folder, table_path), {(r.species, r.gene): r.status for r in t.itertuples()}, dict(zip(sp.species, sp.assembly)))
