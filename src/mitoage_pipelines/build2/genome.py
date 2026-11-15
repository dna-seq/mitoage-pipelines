"""GenBank records of mitogenomes -> named features, region base counts and codon counts, with the rules of
build 2 (site repository: pipeline/lib/genbank.ts, pipeline/lib/names.ts, pipeline/steps/genome.ts,
pipeline/steps/regions.ts, pipeline/lib/genetic-code.ts):

- A feature's bases are those of its annotated segments: a join() feature (the ND3 of birds and turtles,
  with an untranslated frameshift nucleotide) counts its segments only, not its start..end span.
- A region is the union of its features' positions, on the RefSeq strand, in genome order: total_mtDNA
  is the whole sequence, total_pc every CDS, total_rRNA every rRNA, total_tRNA every tRNA, D_Loop every
  D-loop feature (a misc_feature called control region only when the record has no D-loop feature), a
  single gene its first copy in genome order.
- Codons: the coding sequence in its own sense (minus-strand segments reverse-complemented), from
  codon_start, complete codons only; total_pc is the sum of the gene rows.
- The genetic code of a record is the transl_table most CDS state; without one, 2 for chordates, 9 for
  echinoderms, 5 otherwise.
- tRNA-Leu and tRNA-Ser copies whose isoacceptor the annotation does not give are placed by the
  vertebrate gene order (Leu(UUR) follows 16S rRNA, Ser(UCN) follows COX1), in chordates only.
"""
import math
import re
from dataclasses import dataclass, field

PROTEIN_GENES = ("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB")
AGGREGATE_REGIONS = ("total_mtDNA", "total_pc", "total_rRNA", "total_tRNA", "rRNA-12S", "rRNA-16S", "D_Loop")
REGIONS = AGGREGATE_REGIONS + PROTEIN_GENES
CODON_REGIONS = PROTEIN_GENES + ("total_pc",)

BASES = "UCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
STANDARD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE_DIFFS = {
    1: {}, 2: {"AGA": "*", "AGG": "*", "AUA": "M", "UGA": "W"}, 4: {"UGA": "W"},
    5: {"AGA": "S", "AGG": "S", "AUA": "M", "UGA": "W"}, 9: {"AAA": "N", "AGA": "S", "AGG": "S", "UGA": "W"},
    13: {"AGA": "G", "AGG": "G", "AUA": "M", "UGA": "W"}, 14: {"AAA": "N", "AGA": "S", "AGG": "S", "UAA": "Y", "UGA": "W"},
}
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


def genetic_code(table: int) -> dict[str, str]:
    if table not in CODE_DIFFS:
        raise ValueError(f"unsupported genetic code table {table}")
    return {c: CODE_DIFFS[table].get(c, STANDARD[i]) for i, c in enumerate(CODONS)}


def amino_acid_counts(counts: dict[str, int], table: int) -> dict[str, int]:
    code = genetic_code(table)
    out = {aa: 0 for aa in AMINO_ACIDS + ["*"]}
    for c in CODONS:
        out[code[c]] += counts.get(c, 0)
    return out


def js_round(x: float, digits: int) -> float:
    """Math.round(x * 10^d) / 10^d as the site's TypeScript computes it (halves go up, not to even)."""
    f = 10 ** digits
    return math.floor(x * f + 0.5) / f


def fallback_table(phylum: str) -> int:
    return 2 if phylum == "Chordata" else 9 if phylum == "Echinodermata" else 5


# --- GenBank flat files -----------------------------------------------------------------------
@dataclass
class Segment:
    start: int          # 1-based inclusive
    end: int
    strand: int         # 1 or -1


@dataclass
class GbFeature:
    key: str
    location: str
    segments: list[Segment]
    partial: bool
    qualifiers: dict[str, list[str]]


@dataclass
class GbRecord:
    accession_version: str
    length: int
    circular: bool
    definition: str
    organism: str
    features: list[GbFeature]
    sequence: str


def _split_top(s: str) -> list[str]:
    out, depth, cur = [], 0, ""
    for ch in s:
        if ch == "(": depth += 1
        if ch == ")": depth -= 1
        if ch == "," and depth == 0:
            out.append(cur); cur = ""
        else:
            cur += ch
    if cur: out.append(cur)
    return out


def parse_location(loc: str) -> tuple[list[Segment], bool]:
    """complement(join(1..10,20..30)) -> segments in reading order. Remote references are skipped."""
    partial = bool(re.search(r"[<>]", loc))
    clean = re.sub(r"[<>\s]", "", loc)
    segs: list[Segment] = []

    def walk(expr: str, strand: int) -> None:
        m = re.match(r"^(complement|join|order)\((.*)\)$", expr)
        if m:
            parts = _split_top(m.group(2))
            if m.group(1) == "complement":
                before = len(segs)
                for p in parts: walk(p, -strand)
                added = segs[before:]; del segs[before:]
                segs.extend(reversed(added))
            else:
                for p in parts: walk(p, strand)
            return
        if ":" in expr: return
        r = re.match(r"^(\d+)(?:\.\.|\^)(\d+)$", expr)
        if r:
            segs.append(Segment(int(r.group(1)), int(r.group(2)), strand)); return
        r = re.match(r"^(\d+)$", expr)
        if r: segs.append(Segment(int(r.group(1)), int(r.group(1)), strand))

    walk(clean, 1)
    return segs, partial


_FEAT = re.compile(r"^ {5}(\S+)\s+(\S.*)$")
_CONT = re.compile(r"^ {21}[^/]")
_QUAL = re.compile(r"^ {21}/")


def parse_genbank(text: str) -> GbRecord:
    lines = re.split(r"\r?\n", text)
    locus = next((l for l in lines if l.startswith("LOCUS")), "")
    m = re.search(r"(\d+)\s+bp", locus)
    length = int(m.group(1)) if m else 0
    circular = bool(re.search(r"\bcircular\b", locus))
    definition = accession = organism = ""
    features: list[GbFeature] = []
    seq: list[str] = []
    i, section, n = 0, "", len(lines)
    while i < n:
        line = lines[i]
        if line[:1].isupper() and "A" <= line[0] <= "Z":
            section = line.split()[0] if line.split() else ""
            if section == "DEFINITION":
                definition = line[12:].strip()
                while i + 1 < n and lines[i + 1].startswith(" " * 12):
                    i += 1; definition += " " + lines[i].strip()
            elif section == "VERSION":
                parts = line[12:].strip().split()
                accession = parts[0] if parts else ""
            i += 1
            continue
        if section == "SOURCE" and line.startswith("  ORGANISM"):
            organism = line[12:].strip()
        elif section == "FEATURES":
            fm = _FEAT.match(line)
            if fm:
                location = fm.group(2).strip()
                while i + 1 < n and _CONT.match(lines[i + 1]):
                    i += 1; location += lines[i].strip()
                quals: dict[str, list[str]] = {}
                while i + 1 < n and _QUAL.match(lines[i + 1]):
                    i += 1; q = lines[i].strip()[1:]
                    while i + 1 < n and _CONT.match(lines[i + 1]):
                        i += 1; q += " " + lines[i].strip()
                    eq = q.find("=")
                    name = q if eq == -1 else q[:eq]
                    value = "" if eq == -1 else q[eq + 1:]
                    if len(value) >= 1 and value.startswith('"') and value.endswith('"'):
                        value = value[1:-1]
                    if name == "translation":
                        value = re.sub(r"\s+", "", value)
                    quals.setdefault(name, []).append(value)
                segs, partial = parse_location(location)
                features.append(GbFeature(fm.group(1), location, segs, partial, quals))
        elif section == "ORIGIN":
            seq.append(re.sub(r"[^A-Za-z]", "", line))
        i += 1
    return GbRecord(accession, length, circular, definition, organism, features, "".join(seq).upper())


COMPLEMENT = dict(zip("ATGCRYKMSWBVDHN", "TACGYRMKSWVBHDN"))


def reverse_complement(s: str) -> str:
    """Reverse complement; a letter outside the IUPAC set becomes N (as the site's code does)."""
    return "".join(COMPLEMENT.get(c, "N") for c in reversed(s))


def extract(seq: str, segments: list[Segment]) -> str:
    """A feature's sequence in its own sense."""
    out = []
    for s in segments:
        piece = seq[s.start - 1:s.end]
        out.append(piece if s.strand == 1 else reverse_complement(piece))
    return "".join(out)


# --- feature names ----------------------------------------------------------------------------
ROMAN = {"i": "1", "ii": "2", "iii": "3", "iv": "4", "v": "5", "vi": "6"}
GENE_TABLE = {
    **dict.fromkeys(["nd1", "nad1", "nadh1", "ndh1", "mtnd1", "nd1p"], "ND1"),
    **dict.fromkeys(["nd2", "nad2", "nadh2", "ndh2", "mtnd2"], "ND2"),
    **dict.fromkeys(["nd3", "nad3", "nadh3", "ndh3", "mtnd3"], "ND3"),
    **dict.fromkeys(["nd4", "nad4", "nadh4", "ndh4", "mtnd4"], "ND4"),
    **dict.fromkeys(["nd4l", "nad4l", "nadh4l", "ndh4l", "mtnd4l"], "ND4L"),
    **dict.fromkeys(["nd5", "nad5", "nadh5", "ndh5", "mtnd5"], "ND5"),
    **dict.fromkeys(["nd6", "nad6", "nadh6", "ndh6", "mtnd6"], "ND6"),
    **dict.fromkeys(["cox1", "co1", "coi", "coxi", "cox", "mtco1", "cox1p"], "COX1"),
    **dict.fromkeys(["cox2", "co2", "coii", "coxii", "mtco2"], "COX2"),
    **dict.fromkeys(["cox3", "co3", "coiii", "coxiii", "mtco3"], "COX3"),
    **dict.fromkeys(["atp6", "atpase6", "atp6p", "a6", "mtatp6", "atpase6p"], "ATP6"),
    **dict.fromkeys(["atp8", "atpase8", "atp8p", "a8", "mtatp8"], "ATP8"),
    **dict.fromkeys(["cytb", "cob", "cyb", "mtcyb", "cytbp", "cytochromeb"], "CYTB"),
}


def protein_gene(raw: str) -> str | None:
    """A gene or product text -> one of the 13 protein genes."""
    if not raw: return None
    s = raw.lower().strip()
    s = re.sub(r"^(putative|probable)\s+", "", s, count=1)
    s = re.sub(r"\s+protein$", "", s, count=1)
    m = re.search(r"nadh\b.*?(?:subunit|chain)\s*([1-6]l?|i{1,3}|ivl?|vi?)\b", s) or re.search(r"^nadh\s*-?\s*([1-6]l?)$", s)
    if m:
        num = m.group(1); ell = num.endswith("l")
        if ell: num = num[:-1]
        num = ROMAN.get(num, num)
        if re.match(r"^[1-6]$", num):
            if ell and num == "4": return "ND4L"
            if not ell: return f"ND{num}"
    m = re.search(r"cytochrome\s*c\s*oxidase\s*(?:subunit|polypeptide|chain)?\s*(i{1,3}|[1-3])\b", s)
    if m: return f"COX{ROMAN.get(m.group(1), m.group(1))}"
    if re.search(r"^cytochrome\s*b\b", s) or re.search(r"apocytochrome\s*b", s): return "CYTB"
    m = re.search(r"atp\s*(?:synthase|ase)?\s*(?:f0|fo)?\s*(?:subunit|chain)?\s*(6|8)\b", s)
    if m: return "ATP6" if m.group(1) == "6" else "ATP8"
    return GENE_TABLE.get(re.sub(r"[\s_\-.'()]+", "", s))


def rrna(raw: str) -> str | None:
    if not raw: return None
    s = raw.lower()
    if re.search(r"\b12\s*s\b|12s|s-rrna|\bsrrna\b|rrns\b|\brns\b|small\s+subunit|\bmt-?rnr1\b|\brnr1\b|s-rna", s): return "rRNA-12S"
    if re.search(r"\b16\s*s\b|16s|l-rrna|\blrrna\b|rrnl\b|\brnl\b|large\s+subunit|\bmt-?rnr2\b|\brnr2\b|l-rna", s): return "rRNA-16S"
    return None


AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
       "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val"}
AA_NAMES = {"alanine": "Ala", "arginine": "Arg", "asparagine": "Asn", "aspartate": "Asp", "aspartic acid": "Asp",
            "cysteine": "Cys", "glutamine": "Gln", "glutamate": "Glu", "glutamic acid": "Glu", "glycine": "Gly",
            "histidine": "His", "isoleucine": "Ile", "leucine": "Leu", "lysine": "Lys", "methionine": "Met",
            "phenylalanine": "Phe", "proline": "Pro", "serine": "Ser", "threonine": "Thr", "tryptophan": "Trp",
            "tyrosine": "Tyr", "valine": "Val"}
THREE = set(AA3.values())


def trna_amino_acid(raw: str) -> str | None:
    s = raw.strip()
    m = re.search(r"(?:t-?rna|transfer\s+rna)[\s-]*([A-Za-z]{3})\b", s, flags=re.I)
    if m:
        t = m.group(1)[0].upper() + m.group(1)[1:].lower()
        if t in THREE: return t
    m = re.search(r"^(?:mt-?)?trn([A-Za-z])(?:\d|\b|-|$|\()", s, flags=re.I)
    if m: return AA3.get(m.group(1).upper())
    low = s.lower()
    for name, three in AA_NAMES.items():
        if name in low: return three
    return None


def isoacceptor_from_anticodon(aa: str, anticodon: str) -> str | None:
    ac = anticodon.upper().replace("U", "T")
    if aa == "Leu":
        if re.match(r"^.AA$", ac) or re.match(r"^[TC]AA$", ac): return "UUR"
        if re.match(r"^.AG$", ac): return "CUN"
    if aa == "Ser":
        if re.match(r"^.GA$", ac): return "UCN"
        if re.match(r"^.CT$", ac): return "AGY"
    return None


def isoacceptor_from_text(aa: str, text: str) -> str | None:
    t = text.upper().replace("T", "U")
    if aa == "Leu":
        if re.search(r"UUR|UUA|UUG|UUN", t): return "UUR"
        if re.search(r"CUN|CUA|CUU|CUG|CUC", t): return "CUN"
    if aa == "Ser":
        if re.search(r"UCN|UCA|UCU|UCG|UCC", t): return "UCN"
        if re.search(r"AGY|AGN|AGC|AGU|AGR|AGA", t): return "AGY"
    return None


def feature_name(f: GbFeature) -> tuple[str, str, str | None] | None:
    """(name, type, how a Leu/Ser isoacceptor was read) of a feature, or None when it is not one we name."""
    q = lambda k: " ".join(f.qualifiers.get(k, []))
    gene, product, note, std, tag = q("gene"), q("product"), q("note"), q("standard_name"), q("locus_tag")
    if f.key == "CDS":
        for s in (gene, product, std, note, tag):
            g = protein_gene(s)
            if g: return g, "CDS", None
        return None
    if f.key == "rRNA":
        for s in (product, gene, std, note):
            r = rrna(s)
            if r: return r, "rRNA", None
        return None
    if f.key == "tRNA":
        aa = None
        for s in (product, gene, std, note):
            if s:
                aa = trna_amino_acid(s)
                if aa: break
        anticodon = q("anticodon")
        if not aa:
            m = re.search(r"aa:([A-Za-z]{3})", anticodon)
            if m:
                t = m.group(1)[0].upper() + m.group(1)[1:].lower()
                if t in THREE: aa = t
        if not aa: return None
        if aa not in ("Leu", "Ser"): return f"tRNA-{aa}", "tRNA", None
        m = re.search(r"seq:([a-z]{3})", anticodon, flags=re.I)
        iso = isoacceptor_from_anticodon(aa, m.group(1)) if m else None
        via = "anticodon"
        if not iso:
            via = "codons"
            m = re.search(r"codons?\s+recogni[sz]ed:?\s*([ACGTUNRY]{3})", note, flags=re.I)
            iso = isoacceptor_from_text(aa, q("codon_recognized") or (m.group(1) if m else ""))
        if not iso:
            via = "name"
            for s in (product, gene):
                ac = re.search(r"\(([acgtu]{3})\)", s)
                if ac: iso = isoacceptor_from_anticodon(aa, ac.group(1))
                cf = re.search(r"\(([ACGTUNRY]{3})\)", s) or re.search(r"[-\s]([ACGTUNRY]{3})$", s)
                if not iso and cf: iso = isoacceptor_from_text(aa, cf.group(1))
                if iso: break
        if not iso:
            via = "numbering"
            m = re.search(r"trn([LS])([12])\b|trna-?(leu|ser)[-\s]?([12])\b|t-?rna-?(leu|ser)\s*\(?([12])\)?$", gene or product, flags=re.I)
            if m:
                num = m.group(2) or m.group(4) or m.group(6)
                iso = ("UUR" if num == "1" else "CUN") if aa == "Leu" else ("UCN" if num == "1" else "AGY")
        if not iso: return f"tRNA-{aa}", "tRNA", "unresolved"
        return f"tRNA-{aa}({iso})", "tRNA", via
    if f.key == "D-loop":
        return "D_Loop", "D-loop", None
    if f.key in ("misc_feature", "misc_RNA"):
        text = f"{note} {std} {product} {gene}".lower()
        if re.search(r"control\s+region|d-?loop|\bcr\b", text) and not re.search(r"pseudo|partial control|repeat", text):
            return "D_Loop", "D-loop", None
    return None


# --- a genome with its named features ---------------------------------------------------------
@dataclass
class Feature:
    name: str
    type: str
    start: int
    end: int
    strand: int
    segments: list[Segment]
    origin: str                      # "annotation", or "control-region" for a misc_feature
    codon_start: int = 1
    via: str | None = None


@dataclass
class Genome:
    record: GbRecord
    features: list[Feature]
    transl_table: int
    taxid: int | None
    unrecognised: list[str] = field(default_factory=list)
    patched: list[str] = field(default_factory=list)   # protein genes added by a hand-checked annotation patch


def feature_span(segments: list[Segment], length: int) -> tuple[int, int, int]:
    """(start, end, strand) of a feature on the circular genome; start > end when it wraps the origin."""
    strand = segments[0].strand if segments else 1
    genomic = segments if strand == 1 else list(reversed(segments))
    if any(genomic[i].start < genomic[i - 1].start for i in range(1, len(genomic))):
        return genomic[0].start, genomic[-1].end, strand
    start = min(s.start for s in segments); end = max(s.end for s in segments)
    return max(1, start), min(length, end), strand


def _forward(a: int, b: int, n: int) -> int:
    return ((b - a) % n + n) % n


def _resolve_by_position(features: list[Feature], n: int, use_gene_order: bool) -> None:
    for aa, anchor_name, near, far in (("Leu", "rRNA-16S", "UUR", "CUN"), ("Ser", "COX1", "UCN", "AGY")):
        numbered = [f for f in features if f.type == "tRNA" and f.via == "numbering" and f.name.startswith(f"tRNA-{aa}(")]
        original = [(f, f.name) for f in numbered]
        for f in numbered: f.name = f"tRNA-{aa}"
        plain = [f for f in features if f.name == f"tRNA-{aa}"]
        resolved = [f for f in features if f.name.startswith(f"tRNA-{aa}(")]

        def restore():
            for f, name in original:
                if f.name == f"tRNA-{aa}": f.name = name
        if not plain: continue
        anchor = next((f for f in features if f.name == anchor_name), None)
        if anchor is None or not use_gene_order:
            restore(); continue
        a0 = anchor.start if anchor.strand == 1 else anchor.end
        if len(plain) == 2 and not resolved:
            d = []
            for f in plain:
                pos = f.start if anchor.strand == 1 else f.end
                d.append(_forward(a0, pos, n) if anchor.strand == 1 else _forward(pos, a0, n))
            k = 0 if d[0] <= d[1] else 1
            plain[k].name = f"tRNA-{aa}({near})"; plain[1 - k].name = f"tRNA-{aa}({far})"
            plain[0].via = plain[1].via = "position"
        elif len(plain) == 1 and len(resolved) == 1:
            other = far if near in resolved[0].name else near
            plain[0].name = f"tRNA-{aa}({other})"; plain[0].via = "complement"
        else:
            restore()


# --- hand-checked annotation patches (gene_annotation.py) --------------------------------------
PRODUCT = {"ATP6": "ATP synthase F0 subunit 6", "ATP8": "ATP synthase F0 subunit 8", "COX1": "cytochrome c oxidase subunit I",
           "COX2": "cytochrome c oxidase subunit II", "COX3": "cytochrome c oxidase subunit III", "CYTB": "cytochrome b",
           "ND1": "NADH dehydrogenase subunit 1", "ND2": "NADH dehydrogenase subunit 2", "ND3": "NADH dehydrogenase subunit 3",
           "ND4": "NADH dehydrogenase subunit 4", "ND4L": "NADH dehydrogenase subunit 4L", "ND5": "NADH dehydrogenase subunit 5",
           "ND6": "NADH dehydrogenase subunit 6"}


def record_table(rec: GbRecord, fallback: int | None) -> int:
    """The genetic code the record's CDS features name most often (the lower number on a tie), else `fallback`."""
    counts: dict[int, int] = {}
    for f in rec.features:
        if f.key != "CDS": continue
        try:
            t = int(f.qualifiers.get("transl_table", ["x"])[0])
        except ValueError:
            continue
        if t > 0: counts[t] = counts.get(t, 0) + 1
    best, n = fallback, 0
    for t in sorted(counts):
        if counts[t] > n: best, n = t, counts[t]
    if best is None: raise ValueError(f"{rec.accession_version}: the record names no genetic code and none was given")
    return best


def translate_patch(cds: str, table: int, label: str) -> str:
    """The protein of a patched coding sequence, refused unless it is one open reading frame: whole codons, a stop
    codon at the end and none before it. The first codon reads as methionine (mitochondrial initiators)."""
    if len(cds) % 3: raise ValueError(f"{label}: {len(cds)} nt is not whole codons")
    code = genetic_code(table)
    aa = []
    for i in range(0, len(cds), 3):
        c = cds[i:i + 3].replace("T", "U")
        if c not in code: raise ValueError(f"{label}: codon {cds[i:i + 3]} at {i // 3 + 1} is not a plain codon")
        aa.append(code[c])
    protein = "".join(aa)
    stop = protein.find("*")
    if stop != len(protein) - 1:
        raise ValueError(f"{label}: no stop codon at the end" if stop < 0 else f"{label}: stop codon at codon {stop + 1} of {len(protein)}")
    return "M" + protein[1:-1]


def apply_gene_patches(rec: GbRecord, fallback: int | None, patches=None) -> list[str]:
    """Add the patches for this record's accession.version as CDS features; returns the genes added."""
    from .gene_annotation import ANNOTATION_PATCHES
    mine = [p for p in (ANNOTATION_PATCHES if patches is None else patches) if p.accession == rec.accession_version]
    if not mine: return []
    table = record_table(rec, fallback)
    added = []
    for p in mine:
        label = f"annotation patch {p.accession} {p.gene}"
        if any(f.key == "CDS" and (feature_name(f) or (None,))[0] == p.gene for f in rec.features):
            raise ValueError(f"{label}: the record now annotates {p.gene} itself; remove the patch")
        segs, partial = parse_location(p.location)
        if not segs or partial: raise ValueError(f"{label}: location {p.location} is empty or partial")
        if any(s.start < 1 or s.end > len(rec.sequence) or s.start > s.end for s in segs):
            raise ValueError(f"{label}: location {p.location} is outside the {len(rec.sequence)} bp record")
        translation = translate_patch(extract(rec.sequence, segs), table, label)
        rec.features.append(GbFeature("CDS", p.location, segs, False, {
            "gene": [p.gene], "product": [PRODUCT[p.gene]], "transl_table": [str(table)], "codon_start": ["1"],
            "note": [f"annotated by MitoAge, not by the record: {p.evidence}"], "translation": [translation]}))
        added.append(p.gene)
    return added


def annotate(text: str, fallback: int, vertebrate_gene_order: bool = True) -> Genome:
    rec = parse_genbank(text)
    patched = apply_gene_patches(rec, fallback)
    n = len(rec.sequence)
    features: list[Feature] = []
    control: list[Feature] = []
    unrecognised: list[str] = []
    tables: dict[int, int] = {}
    taxid = None
    for f in rec.features:
        if f.key == "source":
            tx = next((x for x in f.qualifiers.get("db_xref", []) if x.startswith("taxon:")), None)
            if tx: taxid = int(tx[6:])
            continue
        if f.key not in ("CDS", "rRNA", "tRNA", "D-loop", "misc_feature") or not f.segments:
            continue
        nm = feature_name(f)
        if nm is None:
            if f.key != "misc_feature":
                unrecognised.append(f"{f.key} {f.qualifiers.get('gene', [''])[0]} {f.qualifiers.get('product', [''])[0]}".strip())
            continue
        if f.key == "CDS":
            try:
                t = int(f.qualifiers.get("transl_table", ["x"])[0])
                if t > 0: tables[t] = tables.get(t, 0) + 1
            except ValueError:
                pass
        start, end, strand = feature_span(f.segments, n)
        try:
            cs = int(f.qualifiers.get("codon_start", ["1"])[0]) or 1
        except ValueError:
            cs = 1
        af = Feature(nm[0], nm[1], start, end, strand, f.segments,
                     "control-region" if f.key == "misc_feature" else "annotation", cs, nm[2])
        (control if af.origin == "control-region" else features).append(af)
    if not any(f.type == "D-loop" for f in features):
        features.extend(control)
    _resolve_by_position(features, n, vertebrate_gene_order)
    features.sort(key=lambda f: (f.start, f.end))
    table, best = fallback, 0
    for t, c in tables.items():
        if c > best: table, best = t, c
    return Genome(rec, features, table, taxid, unrecognised, patched)


# --- regions, codons --------------------------------------------------------------------------
def count_bases(seq: str) -> dict:
    s = seq.lower()
    a, c, g = s.count("a"), s.count("c"), s.count("g")
    t = s.count("t") + s.count("u")
    return {"a": a, "c": c, "g": g, "t": t, "other": len(s) - a - c - g - t}


def coverage(features, n: int) -> bytearray:
    mask = bytearray(n)
    for f in features:
        for s in f.segments:
            lo, hi = s.start - 1, min(s.end, n)
            if hi > lo: mask[lo:hi] = b"\x01" * (hi - lo)
    return mask


def masked(seq: str, mask: bytearray) -> str:
    return "".join(ch for ch, m in zip(seq, mask) if m)


def region_features(g: Genome, region: str):
    f = g.features
    if region == "total_mtDNA": return None
    if region == "total_pc": return [x for x in f if x.type == "CDS"]
    if region == "total_rRNA": return [x for x in f if x.type == "rRNA"]
    if region == "total_tRNA": return [x for x in f if x.type == "tRNA"]
    if region == "D_Loop": return [x for x in f if x.type == "D-loop"]
    first = next((x for x in f if x.name == region), None)
    return [first] if first else []


def compute_region(g: Genome, region: str) -> tuple[dict, int] | None:
    seq = g.record.sequence
    feats = region_features(g, region)
    if feats is None: return count_bases(seq), len(seq)
    if not feats: return None
    s = masked(seq, coverage(feats, len(seq)))
    return count_bases(s), len(s)


def non_coding(g: Genome) -> int:
    """Genome positions in no CDS, rRNA or tRNA (the D-loop counts as non-coding)."""
    n = len(g.record.sequence)
    return n - sum(coverage([f for f in g.features if f.type != "D-loop"], n))


def skew(x: int, y: int) -> float | None:
    return None if x + y == 0 else (x - y) / (x + y)


def make_stats(species_id: int, region: str, b: dict, source: str) -> dict:
    n = b["a"] + b["c"] + b["g"] + b["t"]
    gc = 0 if n == 0 else 100 * (b["g"] + b["c"]) / n
    gs, at = skew(b["g"], b["c"]), skew(b["a"], b["t"])
    return {"species_id": species_id, "region": region, "length_bp": n + b["other"], "a": b["a"], "c": b["c"], "g": b["g"],
            "t": b["t"], "other": b["other"], "gc_percent": js_round(gc, 3),
            "gc_skew": None if gs is None else js_round(gs, 5), "at_skew": None if at is None else js_round(at, 5),
            "source": source}


def count_codons(seq: str, codon_start: int = 1) -> dict[str, int]:
    counts = dict.fromkeys(CODONS, 0)
    s = seq.upper().replace("T", "U")
    for i in range(codon_start - 1, len(s) - 2, 3):
        c = s[i:i + 3]
        if c in counts: counts[c] += 1
    return counts


def gene_codons(g: Genome, gene: str) -> dict[str, int] | None:
    f = next((x for x in g.features if x.type == "CDS" and x.name == gene), None)
    if f is None: return None
    return count_codons(extract(g.record.sequence, f.segments), f.codon_start)


def usage_row(species_id: int, region: str, counts: dict, table: int, source: str) -> dict:
    return {"species_id": species_id, "region": region, "transl_table": table, "source": source,
            "counts": counts, "amino_acids": amino_acid_counts(counts, table)}


def genome_rows(species_id: int, g: Genome) -> dict:
    """Region stats, codon usage, features and non-coding length of one genome."""
    stats = []
    for region in REGIONS:
        c = compute_region(g, region)
        if c: stats.append(make_stats(species_id, region, c[0], "genbank"))
    usage, parts = [], []
    for gene in PROTEIN_GENES:
        counts = gene_codons(g, gene)
        if counts is None: continue
        parts.append(counts)
        usage.append(usage_row(species_id, gene, counts, g.transl_table, "genbank"))
    if parts:
        total = {c: sum(p[c] for p in parts) for c in CODONS}
        usage.append(usage_row(species_id, "total_pc", total, g.transl_table, "genbank"))
    feats = [{"name": f.name, "type": f.type, "start": f.start, "end": f.end, "strand": f.strand} for f in g.features]
    return {"stats": stats, "usage": usage, "features": feats, "non_coding": non_coding(g)}


def feature_keys(features: list[dict]) -> list[str]:
    """name#n, the n-th feature of that name and type in genome order: the row key of a feature."""
    seen: dict[str, int] = {}
    out = []
    for f in features:
        k = f"{f['name']}|{f['type']}"
        seen[k] = seen.get(k, 0) + 1
        out.append(f"{f['name']}#{seen[k]}")
    return out
