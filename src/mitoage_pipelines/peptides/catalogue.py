"""The published peptides encoded inside mtDNA (MDPs), each with its sequence, its placement in the human genome
and the paper both were read from, and their placement in the human row of the build's alignments.

Every sequence below was copied from the primary paper named next to it (a table, a figure or the text, as
`where` says), never written from memory; each PMID was checked against its title on PubMed (2026-11-15). The
human row of the build (NC_012920.1, the revised Cambridge reference) is the check: `place` refuses a peptide
whose published sequence is not read, exactly once, from the human row under the published reading.

Readings (`reading`): "standard" is the standard genetic code (cytosolic translation, ATG start); "mito" is the
vertebrate mitochondrial code as NCBI table 2 has it (AGA and AGG are stops); "mito_agr_r" is table 2 with AGA and
AGG read as arginine, as the authors of MTALTND4 and MTALTCO1 read them.
"""
from dataclasses import dataclass
from pathlib import Path

import polars as pl
from Bio.Seq import Seq

from ..seqio import read_fasta
from . import codes


@dataclass(frozen=True)
class Peptide:
    name: str
    host: str             # region name in genome_features.parquet and in the alignment file's name
    alignment: str        # file under alignments/, without .fasta.gz
    strand: str           # "+": read on the host gene's sense strand; "-": on its complement
    reading: str          # the published reading: standard, mito, mito_agr_r
    sequence: str         # the published peptide, start residue first, no stop
    published_coordinates: str   # as the paper writes them (rCRS numbering unless the note says otherwise)
    nucleotides: str      # the published nucleotide sequence when the paper prints one, else ""
    pmid: str
    doi: str
    title: str
    where: str            # where in the paper the sequence and the placement are
    note: str = ""


COBB = dict(pmid="27070352", doi="10.18632/aging.100943",
            title="Naturally occurring mitochondrial-derived peptides are age-dependent regulators of apoptosis, "
                  "insulin sensitivity, and inflammatory markers",
            where="Figure 1A (name, length, location, sequence)")
KIENZLE = dict(pmid="37198654", doi="10.1186/s12915-023-01609-y",
               title="A small protein coded within the mitochondrial canonical gene nd4 regulates mitochondrial "
                     "bioenergetics")

PEPTIDES: tuple[Peptide, ...] = (
    Peptide("humanin", "rRNA-16S", "rRNA-16S.nt", "+", "standard", "MAPRGFSCLLLLTSEIDLPVKRRA", "2634-2707", "",
            note="Discovered by Hashimoto et al. 2001 (PMID 11371646, doi:10.1073/pnas.101133498, 'A rescue factor "
                 "abolishing neuronal cell death by a wide spectrum of familial Alzheimer's disease genes and Abeta'); "
                 "sequence and placement as Cobb et al. 2016 list them. Also read as 24 residues under mito_agr_r.",
            **COBB),
    Peptide("SHLP1", "rRNA-16S", "rRNA-16S.nt", "-", "standard", "MCHWAGGASNTGDARGDVFGKQAG", "2561-2490", "",
            note="On the complementary (light) strand; the published coordinates are 2 above rCRS, as for every SHLP.",
            **COBB),
    Peptide("SHLP2", "rRNA-16S", "rRNA-16S.nt", "-", "standard", "MGVKFFTLSTRFFPSVQRAVPLWTNS", "2170-2092", "", **COBB),
    Peptide("SHLP3", "rRNA-16S", "rRNA-16S.nt", "-", "standard", "MLGYNFSSFPCGTISIAPGFNFYRLYFIWVNGLAKVVW", "1821-1707",
            "", **COBB),
    Peptide("SHLP4", "rRNA-16S", "rRNA-16S.nt", "-", "standard", "MLEVMFLVNRRGKICRVPFTFFNLSL", "2524-2446", "",
            note="Its stop is TGA, tryptophan in the mitochondrial code.", **COBB),
    Peptide("SHLP5", "rRNA-16S", "rRNA-16S.nt", "-", "standard", "MYCSEVGFCSEVAPTEIFNAGLVV", "2856-2785", "", **COBB),
    Peptide("SHLP6", "rRNA-16S", "rRNA-16S.nt", "+", "standard", "MLDQDIPMVQPLLKVRLFND", "2992-3051", "", **COBB),
    Peptide("MOTS-c", "rRNA-12S", "rRNA-12S.nt", "+", "standard", "MRWQEMGYIFYPRKLR", "",
            "ATGAGGTGGCAAGAAATGGGCTACATTTTCTACCCCAGAAAACTACGATAG",
            pmid="25738459", doi="10.1016/j.cmet.2015.02.009",
            title="The mitochondrial-derived peptide MOTS-c promotes metabolic homeostasis and reduces obesity and "
                  "insulin resistance",
            where="Figure 1A (nucleotides; 'Cytoplasmic Translation'); the text gives no coordinates",
            note="Under the mitochondrial code the second codon, AGG, is a stop (the paper's 'M*')."),
    Peptide("MTALTND4", "ND4", "ND4.codon", "+", "mito_agr_r",
            "MRHNYNKLHLPTTNRPKIAHCMLFNQPHSPRSNSHSHPNPLKLHRRSHSHNRPRAYILITILPSKLKLRTHSQSHHNPLSRTSNSTPTNSFLMTSSKPR",
            "11,557 to 11,854", "", where="Table 1 (row MS/MS, footnote e), sequence and position",
            note="In the +3 frame of ND4 ('frameshift' in the table); the table prints AGA/AGG as arginine.", **KIENZLE),
    Peptide("MTALTCO1", "COX1", "COX1.codon", "+", "mito_agr_r",
            "ICNNLLHSNTHHNRRLWQLTSSPNNRCPRYGVSPHKQHKLLTLTSLSPTPARICYSGGRSRNRLNSLPSLSRELLPPWSLRRPNHLLLTPSRCLLYLRGHQFHHNNYQY"
            "KTPCHNPMPNAPLRLIRPNHSSPTSPISPSPSCWHHYTTNRPQPQHHLLRPRRRRRPHSMPTPILIFRSPWSLYSYPTRLRNNLPYCNLLLRKKRTIWMHRYGLSYDIN"
            "WLPRVYRVSTPYIYSRNRRRHTSMFHLRYHNHRYPHRRQSI",
            "6089 to 6866", "", pmid="41923604", doi="10.1098/rstb.2025.0078",
            title="MTALTCO1: a 259 amino-acid long mtDNA-encoded alternative protein that challenges conventional "
                  "understandings of mitochondrial genomics",
            where="identity (259 residues, +3 frame of COX1, ATT start, TAG stop, AGR read as arginine) from this "
                  "paper (preprint doi:10.1101/2024.09.22.614361, 'MTALTCO1 is a 259AA protein ...'); the residues from "
                  "Kienzle et al. 2023 (PMID 37198654) Table 1, OpenProt IP_306389, 'Within cox1 sense, frameshift "
                  "6089 to 6866'",
            note="The table prints the first residue as I (its codon ATT); here it is the start."),
    Peptide("MENTSH", "ND2", "ND2.codon", "+", "mito", "MKPNPATQNLSMLLNYPH", "start codon ATT holds m.4977 (T4977C)",
            "ATTAAACCAAACCCAGCTACGCAAAATCTTAGCATACTCCTCAATTACCCACATAGG",
            pmid="42559346", doi="10.7150/thno.134637",
            title="MENTSH: A novel mitochondrial microprotein linked to a SNP associated with type 2 diabetes",
            where="Figure 1F (codons and residues; the SNP is the start codon's second base); text: 'This start codon "
                  "(ATT)', the antibody peptide MKPNPATQNLSMLLNYPH... and T4977C (rs28357981)",
            note="Ends at AGG, a stop in the mitochondrial code; read through as arginine the frame runs on "
                 "(the paper's antibody covers MKPNPATQNLSMLLNYPHRMNNSSSTVQP)."),
)

# Published MDPs not placed, with the reason. Recorded so that a reader sees what was left out and why.
LEFT_OUT: tuple[tuple[str, str, str], ...] = (
    ("MASL", "PMID 42337615, doi:10.1186/s12979-026-00576-6, 'An immune-associated mitochondrial DNA variant with sex "
     "differences reveals a putative novel microprotein called MASL'",
     "31 residues in the D-loop; the build has no alignment of the control region, and the paper gives the sequence "
     "only in a figure without coordinates beyond the SNP (G16393A)"),
    ("SHMOOSE", "PMID 36127429, doi:10.1038/s41380-022-01769-3, 'Mitochondrial DNA variation in Alzheimer's disease "
     "reveals a unique microprotein called SHMOOSE'",
     "58 residues starting over the serine tRNA and running into ND5 (rs2853499, m.12372, is residue 47); the build "
     "aligns no tRNA, so its start and most of its length have no columns"),
    ("Gau", "Faure et al. 2011, Biology Direct (antisense strand of COX1)",
     "evidence mostly computational; antisense frames of the protein genes are not part of this test"),
    ("Kienzle 2023 Table 1 other candidates", "PMID 37198654",
     "IP_306387 (172 residues in ND2; MENTSH is its C-terminal end), IP_306398 (ND4, 35), IP_306403 (CYTB, 49) and "
     "three antisense frames found by mass spectrometry: candidates with no characterised product (the western blots "
     "show other weights or survive loss of mtDNA)"),
    ("SHLP2b", "PMID 37644144, doi:10.1038/s41598-023-41053-0, 'Evidence of natural selection in the "
     "mitochondrial-derived peptides humanin and SHLP6'", "a computational prediction, no product shown"),
)


def place(pep: Peptide, data: Path) -> dict:
    """Find `pep` in the human row of its alignment. -> the placement: index of the start codon's first base in
    the human gene's sense sequence, the rCRS coordinate of that base, the stop that follows. Refuses (ValueError)
    unless the published sequence is read exactly once."""
    aln = read_fasta(data / "alignments" / f"{pep.alignment}.fasta.gz")
    h = aln["homo-sapiens"].upper().replace("-", "")
    s = h if pep.strand == "+" else str(Seq(h).reverse_complement())
    table = codes.tables(2, pep.reading)
    n = len(pep.sequence)
    hits = []
    for i in range(0, len(s) - 3 * n + 1):
        nt = s[i:i + 3 * n]
        if not table.is_start(nt[:3]):
            continue
        if codes.translate(nt[3:], table) == pep.sequence[1:]:      # the start residue is the start codon's, M or not
            hits.append(i)
    if len(hits) != 1:
        raise ValueError(f"{pep.name}: the published sequence is read {len(hits)} times from the human {pep.host}, not once")
    i = hits[0]
    idx = i if pep.strand == "+" else len(h) - 1 - i          # the start codon's first base, in the gene's sense sequence
    gf = pl.read_parquet(data / "genome_features.parquet").filter(
        (pl.col("slug") == "homo-sapiens") & (pl.col("name") == pep.host))
    if gf.height != 1 or gf["strand"][0] != 1:
        raise ValueError(f"{pep.name}: the human {pep.host} is not annotated once on the forward strand")
    nt = s[i:i + 3 * n + 3]
    if pep.nucleotides and not nt.startswith(pep.nucleotides[:3 * n]):
        raise ValueError(f"{pep.name}: the human row differs from the published nucleotides")
    return {"peptide": pep.name, "host": pep.host, "strand": pep.strand, "reading": pep.reading, "codons": n,
            "sequence": pep.sequence, "human_index": idx, "rcrs_start": int(gf["start"][0]) + idx,
            "published_coordinates": pep.published_coordinates, "start_codon": nt[:3], "stop_after": nt[3 * n:3 * n + 3],
            "stop_after_is_stop": table.is_stop(nt[3 * n:3 * n + 3]), "pmid": pep.pmid, "doi": pep.doi,
            "title": pep.title, "where": pep.where, "note": pep.note}


def main(data: Path, out: Path) -> pl.DataFrame:
    out.mkdir(parents=True, exist_ok=True)
    rows = [place(p, data) for p in PEPTIDES]
    cat = pl.DataFrame(rows).sort("peptide")
    cat.write_csv(out / "catalogue.tsv", separator="\t")
    pl.DataFrame([{"name": n, "source": s, "reason": r} for n, s, r in LEFT_OUT]).sort("name").write_csv(
        out / "left_out.tsv", separator="\t")
    for r in rows:
        print(f"{r['peptide']}: {r['host']} {r['strand']} index {r['human_index']}, rCRS {r['rcrs_start']} "
              f"(published {r['published_coordinates'] or 'none'}), {r['codons']} codons, start {r['start_codon']}, "
              f"then {r['stop_after']}")
    return cat
