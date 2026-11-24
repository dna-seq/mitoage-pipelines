# H36: RefSeq orthologs of the nuclear OXPHOS subunits for the cybrid hosts and donors the alignments lack

Fetched on 2026-11-24 from NCBI (Datasets v2 `gene/.../orthologs` and `product_report`, E-utilities `efetch` for the protein
sequences); code at commit `eb12408` (`mitoport/donor_orthologs.py`, `mitoage download donor-orthologs`). The files here are
copies of the pinned inputs under `src/mitoage_pipelines/mitoport/` (a test checks that the copy is the pin). Nothing of the
swap was run, and the results of its first run (`results/h36/step2/`) are not touched.

## What exists

The genes are the 72 nuclear genes of `complexes.COMPLEXES` (I 37, III 9, IV 11, V 15); the vertebrate alignments hold 50 of
them, so 22 have no alignment at all (ATP5F1A, ATP5F1B, ATP5F1C, ATP5F1D, ATP5F1E, ATP5IF1, ATP5ME, ATP5MG, ATP5PO, NDUFA2,
NDUFA5, NDUFA7, NDUFC1, NDUFS1, NDUFS4, NDUFS6, NDUFV1, NDUFV2, NDUFV3, UQCR10, UQCRC2, UQCRH).

**The 22 genes without an alignment are exactly the nuclear genes of the four structures that touch none of the mtDNA-encoded
chains** (complex I: NDUFA2, NDUFA5, NDUFA7, NDUFC1, NDUFS1, NDUFS4, NDUFS6, NDUFV1, NDUFV2, NDUFV3; III: UQCR10, UQCRC2, UQCRH;
V: ATP5F1A to ATP5F1E, ATP5IF1, ATP5ME, ATP5MG, ATP5PO; none of IV). The 50 genes with an alignment are the 27, 6, 11 and 6
that touch them (4.5 A, the swap's rule). The swap's unit is a gene whose chain touches an mtDNA-encoded chain, so for the
swap as coded the 22 extra proteins are fetched but are never a task: they end at "no chain of the gene touches the mitochondrial
chains". They are in the table because the brief asked for every nuclear gene of the complexes.

| species | role in the cybrids | RefSeq reference genome | orthologs of the 72 genes | proteins with more than one isoform (of the 72) |
| --- | --- | --- | --- | --- |
| `mus-musculus` | host | GCF_000001635.27 (GRCm39) | 72 | 33 |
| `rattus-norvegicus` | donor of the mouse host (added: it has the most stated outcomes) | GCF_036323735.1 (GRCr8) | 72 | 21 |
| `mus-caroli` | donor of the mouse host | GCF_978045245.1 (Mus_caroli_1) | 70 | 15 |
| `mus-pahari` | donor of the mouse host | GCF_900095145.1 (PAHARI_EIJ_v1.1) | 72 | 14 |
| `cricetulus-griseus` | the Chinese hamster (host of the hamster cybrids; donor in two mouse-host rows) | GCF_003668045.3 (CriGri-PICRH-1.0) | 70 | 16 |
| `mus-spretus` | donor of the mouse host | none | not assessed: no ortholog source | |
| `mus-terricolor` | donor of the mouse host (the paper's *Mus dunni* line, taxon 254704) | none | not assessed: no ortholog source | |
| `otomys-irroratus` | donor of the mouse host | none | not assessed: no ortholog source | |

- **No source.** NCBI holds no RefSeq annotation of a nuclear genome for these three. *Mus spretus* has GenBank assemblies
  (SPRET_EiJ v1 and v3) but no RefSeq one, and NCBI Gene holds only its 13 mitochondrial genes; *Mus terricolor* and
  *Otomys irroratus* have no genome assembly at NCBI and only their mitochondrial genes in NCBI Gene. Nothing stands in for
  them: not the other *Mus* species, not a GenBank-only assembly. (Another database annotates the *M. spretus* genome; that was
  not looked at and is not used.)
- **Hamster host.** `cybrid_measurements.tsv` holds no row with the Chinese hamster as host: `cricetulus-griseus` and
  `mesocricetus-auratus` appear only as donors into the mouse, both "no respiring cybrid". So there is no stated outcome for a
  hamster-host cybrid to run the swap on. The Chinese hamster's orthologs are fetched all the same, as the brief asked.
  The Syrian hamster (`mesocricetus-auratus`) is in the vertebrate alignments and was not fetched.
- **Gaps (4 of 360).** No gene of the species is in NCBI's ortholog group of the human gene, so none is taken:
  `mus-caroli` ATP5F1C (human gene 509) and ATP5MF (9551); `cricetulus-griseus` NDUFB1 (4707) and NDUFB9 (4715). NCBI Gene does
  hold genes of those species that look like them (looked up by symbol or description, **not taken and not checked as
  orthologs**): hamster LOC100754039 "NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 1", on the reference assembly;
  hamster LOC100759795 (NDUFB9), annotated only on the older CriGri_1.0; *M. caroli* Atp5mf (110294250), annotated on no current
  assembly; no *M. caroli* gene was found for ATP5F1C. Whether to take any of them is the lead's call.
- **No paralog ties.** Of the 360 (species, gene) pairs none had two genes of the species in the group: NCBI's ortholog group
  holds at most one gene per species (a human gene with several paralogs in the mouse, such as SERPINA1, has no mouse gene in
  its group at all). The "ambiguous (paralogs)" status is coded and tested on two real mouse genes put in one group, and
  occurs in no row.

## Rules

1. A species has a source when its RefSeq assembly with category "reference genome" exists (`donor_species.tsv`).
2. The ortholog of the human gene is the one gene of the species in NCBI's ortholog group of it (NCBI Datasets
   `gene/id/<id>/orthologs`, the table `gene_orthologs` of NCBI Gene). No gene: `no ortholog`; two or more: `ambiguous (paralogs)`
   and none is taken, with the candidates listed; a gene that is not protein coding or not annotated on the reference assembly:
   not taken (no row has any of these either).
3. The protein is the **longest** RefSeq protein of that gene on the reference assembly (`proteins` counts them). Proteins of
   equal length: `NP_` before `XP_`, then the lower accession number, then the lower version, so the choice does not depend on the
   order NCBI lists them in; `tied_with` names the ones passed over and `tied_differ` counts those with another sequence (1 of
   356: mouse UQCR10, NP_932096.1 chosen over NP_001349808.1, both 64 residues). The Chinese hamster gene records carry the proteins of two
   assemblies; only those of GCF_003668045.3 count.
4. Every sequence is checked against the length NCBI's product report gives, its SHA-256 is in the table, and
   `mitoage download donor-orthologs` refuses a sequence that is not the pinned one and writes nothing.

## Files

| file | content |
| --- | --- |
| `donor_orthologs.tsv` | the accession table: one row per (species, gene), 360 rows: status, human gene id, the species' gene id and symbol, RefSeq protein accession.version, length, SHA-256 of the sequence, how many proteins the gene has on the assembly, the ties, the candidates of an ambiguous gene |
| `donor_species.tsv` | the 8 species, their taxon, the reference assembly or none |
| `placement_check.tsv` | the reading of the swap (`cybrid_swap.protein_subs`) on every nuclear chain of the four solved structures, for each species: residues placed, not placed, substituted; nothing was scored |

Sequences are in `data/input/donor-orthologs/<species>.faa` (not in git; `mitoage download donor-orthologs`), with these
SHA-256 digests (second download into another folder: identical bytes; second `download donor-orthologs-select`: the same table
and species table byte for byte):

| file | proteins | SHA-256 |
| --- | --- | --- |
| `cricetulus-griseus.faa` | 70 | `6e7b6bc7d007496ec85cce53d887555e4a0b5a62d00fcda1acf258830032e889` |
| `mus-caroli.faa` | 70 | `bca5a5d8a5533c8248c7b6c7494231cf1f3eb3ea9b66fdaebc68cb141c39e249` |
| `mus-musculus.faa` | 72 | `cf67c7956e3963350884e39ccaa3fc88004b34dbb9241883f0892ef94b548fcf` |
| `mus-pahari.faa` | 72 | `9d5049c29cf57fcc01e0278d0c9e9e27c722fd4009aa9896a222a7719e013caf` |
| `rattus-norvegicus.faa` | 72 | `98c8d4a58da7d7974052ac841d40143bf166707a77726e02084aef76e205689e` |

## Checks

- Against the alignments: for the mouse and the rat, 45 of the 50 genes with an alignment have a RefSeq protein identical to the
  alignment row (gaps removed); the other five are another isoform of the same gene (the row is a shorter one, or a slightly
  different one: mouse NDUFS7, NDUFA12, NDUFB5, NDUFB9, COX4I1; rat NDUFS5, NDUFS7, UQCRQ, COX4I1, ATP5PD). The longest-isoform
  rule takes the longer one.
- The reading of the swap on the real structures (`placement_check.tsv`, 471 chain-species pairs over every nuclear chain of the
  four solved structures): of 75,368 residues 75,320 are placed and 48 are not (5 to 15 per species, at most 14 in one chain: a
  stretch of the structure the donor's protein lacks), and 10.6% of the placed residues differ from the human chain's.
- Tests and breaks (job 6629122 at commit `3034258`; the break scripts of `2a18c31`): the tests of the swap, its model and the donor
  orthologs, 105, all pass (none skipped; they read the real structures, records, build, alignments and the downloaded proteins). The
  rest of the suite: 407 pass, 177 skip, and 43 fail or error on inputs this checkout lacks (`tests/test_peptides.py`: `genome_features.parquet`;
  `tests/test_partners.py`: `composition.tsv`), not on code. Breaks: `slurm/h36-orthologs-mutations.py` 60 of 60 caught (the guards of the
  parsing, the choice, the table, the download and `read`, and of the swap's reading of a protein), `slurm/h36-swap-mutations.py` 64 of 64
  caught (the 63 that applied in the job, and the one whose text had become two lines, run alone after fixing it). The first run of the
  new script let three breaks survive (an unknown of the chain counted as a difference; the given number of differences ignored; only the first
  copy of a subunit given the donor's protein); each got a test and the rerun caught them.

## Not done

The swap was not run, for any host or donor. `cybrid_swap.run` still places the human structure and the human chain as the baseline:
a mouse or hamster host needs its own nuclear ortholog on the human chain as the baseline (B0), with the donor's as B1; the reading
function (`protein_subs`) serves both, the host's proteins being in the table, but the baseline and the host's tasks are not coded.
