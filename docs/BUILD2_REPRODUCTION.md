# Reproducing MitoAge build 2

This document says how to go from MitoAge build 1 (the 2015 data) to MitoAge build 2 (2026-11-13) with the
code of this repository, which sources that takes, which rules turn them into the build 2 tables, and how
far the result agrees with the published build 2. Every number below comes from commands run on
2026-11-14 with the code at the commit this document is in.

- Build 2 was made by the TypeScript pipeline of the site repository,
  [dna-seq/mitoage](https://github.com/dna-seq/mitoage) (`bun run data`, `bun run align`), at commit
  `8ff350c504bb59c9171eaed729eb6cddaada0f0b` of `main` for the files cited here. Its own accounts:
  [README](https://github.com/dna-seq/mitoage/blob/main/README.md),
  [docs/BUILD2_REPORT.md](https://github.com/dna-seq/mitoage/blob/main/docs/BUILD2_REPORT.md),
  [docs/BUILD2_OPEN_QUESTIONS.md](https://github.com/dna-seq/mitoage/blob/main/docs/BUILD2_OPEN_QUESTIONS.md),
  [docs/DECISIONS.md](https://github.com/dna-seq/mitoage/blob/main/docs/DECISIONS.md),
  [plans/DATA_BUILD2_DECISIONS.md](https://github.com/dna-seq/mitoage/blob/main/plans/DATA_BUILD2_DECISIONS.md).
- The published build is the Hugging Face dataset
  [longevity-genie/mitoage](https://huggingface.co/datasets/longevity-genie/mitoage): build 2 at revision
  `653afe2a3b79e18ae892e579ded9acdb17e0c334` (top-level tables and `alignments/`), build 1 in its `build1/` folder.
- The inputs no upstream archive keeps are in the snapshot dataset
  [longevity-genie/mitoage-build2-sources](https://huggingface.co/datasets/longevity-genie/mitoage-build2-sources),
  revision `9d55d72648d0838c3b5f37fca6185f8b938b3272` (see "Sources").

The code here is Python, written for this purpose from the rules; it is not a translation of the TypeScript
line by line, and it reads the published build 2 only in `mitoage build2 compare` ("never fill a value from
the same source that checks it").

## Build 1 and build 2

**Build 1** is the database of Toren et al. 2016 (Nucleic Acids Res 44:D1262-D1265), released on
2015-07-01: 922 animal species with a complete RefSeq mitogenome and a maximum lifespan from AnAge (build
12 era), with base composition per region and codon usage per protein gene. It is published, frozen, as the
`build1/` folder of the Hugging Face dataset: the full CSV export of 2015 converted to parquet without
changing a value (34 files, one per CSV). The export has 923 species; *Loxioides bailleui* (HAGRID 1598,
NC_025626.1, 9,243 bp) is not in the release. Its base columns are complement-labelled (see the rules).

**Build 2** re-checks the 922 against current sources and adds species: 2,335 species = 922 from build 1
(392 original, 465 updated, 46 renamed, 19 retired) + 1,413 new. Its tables: `species`, `region_stats`,
`codon_usage`, `genome_features`, `lifespan_observations` (every lifespan value of every source),
`changelog` (one row per changed value of a build 1 species, with what explains it), `review_list.csv`,
and, outside the tables, the gene alignments and the dated species tree (`src/data/generated/phylogeny.json`
in the site repository).

## Sources

Each is downloaded by one command into `data/input/build2/` and recorded in `data/input/MANIFEST.json`
(source, identifier, revision or md5, UTC time, files, bytes, and whether it matches build 2's pin).
`--from-snapshot` reads the snapshot dataset's copy instead (no token; checked against the snapshot's
`MANIFEST.json` and against build 2's pin). Pins are in `src/mitoage_pipelines/build2/pins.py`.

| source | what build 2 read (pin) | licence | command | matched the pin on 2026-11-14 |
| --- | --- | --- | --- | --- |
| MitoAge build 1 and build 2 | `longevity-genie/mitoage` @ `653afe2a3b79e18ae892e579ded9acdb17e0c334` | build 1: free with citation of Toren et al. 2016; build 2: mixed (FishBase/SeaLifeBase values CC BY-NC 4.0) | `download hf` | yes (revision) |
| NCBI RefSeq mitochondrion | release 237, `mitochondrion.1.genomic.gbff.gz` of 2026-09-03, md5 `4491799651740f49c8f6a011cf9a96e5` | no use restrictions (NCBI) | `download refseq --from-snapshot` (or `--from-file`, `--efetch`); `download refseq` takes the current release | NCBI's current release is **238** (2026-11-13, md5 `11660db0524ea3a09468b969bb7cf31f`): release 237 exists only in the snapshot |
| NCBI Taxonomy | `new_taxdump.tar.gz` current on 2026-11-13 (Last-Modified 2026-11-13T00:42:20Z), md5 `91f95787dd7c7cfd7042ab7610f2b4b2` | no use restrictions; cite Schoch et al. 2020 | `download taxonomy --from-snapshot`; `--archive 2026-11-01` (the closest monthly archive) or `--current` | the dump of 2026-11-13 is in no NCBI archive (archives keep the 1st of each month); snapshot copy has the pinned md5 |
| AnAge | build 15 (2023-07-03), `dataset.zip` md5 `7ec4d19f6c8f57b22595a9ffd50ea66b` | HAGR conditions of use; cite Tacutu et al. 2013 | `download lifespan` | yes (md5); the URL is unversioned, a change is reported as a warning |
| FishBase, SeaLifeBase | v26.06 parquet on source.coop `cboettig/fishbase` (`fb/`, `slb/`): `species`, `popchar`, `popgrowth`, `synonyms` | CC BY-NC 4.0; cite Froese & Pauly; Palomares & Pauly | `download lifespan` | yes (md5 of all 8 files) |
| USGS Bird Banding Laboratory | NABBP Longevity Records v2025.1, the table of `Longevity.php` | U.S. Government work | `download lifespan` (page saved, table written as CSV) | table: yes (CSV md5 `a84743335e7f080d5891c2f719f53702`); the page's HTML differs from the copy build 2 saved, its table does not |
| EURING longevity list | requested list of 2023-09-01: **not retrieved** by build 2 (euring.org in maintenance) | not checked | `download lifespan` tries it and writes `euring/status.json` | still redirected to maintenance.euring.org on 2026-11-14; no parser exists, nothing is substituted |
| RefSeq records cited by build 1 | 922 accession.versions, efetch `gbwithparts` | NCBI | `download cited` | byte-identical to build 2's cached copies (922/922) |
| CONTIG-only records | release records with no sequence (6 used); efetch of the same version, checked against the release length | NCBI | `download contig` | used the same 6 |
| INSDC source records | records named in the COMMENT of a cited record, for the re-annotation test | NCBI | `build2 primary-list`, then `download primary` | see "Rules: changelog" |
| Old site ids | `www.mitoage.info` export pages (24) | MitoAge | `download legacy` | all 922 ids recovered |
| Open Tree of Life | synthetic tree opentree16.1: TNRS and induced subtree | CC0 1.0; cite Hinchliff et al. 2015 | `download tree-sources --from-snapshot` (build 2's answers) or `download tree-sources` (live) | snapshot: the 2,241 answers build 2's tree read |
| TimeTree 5 | public API `mrca/id/<ids>/json` | free public resource; cite Kumar et al. 2022 | same | same |
| Build 2's dated tree (for the comparison only) | `src/data/generated/phylogeny.json` of dna-seq/mitoage @ `8ff350c5` | as the site | `download tree --ref <sha>` | — |

**The snapshot dataset** `longevity-genie/mitoage-build2-sources` (public, revision
`9d55d72648d0838c3b5f37fca6185f8b938b3272`, 18 files, 602,779,343 bytes) holds: RefSeq release 237 as
build 2 read it; the taxonomy dump of 2026-11-13; the cited, CONTIG-only and INSDC records as gzipped
GenBank files; AnAge build 15 `dataset.zip`; the USGS page and table; the 2,241 OTL/TimeTree answers of the
tree (`tree/answers.tar.gz`, file names are the hashes the site's pipeline gives each request); the old
site's pages; and a mirror of the 8 FishBase/SeaLifeBase tables read. Its card lists for every file the
original URL, the retrieval date, md5 and SHA-256, the licence and the citation. The NCBI monthly taxonomy
archive is not mirrored. The taxonomy dump of 2026-11-13 is, because no archive has it and the closest
archive changes the species set (below). Building and uploading the snapshot is a maintainer command
(`mitoage snapshot build`, then `HF_TOKEN=... mitoage snapshot upload`); reviewers only download.

## Rules, in plain words

The pointer after each rule is the file of the site repository that defines it; the Python is in
`src/mitoage_pipelines/build2/`.

**Build 1 as the starting point** (`pipeline/build.ts`, `pipeline/steps/regions.ts`; here `build1.py`,
`assemble.py:build1_view`).
- 922 species: the export without *Loxioides bailleui*.
- Base-label fix: in every `*_base_composition` file of build 1 the column `G` holds the number of C on
  the RefSeq strand, `C` holds G, `A` holds T and `T` holds A. Build 2 reads a = column T, c = column G,
  g = column C, t = column A. GC% is unaffected; GC and AT skews change sign against the old site.
- A species' build 1 AnAge entry: the AnAge build 15 entry of the same scientific name, else the only
  entry with the same common name.

**Complete mitogenomes, one per species** (`pipeline/sources/refseq-release.ts`, `pipeline/steps/refseq.ts`;
here `refseq.py`).
- Metazoa only. Complete: an `NC_` record whose definition is not partial, nearly complete or one of
  several molecules, and says complete genome or complete mitochondrial genome, complete sequence, or is a
  circular organelle assembly.
- The species of a record: its taxid made current through `merged.dmp`, then its nearest ancestor of rank
  species. Several complete records of one species: newest LOCUS date, then full length (not shorter than
  75 % of the median of its NCBI class), then the higher accession.version.

**Status of the 922** (`pipeline/steps/build2.ts`; here `assemble.py:plan`).
- Same accession.version in the release: live (884). Same accession, another version: the new version (2).
  Accession not in the release: the complete record of the same NCBI species, one at species rank before a
  subspecies one, never one another build 1 species cites (18 successors); none: **retired** (17), kept
  with its build 1 rows. *Carduelis pinus* is retired: its only candidate is *Spinus spinus*' record, cited
  by another row. *Andrias davidianus* is retired: AnAge build 15 has no maximum longevity for it.
- Two build 1 rows on one NCBI species: the row whose AnAge name gives another species takes that species.
- Renames come from NCBI only: the NCBI scientific name of the species replaces the build 1 name when the
  change is a genus transfer, an epithet spelling or gender change, or both; the old name goes to
  `previous_scientific_names` (46 renamed). Order, family, phylum and common name follow AnAge build 15;
  `class` keeps build 1's scheme for every species (AnAge build 12-era labels: Actinopterygii, Reptilia,
  Cephalaspidomorphi for lampreys, Sarcopterygii for lungfish and coelacanths).

**Lifespan sources and matching** (`pipeline/sources/{anage,fishbase,usgs-bbl,euring,observation}.ts`,
`pipeline/steps/lifespan.ts`; here `lifespan.py`, `taxonomy.py`).
- Values: AnAge maximum longevity (record type from specimen origin: wild, captivity -> captive, else
  observed; quality from data quality); FishBase/SeaLifeBase `popchar.tmax` (observed; quality questionable
  when tmaxQuality is -1), `species.LongevityWild` (wild), `species.LongevityCaptive` (captive),
  `popgrowth.tmax` (estimated, quality low); USGS "Minimum age at encounter" (wild; years + months/12).
- Matching a source name to an NCBI species, in order: NCBI scientific name, NCBI synonym (synonym,
  equivalent name, genbank synonym), the source's own synonym table (FishBase and SeaLifeBase; misapplied
  names never), the binomial of a trinomial. Names compare after removing parenthesised parts and quotes,
  collapsing spaces and lower-casing. A name that gives several species is ambiguous and not used; several
  records of one source on one species keep the one with the strictly best method, a tie makes them all
  ambiguous. Unmatched and ambiguous names go to the review list.
- Precedence, per taxon: AnAge is always first (a build 1 species takes AnAge through its build 1 link).
  Without an AnAge value: fish take FishBase; aquatic invertebrates SeaLifeBase (only species flagged
  fresh, brackish or salt water); birds EURING, then USGS BBL; other vertebrates nothing. Within a source:
  reported values that are not questionable, then questionable ones, then estimates; the maximum of the
  first tier that has one. **Tie-break:** at equal maxima a value with a stated specimen origin (wild,
  captive) beats one without (observed), then the first listed.
- Every value of every matched source is a row of `lifespan_observations`, plus the build 1 value of each
  build 1 species; `is_chosen` marks the one in `species.max_longevity_years`.

**New species** (`pipeline/steps/build2.ts`). Every NCBI species with a complete mitogenome, not linked to a
build 1 row (by its record, its species or its AnAge entry), with a chosen lifespan; an AnAge entry already
linked to a build 1 species is not added again. Id = 10,000,000 + NCBI species taxid (AnAge renumbered
HAGRIDs since 2015). Slug from the NCBI name, with the accession appended on a clash. Class by the build 1
scheme from the record's lineage, else AnAge's class.

**Genome values** (`pipeline/lib/genbank.ts`, `pipeline/lib/names.ts`, `pipeline/steps/genome.ts`,
`pipeline/steps/regions.ts`, `pipeline/lib/genetic-code.ts`; here `genome.py`).
- Every live species is recomputed from its record in the release (CONTIG-only records from efetch of the
  same version); retired species keep their build 1 rows.
- Feature names are normalised (ND1..CYTB from gene, product, standard name, note or locus tag;
  rRNA-12S/16S; tRNA-Xxx, Leu and Ser by isoacceptor from the anticodon, the codons recognised, the name or
  the numbering, and in chordates by the gene order: Leu(UUR) follows 16S rRNA, Ser(UCN) follows COX1).
- **join() rule:** a feature's bases are those of its annotated segments only (the ND3 of birds and turtles
  with its untranslated frameshift nucleotide is not counted over its start..end span).
- **D-loop:** the D-loop features; a misc_feature called control region (not a pseudo, partial or repeat
  one) only when the record has no D-loop feature. `non_coding_bp` counts the positions in no CDS, rRNA or
  tRNA, so the D-loop is non-coding.
- Regions: total_mtDNA is the sequence; total_pc, total_rRNA, total_tRNA, D_Loop are the union of their
  features' positions; a single gene is its first copy in genome order. Counts on the RefSeq strand;
  GC% = 100 (G+C)/(A+C+G+T), skews (G-C)/(G+C) and (A-T)/(A+T), rounded half up as JavaScript does.
- Codons: the coding sequence in its own sense from codon_start, complete codons only; total_pc is the sum
  of the gene rows; amino acids with the record's genetic code (the transl_table most CDS state; else 2 for
  chordates, 9 for echinoderms, 5 otherwise).

**Changelog and provenance** (`pipeline/steps/diff.ts`, `pipeline/steps/build2.ts`; here `diff.py`,
`assemble.py`, `validate.py`).
- Every value of a build 1 species that differs from build 1 is one changelog row with the old and new
  value and its `diff_class`. For a single region, tried in order: `rule_full_span` (build 1 equals the
  full span of a join() feature), `upstream_reannotation_same_version` (build 1 equals the INSDC record
  the RefSeq COMMENT names), `..._probable` (build 1 equals the gene feature's span), 
  `coordinates_differ_origin_unresolved` (build 1 equals a window <= 60 bp away), else
  `build1_disagrees_with_cited_record`; aggregates and total_pc codons inherit (`derived_from_...`); a
  changed accession makes every difference `refseq_successor` or `upstream_version_change`. The rule-wide
  corrections (base labels, join() spans, the excluded 923rd species) are logged once with `scope = rule`.
- **INSDC records:** the build 2 pipeline used the INSDC records its audit had cached (160). Here the choice is
  written down: `build2 primary-list` assembles once without any and lists the build 1 species that keep
  their accession.version and whose region or codon values differ from build 1 (151 on the pinned inputs);
  `download primary` fetches their INSDC records. With these 151 every changelog class equals build 2's.
- `record_status`: `new` for new species; `retired`; `renamed` when the name changed; `updated` when any
  changelog row names the species; else `original`. Child rows (regions, codons, features) take `new`,
  `retired`, `updated` or `original` from their own changelog rows. `assemble` stops if the species'
  status and the changelog disagree.

**Alignments** (`pipeline/align.ts`, `pipeline/steps/alignments.ts`, `pipeline/lib/align.ts`; here
`align.py`). The 13 proteins as RefSeq translates them and the two rRNAs, first copy per record, species
sorted by slug, MAFFT v7.526 FFT-NS-2 (`--retree 2 --maxiterate 0 --quiet --thread 4`), checked to hold
exactly the input sequences; the codon alignment the protein alignment implies, the stop codon left out.

**Dated tree** (`pipeline/sources/tree.ts`; here `tree.py`). Names to Open Tree ids by TNRS (exact, context
Animals); the induced subtree of the synthetic tree, unknown ids dropped; unary nodes collapsed,
polytomies kept; each internal node dated by TimeTree (up to three taxids from each of its first 12 child
clades), ages made monotone, undated nodes midway between parent and oldest child.

## Commands

```sh
uv sync && conda env create -f environment.yml      # MAFFT 7.526 for `build2 align`

# inputs (pinned: the snapshot where no upstream archive keeps the file)
uv run mitoage download hf                           # build 1 and the published build 2, revision 653afe2a3b
uv run mitoage download refseq --from-snapshot       # RefSeq release 237 (NCBI keeps only the current release)
uv run mitoage download taxonomy --from-snapshot     # the taxonomy dump of 2026-11-13
uv run mitoage download lifespan                     # AnAge, FishBase, SeaLifeBase, USGS BBL; tries EURING
uv run mitoage download cited                        # the 922 RefSeq records build 1 cites (efetch)
uv run mitoage download contig                       # records with only a CONTIG line in the release (efetch)
uv run mitoage download legacy                       # the old site's species ids

# steps
uv run mitoage build2 refseq                         # complete mitogenomes, one per species
uv run mitoage build2 match                          # lifespan-source names -> NCBI species
uv run mitoage build2 status                         # the 922: live / new version / successor / retired; candidates
uv run mitoage build2 primary-list && uv run mitoage download primary
uv run mitoage build2 assemble                       # the seven tables -> data/output/build2
uv run mitoage download tree-sources --from-snapshot && uv run mitoage build2 tree
uv run mitoage build2 align
uv run mitoage download tree --ref 8ff350c504bb59c9171eaed729eb6cddaada0f0b   # build 2's tree, for the comparison
uv run mitoage build2 compare                        # -> data/output/build2/compare
```

`build2` steps take `--release` (237, another release number, or `efetch`) and `--taxonomy` (a folder
under `data/input/build2/taxdump`). Run times of this sequence in a fresh data folder on 2026-11-14 (the machine has 32 cores; the
sandbox reported 8 to `nproc`; MAFFT ran with 4 threads):

| command | wall time |
| --- | ---: |
| `mitoage download hf` | 13.3 s |
| `mitoage download refseq --from-snapshot` | 15.5 s |
| `mitoage download taxonomy --from-snapshot` | 23.7 s |
| `mitoage download lifespan` | 19.8 s |
| `mitoage download cited` | 13.2 s |
| `mitoage download contig` | 8.4 s |
| `mitoage download legacy` | 55.3 s |
| `mitoage build2 refseq` | 21.4 s |
| `mitoage build2 match` | 21.8 s |
| `mitoage build2 status` | 22.8 s |
| `mitoage build2 primary-list` | 51.6 s |
| `mitoage download primary` | 3.3 s |
| `mitoage build2 assemble` | 52.9 s |
| `mitoage download tree-sources --from-snapshot` | 2.4 s |
| `mitoage build2 tree` | 0.4 s |
| `mitoage build2 align` | 100.5 s |
| `mitoage download tree --ref 8ff350c504bb59c9171eaed729eb6cddaada0f0b` | 1.0 s |
| `mitoage build2 compare` | 1.6 s |
| total | 429 s |
| `mitoage download tree-sources` (live OTL and TimeTree instead of the snapshot) | 1209.9 s |
| `mitoage download refseq --efetch` (build 2's 2,316 versions) | 49.4 s |

## How far it reproduces

**Pinned inputs** (the commands above, in a fresh data folder): **every value is identical.**

| table | rows published | rows reproduced | same keys | fields | values compared | identical | within tolerance | differ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| species | 2,335 | 2,335 | 2,335 | 49 | 114,415 | 114,415 | 0 | 0 |
| region_stats | 46,309 | 46,309 | 46,309 | 17 | 787,253 | 787,253 | 0 | 0 |
| codon_usage | 32,664 | 32,664 | 32,664 | 94 | 3,070,416 | 3,070,416 | 0 | 0 |
| genome_features | 88,647 | 88,647 | 88,647 | 12 | 1,063,764 | 1,063,764 | 0 | 0 |
| lifespan_observations | 7,197 | 7,197 | 7,197 | 19 | 136,743 | 136,743 | 0 | 0 |
| changelog | 14,392 | 14,392 | 14,392 | 9 | 129,528 | 129,528 | 0 | 0 |
| review_list | 95 | 95 | 95 | 3 | 285 | 285 | 0 | 0 |
| **all** | 191,639 | 191,639 | 191,639 | | 5,302,404 | 5,302,404 | 0 | 0 |
| alignments (28 FASTA files) | 28 | 28 | 28 | SHA-256 | 28 | 28 | 0 | 0 |
| dated tree (`phylogeny.json`) | 1 | 1 | 1 | bytes | 1 | 1 | 0 | 0 |

Keys: species `id`; region_stats and codon_usage `species_id, region`; genome_features
`species_id, feature_key`; lifespan_observations `species_id, source, source_record_id`; changelog
`table, row_key, field, old_value` (numbered when repeated); review list `review_type, source,
source_name, scientific_name`. Floats are compared exactly (tolerance 1e-9 was not needed).
Alignments: 28 of 28 files have the SHA-256 of the published uncompressed FASTA, and the lists of absent
genes, duplicated genes and codon notes are equal. Dated tree: `phylogeny.json` is byte-identical
(SHA-256 `64c7265e834faa8dfc750d7a1c7ea97004b60ee73177a522fa539aac336c13b7`; opentree16.1, 2,305 leaves,
2,229 clades, 2,045 nodes dated by TimeTree). The new species by class and lifespan source equal
`docs/BUILD2_REPORT.md` (Actinopterygii 236 AnAge + 453 FishBase, Aves 245 + 26 USGS BBL, Mammalia 237, ...;
AnAge 838, FishBase 475, SeaLifeBase 74, USGS BBL 26).

**With other inputs.** The first two variants change one input; the third takes every source as it is today. `compare` writes every disagreement with a reason
class (`species_in_one_build_only`, `input_label` — the text names another input version, and equals
once the names are swapped —, `same_version_record_differs`, or `unexplained`).

| variant | species (published 2,335) | values compared | identical | differ | reason classes of the differences | tree |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| Taxonomy: monthly archive 2026-11-01 instead of the dump of 2026-11-13 | 2,333 (3 only published, 1 only reproduced) | 5,295,636 | 5,295,546 | 90 values; 324 rows on one side | species_in_one_build_only 319, input_label 88, unexplained 7 (explained below) | not run |
| RefSeq: build 2's 2,316 live versions by efetch instead of release 237 (`--release efetch`) | 2,335, same set | 5,302,302 | 5,288,073 | 14,229 values; 65 rows on one side | input_label 14,115, same_version_record_differs 152, unexplained 27 (explained below) | not run |
| Fresh: RefSeq release 238 (2026-11-13), taxonomy dump of 2026-11-14, AnAge/FishBase/SeaLifeBase/USGS as downloaded, live OTL and TimeTree | 2,335, same set | 5,302,368 | 5,288,056 | 14,312 values; 43 rows on one side | input_label 14,203, same_version_record_differs 152 | byte-identical |

What the classes stand for in each variant:

- **Taxonomy archive 2026-11-01.** NCBI changed the taxonomy between 2026-11-01 and 2026-11-13. With the
  archive, *Saurida elongata* (NC_025940.1, FishBase) and *Planiliza lauvergnii* (NC_017887.1, FishBase) are
  not new species, and the record NC_084096.1 is *Notropis nubilus* (taxid 156700, id 10156700) instead of
  *Miniellus nubilus* (taxid 3363638, id 13363638): 319 rows of other tables belong to these four species.
  88 changelog rows name the taxonomy version. The 7 rows left unexplained by the automatic classes follow
  from the same change: the review flag `source_match_in_review:FishBase:Pristis microdon` of species
  11577861 and five review-list rows (*Pristis microdon*, *Miniellus albizonatus*, *Alosa maeotica* twice,
  *Coregonus atterensis*: these names match differently under the older taxonomy), and the row count of `lifespan_observations` in its changelog row (7,190 against 7,197). This is
  why the snapshot keeps the dump of 2026-11-13.
- **RefSeq by efetch.** The species set and every value of every species are the same except *Pongo
  pygmaeus* NC_001646.1 (id 3224): NCBI re-annotated that version in place on 2026-10-08 (LOCUS date
  03-APR-2023 in release 237, 08-OCT-2026 in release 238 and in efetch), so its regions, codons, features
  and changelog rows differ (109 values and 43 changelog rows; build 2's open question 5 names it). 14,115
  values name the input. The 27 unexplained: *Andrias davidianus* (id 944) is retired for want of a lifespan
  in both, but its record NC_004926.1 is not among build 2's live versions, so the reproduction gives the
  "suppressed" reason (5 values); and the review list loses 22 items (20 unmatched source names whose genus
  has a complete mitogenome in the release, 2 RefSeq species excluded through a build 1 link): both need
  the whole release, which efetch of a version list does not give.
- **Fresh.** Release 238 has the same complete mitogenome for every species of build 2 and no new one that
  changes the set (2,335 species, the same 1,413 new); taxonomy changes between 2026-11-13 and 2026-11-14
  change nothing. The only genome difference is again *Pongo pygmaeus*' re-annotation (109 values, 43 changelog rows); every
  other difference is the name of the input version (release 238 and its date, the taxonomy dump of
  2026-11-14) in the changelog, `retired_reason` and provenance text. The live Open Tree and TimeTree
  answers of 2026-11-14 had the same content as build 2's (all 2,241; byte-identical when written as the
  site's code writes them), so the tree is byte-identical. AnAge, FishBase/SeaLifeBase v26.06 and the USGS
  table had build 2's checksums, so this run does not test a change of lifespan sources.

The comparison tables of each run are written by `build2 compare`: `rows.tsv` (row sets), `fields.tsv`
(per field: compared, identical, within tolerance, differs), `disagreements.tsv` (every disagreement with
its key, both values and its reason class), `disagreements_by_reason.tsv` and `summary.json` (with the
alignment and tree results).

## Known gaps

- **RefSeq release 237 and the taxonomy dump of 2026-11-13 are kept only in the snapshot dataset.** Without
  it, a reader has `--efetch` (build 2's versions, today's text of each) and the archive of 2026-11-01; the
  differences that brings are counted above.
- **EURING** was not retrieved by build 2 and is not now; no values are used and no parser exists. The
  26 birds whose lifespan comes from USGS BBL stay provisional (review flag), as in build 2.
- **AnAge and the USGS page are unversioned.** Their md5 matched build 2's on 2026-11-14; a later run
  warns when they change and then differs from build 2.
- **efetch returns NCBI's current text of a version.** The cited records were byte-identical to build 2's
  copies on 2026-11-14; *Pongo pygmaeus* NC_001646.1, re-annotated on 2026-10-08, already shows how a later
  re-annotation changes the rows when the record is read by efetch instead of from release 237.
- **The INSDC records for the re-annotation test** were chosen by hand in build 2's audit (160); here by a
  written rule (151). The changelog is identical with either set on these inputs; on other inputs the rule
  may choose differently from what the audit would have.
- **The possible-duplicate pairs** (three) are a fixed list from build 2's audit, as in the site's code.
- **The live tree** depends on the current Open Tree synthesis and TimeTree's answers (variant above).
- Correlations and `explorer.json` of the site (computed at site build) are not reproduced here.
