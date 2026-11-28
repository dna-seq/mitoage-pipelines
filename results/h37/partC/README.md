# H37 part C — the gene set, resolved before anything is tested (2026-11-28)

Part C asks whether the per-branch rate of a mitochondrial tRNA covaries with the per-branch rate of its own
nuclear synthetase more than with the others', the mispairing being the null H5 defines. **Nothing of that is run
yet.** This folder holds only the first step the entry names: *which* nuclear gene is "its own synthetase", read
from NCBI at run time rather than from the entry. It is written down now, with the table it produced, so that the
pairing cannot be adjusted later to suit a result.

The entry's precondition was met by part A: its rate driver reproduces H7's protein-gene table in all 655 rows.

## What was run

`mitoage mitoport synthetases` (`src/mitoage_pipelines/mitoport/synthetases.py`), on 2026-11-28. The NCBI answer it
read is kept outside the repository (61 MB of gene records); a rerun with `--xml` pointing at a saved answer
repeats the run on the same bytes, and without it the search is made again, which may give a different set as NCBI
re-annotates.

## The rule, and why it is this one

A search for the mitochondrial synthetases *by their own words* does not find them, and that is not a detail: each
of the obvious queries is wrong in a way that would have gone unnoticed in the table.

| query, human genes, alive | returns |
| --- | --- |
| `"tRNA synthetase 2, mitochondrial"[Title]` | 0 — NCBI does not index gene descriptions in `[Title]` |
| `"mitochondrial aminoacyl-tRNA synthetase"` | 6 of them |
| `"aminoacyl-tRNA ligase activity"` | 38, and **HARS2 is not among them**, nor is the gene for histidine found any other way |
| `"tRNA ligase activity"` | 42, which contain every one of them |

So the search is the widest of the four, and the narrowing is done on each record, not in the query:

1. **`mitochondrial`** — NCBI's own description ends in `, mitochondrial`. 16 genes.
2. **`dual`** — it does not, but the gene carries a mitochondrial cellular component *and* a GO process naming a
   mitochondrial aminoacylation. 2 genes (GARS1, YARS2).
3. **`ambiguous`** — a mitochondrial component with no such process. 2 genes (KARS1, QARS1). **Not used**, and
   written into the table as not assessed with its evidence. A mitochondrial component on its own is also what a
   cytosolic synthetase with one stray annotation looks like, and this rule does not tell the two apart.
4. **`other`** — no mitochondrial component. 22 genes.

**The amino acid is taken from the description, never from the GO process terms.** NCBI's process annotation for
DARS2, the aspartyl synthetase, reads `mitochondrial asparaginyl-tRNA aminoacylation`: at least one gene of the set
has a process term naming the wrong amino acid, while description and symbol agree. A test holds that record so the
rule cannot drift back.

## What it gave

**18 of the 20 amino acids have a gene** (`synthetases.tsv`): Ala AARS2, Arg RARS2, Asn NARS2, Asp DARS2, Cys
CARS2, Glu EARS2, Gly GARS1, His HARS2, Ile IARS2, Leu LARS2, Met MARS2, Phe FARS2, Pro PARS2, Ser SARS2, Thr
TARS2, Trp WARS2, Tyr YARS2, Val VARS2.

**Two are not assessed, and both were foreseen in the entry**, which is a check on the rule rather than a
coincidence: the entry says Gln has no direct synthetase and that two of the set are dual-localised.

- **Gln** — the mitochondrial Gln-tRNA is not charged by a synthetase at all but by the GatCAB
  amidotransferase, which has no tRNA ligase activity and so is not in the search. QARS1, the cytosolic
  glutaminyl synthetase, comes out `ambiguous` on its mitochondrial-matrix component and is not used.
- **Lys** — KARS1 is dual-localised in the literature, but NCBI gives it mitochondrial components with no
  mitochondrial aminoacylation process, so this rule cannot separate it from QARS1. It is not used. The cost is
  stated here rather than paid quietly: part C will test 18 pairs, not 20, and tRNA-Lys and tRNA-Gln are left out.

## The orthologs, and a defect found in taking them (2026-11-28)

The 18 ortholog sets are downloaded over Mammalia (`mitoage download orthologs`, 2026-11-28; the files are outside
the repository). Taking one protein per species out of them turned up something that is not part C's own and is
reported rather than worked around.

`mitoport/select.py` picks, for each species, the isoform whose length is nearest the human one — and took "the
human one" to be **the first human RefSeq protein in NCBI's file order**. That order is by accession and means
nothing biologically: a gene whose alternative isoforms were accessioned after its canonical one (`NP_001…` sorts
before `NP_0…`) anchors on an isoform. Measured on this machine's 72 ortholog sets, the first human protein is not
the one NCBI tags `[isoform=1]` in 11 of the 18 synthetases and in 20 of the 54 genes of H5's set. Where the
anchor lands on a short isoform the whole column follows it: CARS2 anchored at 302 residues against a canonical
564, and 161 of its 181 species then came out more than 20% from the anchor.

So `select` now takes an anchor. `first` is the old rule and stays the **default, so that no published number
moves**; `longest` takes the longest human RefSeq protein, which depends neither on file order nor on an isoform
tag — NCBI writes `[isoform=1]`, `[isoform=a]`, `[isoform=X1]` or nothing, so no tag rule covers every gene — and
the fullest annotated form is what an alignment of orthologs is for. **Part C uses `longest`, fixed here before any
rate is computed.** Over the 18 synthetases the species per gene do not change; what changes is which isoform each
species contributes, and the count of species more than 20% from the anchor falls from **396 to 28**, nearly all of
it in three genes (CARS2 161 to 5, RARS2 115 to 0, NARS2 92 to 2).

Whether the earlier mitoport runs that used the old anchor should be done again under the new one is not part C's
to decide and is not assumed either way: it is desk question #100.

## A third amino acid is not assessed, and this was decided before any rate existed (2026-11-28)

The test runs on the mammals H5's run used, so that the coupling the bar compares against is measured on the same
branches. H5's set is 62 mammals, and `mitoport align` wants a complete matrix: a species missing from one gene is
missing from all of them.

**HARS2 does not cover them.** Its NCBI ortholog set holds 164 genes where the other seventeen hold about 270, and
it reaches 38 of H5's 62 mammals where every other synthetase reaches 60 to 62. Keeping it would set the complete
matrix at **35 species**; dropping it gives **57**. That is 22 species taken off all seventeen other pairs to keep
one, and it is the whole weight of a pair that would itself rest on 38.

So **histidine is not assessed for want of orthologs**, and it is recorded here with its number rather than left
to fall out of a later step. Dropping YARS2 as well would give 59 instead of 57, which is not worth an amino acid,
so it stays.

**Part C therefore tests 17 pairs over 57 mammals**: 20 amino acids, less glutamine (charged by GatCAB, not by a
synthetase), lysine (KARS1's annotation cannot be told from a cytosolic synthetase's) and histidine (orthologs).
Each of the three is a different reason, and none of them is a result about rates.

## What this does not yet decide, and must be fixed before the test is run

Written here as open, so that none of it is settled after a number is seen.

- The tRNA side: part A's RNA rates are per sister pair, and part C needs a per-branch rate on the mitoport
  topology. Which of the two is used, and how a tRNA with too few informative columns on a branch is treated.
- The isoacceptors: Leu and Ser each have two tRNA genes and one synthetase. Whether the pair is the synthetase
  against each tRNA separately or against the two together.
- The null: H5's mispairing, drawn over the other 16 synthetases, and how many draws.
- The controls' set: the cytosolic proteins H5 uses as the nuclear background, for these genes too.

The bar itself is already fixed in the entry and is not restated here: the real pair must beat the 95th percentile
of its mispaired nulls, in more than half the classes tested, and the coupling must be at least as strong as the
protein-gene coupling H5 and H14 found on the same mammals. Below any of these it does not count.
