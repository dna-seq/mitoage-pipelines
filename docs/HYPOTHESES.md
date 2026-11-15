# Hypotheses, tests and what came out

Each entry: the hypothesis as we held it, the test, the control that keeps the test honest, the result, and
its limits. Results that did not hold are kept next to those that did: a reader should see what was tried.
Numbers are from the runs of 2026-11-13 on build 2 of MitoAge (H1 and H2's power and second-alignment
checks: 2026-11-14, dataset revision 6efa0959282e, on which `controls` and `composition` give the numbers of
2026-11-13 again) (about 2,300 species with a mitochondrial
genome and a maximum lifespan); the command that produces each table is named. "Lifespan" is maximum
recorded lifespan, log10 years. Nothing here is peer reviewed.

## Data

- The MitoAge build: species table, lifespans, alignments of the 13 mitochondrial proteins (amino acids and
  codons) and of the two rRNA genes. `mitoage download build`.
- Nuclear orthologs of respiratory-complex subunits and of control proteins, from NCBI's ortholog sets.
  `mitoage download orthologs`.
- Solved structures of the human respiratory complexes from RCSB, each file's SHA-256 recorded.
  `mitoage download pdb`.

**Which build.** H1, H2, H7, H8 and H12 give the numbers of the build with the nine European birds from the EURING
list (longevity-genie/mitoage f75bdc9, the data of dna-seq/mitoage fb54cec, 2,344 species), rerun on 2026-11-15 by
`slurm/manuscript-rerun.sbatch` (job 6489877 from 2ee27a1; tables in `results/euring/`). Where a number changed
from the build before (c525adb, 2,335 species; job 6489876 reproduced the published tables), the earlier value
follows in parentheses; `results/euring/README.md` and `comparison.txt` list every change. The nine birds add two
species to the tested sets, and every gene was realigned with them, which moves alignment columns in every class.
Mammal results are byte-identical except H2's column screen, whose columns moved. H12's second round (partner
discordance and ESM-2 embeddings) and H2's power calculation were not rerun and describe the build before.

## H0. The alignments do not depend on the aligner

**Why it matters.** Every site-level result below reads alignment columns.
**Test.** MAFFT (fast, as shipped), MAFFT iterative, FAMSA 2, Muscle 5 on each protein; the share of residue
pairs the shipped alignment puts in one column that another tool also puts in one column.
`mitoage longevity compare-aligners`.
**Result.** At least 0.985 with every tool for COX1, COX2, COX3, CYTB, ND1. Soft: ND6 (0.87 to 0.97), ATP8
(0.92 to 0.99), ND5 (0.96 to 0.99).
**Consequence.** A site-level finding in a soft gene is kept only if it holds on a second alignment.
**Limit.** Agreement between tools is not correctness: there is no reference alignment.

## H1. Amino-acid content of mitochondrial proteins follows lifespan

**Hypothesis.** Published, and contested. Long-lived animals have fewer cysteines in their mitochondrially
encoded proteins (Moosmann & Behl 2008, Aging Cell 7:32) and less methionine (Aledo et al. 2011, Aging Cell
10:198); higher primates gained threonine in the transmembrane regions of these proteins (Kitazoe et al. 2008,
PLoS ONE 3:e3343, doi:10.1371/journal.pone.0003343). Against this, Aledo, Valverde & de Magalhães (2012, J Mol
Evol 74:332, doi:10.1007/s00239-012-9510-7), on 173 mammals, found no effect of longevity on cysteine content
and threonine usage fully accounted for by nucleotide composition. We test all 20 amino acids and read
cysteine, threonine and methionine against this record; none of them is a positive control.
**Test.** Content of each amino acid over the 13 proteins (percent) against log10 lifespan, with log10 adult
body mass, and then also with base composition at fourfold-degenerate sites (shares of A, C and T), as a
phylogenetic regression on the build's dated species tree, with Pagel's lambda estimated for each amino acid
and model. Slopes are percentage points per tenfold lifespan; q is Benjamini–Hochberg over the 20 amino
acids, within class and model. `mitoage longevity aa-pgls --class <Class>` (after `controls`, `fourfold` and
`dated-tree`), into `aa_pgls.<Class>.tsv`, which also holds the earlier mixed model (random intercepts of
order and family) refitted on the same species (method `mixed`).
**Species.** The manuscript's default dataset: not retired, lifespan record of high or acceptable quality and
not estimated; with body mass; a tip of the dated tree. Mammals: 625 with all 13 proteins and a lifespan,
615 in the default dataset (10 questionable lifespans out), 612 with body mass, 604 in the tree (8 species
are not in it). Birds 381 (379), ray-finned fishes 235, reptiles 97.
**Result, mammals (604).** Lambda is 0.96 to 1 for every amino acid: content follows the tree almost as
Brownian motion would. With body mass and composition **no amino acid passes q < 0.05**; the smallest q is
0.11 (alanine −0.13, p = 0.009; aspartate −0.045, p = 0.011; isoleucine −0.16, p = 0.016). Threonine +0.15
(p = 0.03, q = 0.12), cysteine −0.021 (p = 0.05, q = 0.15), methionine +0.004 (p = 0.94); histidine +0.05 and
glycine +0.035 (q = 0.12 and 0.13); the other twelve have p > 0.1. With body mass only five pass: isoleucine
−0.26 (q = 0.004), histidine +0.07, glycine +0.05, phenylalanine −0.09 and threonine +0.18 (q 0.03 to 0.045);
cysteine does not (−0.017, p = 0.10). Composition takes about a sixth of threonine's slope and two fifths of
isoleucine's.
**What changed the result: the tree, not the species set.** The earlier mixed model on the same 604 species
still gives seven amino acids at q < 0.05 with composition: isoleucine −0.35, cysteine −0.052 (q = 0.0006),
leucine +0.24, phenylalanine −0.15, aspartate −0.062, threonine +0.30 (q = 0.002) and alanine −0.17. Random
intercepts of order and family leave most of the shared ancestry in the residual; the dated tree, with lambda
near 1, takes it out: threonine's slope halves and cysteine's falls by more than half.
**Result, other classes, with composition.** Nothing at q < 0.05 in birds (smallest p 0.018; 0.019), ray-finned
fishes (0.005) or reptiles (0.004); cysteine and methionine have p ≥ 0.15 in all three. Threonine goes the
other way in fishes (−0.076, p = 0.02, q = 0.17; with body mass only −0.094, q = 0.03). Whatever is seen in
mammals is not shared across vertebrates.
**Reading.** On the dated tree the published directions are there in mammals (threonine up, cysteine down)
at about half the size the mixed model gave, and neither is beyond chance among 20 amino acids once
composition is in the model; methionine shows nothing. This agrees with Aledo et al. 2012 on cysteine, and
on threonine as far as composition and relatedness together go. It does not test Kitazoe et al.'s claim,
which is about transmembrane regions in higher primates, not whole proteins across mammals.
**Limits.** Composition is three shares of fourfold sites over the 12 heavy-strand genes, a coarse measure of
mutation pressure; the tree's polytomies are kept as they are; one lifespan per species.
**Earlier version (runs of 2026-11-13), kept as a lesson.** 622 mammals with body mass, including the 10 with
questionable lifespans; mixed model of order and family (`mitoage longevity controls`, `composition`):
threonine +0.30 (p = 0.0005) and cysteine −0.05 (p = 0.00008) were reported as positive controls; alanine
(−0.175, p = 0.0005) and aspartate (−0.06, p = 0.0006), as significant, were not reported. Without any control
for relatedness threonine has p = 1e-60 across all species; within families it is 0.9.

## H2. Single residues are associated with lifespan

**Hypothesis.** Some positions of the 13 proteins carry a residue more often in long-lived species,
independently in several lineages.
**Test.** Each alignment column and residue (present in at least 15 species and absent in 15, both states
in at least 5 families) against lifespan, body mass and fourfold composition in the model.
**First version, kept as a warning.** Comparison within families (`mitoage longevity sites`): 226 of 5,017
tests pass a 10 % false-discovery threshold, but the p-values are inflated over the whole range. Genera
inside a family are not independent.
**Fixed before any result.** The tree is the build's dated species tree (Open Tree topology, TimeTree
dates), never one made from these genes; lifespan records of high or acceptable quality only.
**Second version, and a mistake.** A phylogenetic regression with Pagel's lambda
(`mitoage longevity pgls`). It was first run on a tree built from the 13 proteins themselves and on every
lifespan record, against the rule above. That run gave cytochrome b position 43 (methionine, q = 0.03) as
"the one site". Rerun as fixed (dated tree, quality filter), that site has q = 0.29.
**Result as fixed, mammals (604).** P-values are flat. Two columns pass at 5 % false discovery, among 5,262 tests
(5,259 before the realignment; the same two columns pass on both alignments, with the same carriers):
- ND5 position 200 (human numbering): histidine in place of glutamine. It is the marsupials (24 of 28
  carriers): one old change in one clade. Not a longevity site on this evidence.
- COX1 position 419: isoleucine in 24 species with about 40 % longer lifespan given body mass
  (p = 2e-5, q = 0.049): great apes, elephants, bears, African mole-rats and single species in eight other
  families. On the gene-built tree it had q = 0.27.
**Result as fixed, birds (381; 379).** Nothing (smallest p 0.0007; 0.001), among 2,738 tests (2,716).
**Second alignment (H0's rule for the soft gene ND5).** Rerun on the FAMSA, Muscle 5 and iterative MAFFT
alignments of the same sequences (`mitoage download prepared --topic aligners`, then `mitoage longevity
pgls --alignments-dir data/interim/aligners --tool <tool> --out ...`): in each, both columns hold the same
residues in the same species (573 glutamine and 28 histidine at ND5 200, 24 isoleucine at COX1 419) and
come first and second, with q = 0.005 and 0.049. They stand on the alignment; what is said of them above
does not change.
**Power** (`mitoage longevity pgls-power`). With the standard errors the tests themselves have, a lone column
passes 5 % false discovery (p ≤ 0.05/5,259 = 9.5 × 10⁻⁶) with 80 % power only if its residue changes lifespan,
given mass and composition, by 0.17 log10 years (1.49-fold) with 20 to 29 carriers, 0.13 (1.35-fold) with 50
to 99, and 0.10 (1.27-fold) with 200 to 399. Half the tests have 85 carriers or fewer. COX1 419 (0.146 with
24 carriers) is below what the screen finds reliably.
**Conclusion.** No single residue is robustly associated with lifespan in mammals or birds: the one or two
columns near the threshold change place with the tree, and among about 1,650 columns that is what chance
gives. The screen can see only large effects, of about 1.3- to 1.5-fold lifespan. The gene-built tree remains available (`mitoage longevity concat`, `mitoage longevity tree`) only to
show this sensitivity.

## H3. Selection on mitochondrial genes differs between long-lived and short-lived lineages

**Hypothesis.** Purifying selection is stronger (better-protected mitochondria) or weaker (smaller
populations) in long-lived lines.
**Test.** HyPhy RELAX per gene; species long-lived and short-lived for their body mass, in equal numbers
from each end within each order, at most three per family; tips labelled. `mitoage longevity relax-prep`,
`mitoage longevity relax`.
**Control.** A second draw of species disjoint from the first.
**Result.** First draw (57 against 57): ND3 relaxed (K = 0.71, p = 0.0008) and cytochrome b relaxed
(K = 0.51, p = 0.003). Second draw (48 against 48, a weaker contrast by construction): neither repeats
(ND3 K = 0.89, p = 0.38; CYTB K = 1.03, p = 0.52); ATP8 appears there only. **No gene shows a difference
that holds across both draws.** We do not claim the first draw's hits.
**Provisional, not to be quoted: Coevol** (lifespan as a continuous trait on every branch, on the dated
tree; `mitoage longevity coevol-prep`, `longevity coevol --chain 1` and `--chain 2`, `longevity coevol-read`).
First run on the cluster on 2026-11-14 and 15 by another agent: two chains of 106 mammals (one per family),
1,100 cycles each, read after a burn-in of 200. The inputs `coevol-prep` writes are byte-identical to that
run's; the chains themselves were not rerun. They had not mixed when these numbers were read (tracecomp:
effective sample size 7 for dN/dS and 15 for the likelihood, relative difference between the chains 0.52 and
0.79; above 100 and small are needed before any of this is quoted), and with that caveat: dS against longevity −0.52 and −0.49 (posterior probability of a positive correlation 0 in
both chains), dN/dS against longevity +0.45 and +0.44 (0.99); as partial correlations, given the other three
variables (body mass among them), dN/dS against longevity is +0.28 and +0.32 (0.88 and 0.91) and dS against
longevity −0.13 and −0.10 (0.07 and 0.14), so the raw dS association went with body mass. The runs need to be
resumed to convergence; until then the result for lifespan as a continuous trait is not assessed.
**The tree Coevol saw.** `coevol-prep` writes the dated tree's topology without its branch lengths, so
Coevol was not given the ages: it set every branch to 0.001 and sampled the node ages itself, with no
calibration. "On the dated tree" above means its topology only. Sampling about a hundred ages jointly with the
rates is a known slow part of Coevol's sampler and a likely cause of the poor mixing. `coevol-prep
--fixed-ages` keeps the ages and `coevol --fixed-ages` passes `-fixtimes`, which holds them at the dated tree's
(scaled to a root age of 1); two such chains are run beside the resumed ones.
**Two-step cross-check (reading fixed on 2026-11-15, before it was run).** The same 106 mammals, alignment
and topology as Coevol. Step 1, `mitoage longevity branch-rates`: HyPhy fits MG94xREV with a synonymous rate
α and a non-synonymous rate β of its own on every branch (the tree unrooted, as HyPhy needs; the two branches
at the root become one); dS = α·(S̄+N̄) per synonymous site, dN/dS = β/α. Step 2, `mitoage longevity
branch-pgls`: on the 106 terminal branches, whose species' traits are observed, phylogenetic regression
(Pagel's λ estimated per response, the dated tree's shared paths) of log10 dN/dS and of log10 (dS per
million years of the branch's duration on the dated tree) on log10 lifespan and log10 body mass, each alone
and both together. A branch with α or β estimated at 0 has no logarithm and is left out and counted.
- **Positive control:** dN/dS rises with body mass (the nearly-neutral expectation; Popadin et al. 2007,
  mitochondrial genes of mammals). If the slope of log10 dN/dS on log10 mass alone is not positive with
  p < 0.05, the test lacks power on these data and nothing below is read.
- **Secondary control:** dS per million years falls with body mass (Nabholz et al. 2008). Reported, not
  required.
- **The question:** lifespan with mass in the model (lifespan alone is mostly mass). Two tests, dN/dS and dS
  per Myr, each two-sided, Holm-corrected over the two. A corrected p < 0.05 with the sign of Coevol's
  provisional partial correlations (dN/dS +, dS −) is called agreement with Coevol; the opposite sign at
  p < 0.05 is disagreement; anything else is "no lifespan signal beyond body mass in the two-step test".
- **Sensitivity, reported beside it:** the same models without terminal branches whose dS exceeds 2
  (saturated synonymous sites).
- **What it cannot do:** it ignores the uncertainty of each branch's α and β (Coevol does not), and uses tips
  only. It is a check of direction and size, not a replacement for the converged chains.
**Open.** The RELAX runs above used the gene-built tree for their species tree; `relax-prep` now writes the
dated one, and the runs have not been repeated on it.

## H4. Peptides encoded inside mtDNA are conserved across species

**Hypothesis.** The published mitochondrial-derived peptides (humanin, SHLP6, MOTS-c inside the rRNA genes;
MTALTND4 in another frame of ND4) are kept as open reading frames beyond human, and their presence follows
lifespan.
**Test.** Find each peptide in the human row, project its columns onto every species, ask for a start codon
and no stop before 90 % of the length. `mitoage peptides pilot`.
**Control.** Other windows of the same length in the same gene (and frame) that are stop-free in human.
**Result.** All four are found in human with the published sequence. Across species only SHLP6 looks kept
(start codon in 2,189 of the 2,191 vertebrates with the region: all but one ray-finned fish of 891 and one bird
of 418; frame open in 88 to 97 % per class); humanin is open in
89 % of mammals but keeps its ATG in 74 %; MOTS-c is no better than a random open window of 12S; MTALTND4 is
kept in 36 of 116 primates. **None is beyond the null**: 10 to 14 % of random 16S windows open in human are
kept as widely as SHLP6.
**Limits.** A kept frame is not a translated peptide; nuclear copies are out of sight; outside vertebrates
the windows do not line up. The test against the host gene's synonymous codons shuffled is H15's.

## H5. Mitochondrial and nuclear subunits of one complex evolve together (mitoport)

**Hypothesis.** Nuclear subunits that touch mitochondrially encoded ones change rate on the same branches,
more than nuclear proteins that do not; and this pairing differs with lifespan.
**Test.** Rate covariation: one fixed tree, branch lengths per protein, each nuclear protein's rate relative
to cytosolic control proteins, each mitochondrial complex's rate relative to the mean of the 13
mitochondrial proteins, correlation over branches. `mitoage mitoport select`, `align`, `branch-lengths`,
`erc`.
**Controls.** Cytosolic proteins; complex II (no mitochondrial subunit); mitochondrial matrix enzymes.
**Result (62 mammals with all 67 proteins, dated species tree).** Complex IV's eleven nuclear subunits follow
the mitochondrial COX genes (median correlation +0.14 against +0.01 for 23 cytosolic controls, p = 0.004),
and not cytochrome b (−0.17). Complex III's nine nuclear subunits against cytochrome b: nothing (−0.11
against −0.06). On the topology of the tree built from the proteins the same test gave p = 0.0001; the
direction and size are the same, the dated tree is the one that counts. This reproduces a known phenomenon
and serves as the method's positive control. Lifespan: H9.
**Residue level, complex IV only** (`mitoage download pdb 5Z62`, `mitoport contacts`, `mitoport burial`,
`mitoport iface-align`, `mitoport iface`; into `sites.tsv` and `coupling.tsv`). Every
residue of COX1, COX2 and COX3 in the human structure 5Z62 was mapped to its column of the MitoAge alignment;
the three chains match MitoAge's human sequences at identity 1.000. 298 residues touch a nuclear subunit
(4.5 Å heavy-atom cutoff), and each was matched to a non-contact residue of the same gene by its accessibility
in the isolated chain. Contact residues are exposed residues (median relative accessibility 0.41 in the
isolated chain, against 0.05 for the residues that are neither contact nor control) and are less conserved
than their matched controls (median identity to human 0.922 against 0.971, Mann-Whitney p = 1.5e-05). Column
occupancy over the 2,174 species with an alignment, a tip and a usable lifespan has median 1.000; 290 of the
298 contact columns are at 0.99 or above, the lowest (COX1 position 2) at 0.46. Then, per nuclear protein, the
correlation over pairs of mammals between its distances and the mitochondrial COX distances, patristic
distance partialled out, with two distances (differing alignment columns; ESM-2 650M per-residue embeddings
averaged over the column set); 61 proteins, 87 to 182 mammals each (166 to 181 for complex IV's subunits).
**The group contrast reproduces**: complex IV's eleven nuclear subunits +0.58 (sequence) and +0.42 (ESM)
against +0.03 and +0.05 for 29 cytosolic controls (p = 7.2e-05 and 1.1e-04); complex III's nine +0.38
(sequence, p = 0.0011); complex II, which has no mitochondrial subunit, +0.05 (p = 0.94); matrix enzymes
+0.03 (p = 0.56). **The interface prediction fails**: on the mitochondrial side, contact residues couple more
than their burial-matched controls by +0.049 for complex IV against +0.016 for the controls (sequence,
p = 0.43; ESM +0.002, p = 0.37); on the nuclear side, interface columns do no better than the rest with the
sequence distance (−0.022, p = 0.72). With the ESM distance the nuclear interface columns of complex IV do
couple more than its other columns (+0.051 against +0.004, p = 0.0096, one of many contrasts and uncorrected),
so the two distances disagree on the nuclear side and this is not assessed further. **Reading.** The coupling
is a property of whole proteins that share a complex with mtDNA products, not of the residues that touch: it
is no evidence of compensation at contacts. Complexes I, III and V have no contact map yet, so this is said
of complex IV only. **Reproduction.** First run on the cluster on 2026-11-14 by another agent; rerun from this
code on 2026-11-15: the contact map and `sites.tsv` are byte-identical, and `coupling.tsv` is identical when
read on that run's nuclear alignments. On freshly made alignments (MAFFT with two threads is not repeatable)
60 of 62 are identical and one protein's sequence correlations move (UQCRQ, +0.565 to +0.572), which changes
no group median. The `esm_` columns read the per-residue ESM-2 embeddings (`mitoage download prepared
--topic residues`).
**Traps met.** IQ-TREE re-places identical sequences unless told `--keep-ident`, which silently changes the
tree; NCBI's ortholog sets for CYC1, SDHC and SDHD held 21 to 48 sequences of another gene or a broken
model, filtered by identity to human; a test of nuclear rate against overall mitochondrial rate with one
shared denominator is invalid and is not reported.
**Open.** Contact maps for complexes I, III and V; the pairing against lifespan beyond complex IV's contact
residues.

## H6. Partners' protein-language-model embeddings move together across species (mitoport)

**Hypothesis.** If nuclear and mitochondrial subunits co-evolve, the distance between two species in the
embedding of a nuclear subunit follows the distance in the embeddings of its mitochondrial partners, beyond
what the time separating the species explains.
**Test.** ESM-2 650M, mean over residues, for 61 nuclear proteins and the 13 mitochondrial ones in up to 182
mammals (embeddings computed in the kron-fm project). Per nuclear protein: partial correlation over species
pairs of its cosine distances with those of a mitochondrial complex's subunits, given patristic distance on
the dated tree and its logarithm.
**Controls.** Cytosolic proteins; complex II; mitochondrial matrix enzymes. Because species pairs are not
independent, only the comparison between groups of proteins is tested.
**Result.** Complex IV nuclear subunits +0.34 with the mitochondrial COX proteins (cytosolic controls
+0.03, p = 0.0004), complex III nuclear +0.23 (p = 0.005), cytochrome c +0.41; complex II +0.05 and matrix
enzymes +0.02, like the controls. **Not complex-specific**: complex IV's nuclear subunits follow the
mitochondrial subunits of complexes III, I and V almost as well (+0.31, +0.27, +0.20).
**Reading.** Nuclear proteins that share a complex with mtDNA products change more where the mitochondrial
genome changes more; being mitochondrial without such a partner is not enough. The common cause is
probably the lineage's mitochondrial rate, which the tree's time does not remove.
**Open.** The coupling against lifespan; per-residue embeddings at contact and non-contact positions.

## H7. Lifespan predicts the substitution rate of the mitochondrial proteins

**Hypothesis.** Published, and contested in both directions: long-lived lineages have been reported with
slower mitochondrial substitution rates in bivalves (Mortz et al. 2021, Front Mol Biosci,
doi:10.3389/fmolb.2021.626042) and in fishes (Le Cam et al. 2026, Phil Trans R Soc B 381,
doi:10.1098/rstb.2025.0085), while a survey across vertebrates found that whether the association appears
at all depends on how the rate is normalised (Sterling, Zwonitzer & Havird 2026, GBE,
doi:10.1093/gbe/evag067). The caveat is the reason to test it here.
**Test.** The relative-rate test on cherries of the dated tree. Two tips sharing a parent have had the same
time since they split, so measuring each against one common outgroup makes the rate difference a contrast
that needs no divergence-time correction. dN and dS are Nei-Gojobori counts with the Jukes-Cantor
correction, in each species' own mitochondrial code, averaging codon pairs that differ more than once over
every pathway that avoids a stop. The contrast is regressed on the lifespan contrast through the origin,
per gene and over the concatenated 13. `mitoage longevity rates`, into `rates.tsv` and `rate_contrasts.tsv`.
**Control.** Four normalisations, because the choice is what moves results in the published record: the raw
difference, the difference per million years of the cherry's depth, the dN/dS ratio, and the raw difference
with the body-mass contrast in the model. A sign test on the same pairs asks the question directly. Pairs
whose dS exceeds 1 are flagged saturated and every test is reported with and without them.
**Species.** H1's default dataset: 1,371 species (1,369) with a tip in the dated tree, giving 339 cherries;
dS > 1 in 24.2% of the concatenated pairs.
**Result.** **Nothing survives on the concatenated genes, under any normalisation, in any class.** The
closest are mammals with body mass in the model (slope +0.113 per unit log10 lifespan, p = 0.021, q = 0.15,
115 pairs) and birds per million years (p = 0.071, q = 0.28, 95 pairs; p = 0.073, q = 0.26). Four of 655 tests pass
q < 0.05, all single genes, the same four on both builds, and their signs disagree: ND6 across all species
(-0.197, q = 0.0075; -0.219, q = 0.0076) and avian CYTB per million years (-0.022, q = 0.0081; -0.021,
q = 0.011) go one way, mammalian ND5 (+0.220, q = 0.010) the other, and mammalian ND2
on dN (-0.017, q = 0.039) is a fourth direction. Four of 655 at a 5% false discovery rate is what the
correction is designed to allow.
**Reading.** In these data lifespan does not predict the mitochondrial substitution rate, and that answer is
the same under all four normalisations. This does not refute the bivalve and fish reports, which are within
clades this build barely covers; it says the effect is not general across vertebrates.
**Limit.** A cherry gives one contrast, so the test uses 339 of the 1,371 species. Nei-Gojobori counting is
appropriate for close pairs and strained at the divergences a quarter of these cherries reach, which is why
the saturated flag exists rather than a silent filter.

## H8. mtDNA motifs, repeats and the deletion contact zones follow lifespan

**Hypothesis.** Non-B DNA motif counts per mitogenome were proposed as longevity features by the MitoAge
group itself (Toren et al. 2021, Mech Ageing Dev, doi:10.1016/j.mad.2021.111429), repeats have been reported
to track lifespan (Yang, Seluanov & Gorbunova 2013), and the repeats that flank mitochondrial deletions sit
in two defined zones of the genome (Persson et al. 2023, BMC Biol, doi:10.1186/s12915-023-01606-1).
**Test.** The features the site's pipeline publishes in `mtdna_structure.parquet` - G-quadruplex motifs on
each strand, homopurine and homopyrimidine tracts, homopolymer runs, di/tri/tetranucleotide tandem repeats,
direct, inverted and mirror repeats, and the repeat pairs bridging the ATP8-COX3 and ND5-CYTB spans - against
log10 maximum lifespan, in H1's model unchanged: phylogenetic regression on the dated tree with Pagel's
lambda, body mass and fourfold composition in it. Features are standardised, so a slope is per standard
deviation. `mitoage longevity motif-pgls`, into `motif_pgls.tsv`.
**Control.** Every test is run a second time on records whose sequencing technology includes long reads,
since mitoVGP and MitoHiFi both report that short-read reference mitogenomes misrepresent exactly these
features.
**Result.** 132 tests over four classes on the full build; **none passes a 5% false discovery rate**
(smallest q = 0.16, smallest p = 0.005). The strongest raw correlation in the set is the one the
literature reports: perfect inverted repeats against lifespan in 604 mammals, Spearman r = -0.399
(p = 2 x 10^-24), and over the major arc r = -0.334 in 530 (p = 3 x 10^-15). With body mass, fourfold
composition and Pagel's lambda on the dated tree, both are gone (p = 0.67 and p = 0.45). Repeat pairs
bridging the two contact zones behave the same way: r = -0.181 raw (p = 8 x 10^-6, 602 mammals),
p = 0.03 and q = 0.56 after the controls.
**What this replaced, twice.** An exploratory Python run over 2,311 species, with tandem-repeat counts
that did not reject units built from a shorter period, reported mammalian G-quadruplexes at q = 0.048
over 601 species. Rerun on the database's own features while they still covered only 922 species, that
result was gone (q = 0.66) and a reptile G-quadruplex test passed instead (q = 0.027). On the full
build that one is gone too (q = 0.38). Three runs, three different survivors, and none of them
survives the next run: the screen is reporting noise, and this entry reports that rather than any of
the three.
**Limit, and it is the important one.** **The control could not be run.** No record in the build states
a sequencing technology, so the long-read subset is empty and the mitoVGP caveat stays open: every
repeat, motif and contact-zone number above is computed on records of unknown assembly type.

## H9. Mismatch at complex IV's mito-nuclear contacts follows lifespan (mitoport)

**Hypothesis.** The port project's question narrowed to mitochondria: species whose COX1-3 residues at the
contacts with nuclear subunits differ more from the consensus than burial-matched non-contact residues do
have a different maximum lifespan.
**Test.** `mitoage mitoport iface`, into `lifespan.tsv` (with each slope's standard error and 95% interval);
first run on the cluster on 2026-11-14 by another agent, and the rerun from this code gives the same table.
Per species, the share of H5's 298 contact columns where its residue differs from the most common residue
of the column among the species of the set (gaps left out), the same at their 298 matched controls, and the
difference, each standardised; phylogenetic regression of log10 maximum lifespan on it, log10 adult body mass
in the model, Pagel's lambda estimated on the null model, the build's dated tree.
**Control.** The burial-matched non-contact residues, as their own regression and as the subtracted term.
**Result.** All species (1,373, lambda 0.96): slope of contact minus control −0.024 log10 years per standard
deviation, 95% interval [−0.050, +0.002], which bounds any effect to at most about an 11% shorter or a 0.5%
longer maximum lifespan per standard deviation of mismatch. Mammals (604): −0.006, [−0.028, +0.017]. The 176
species that also have nuclear partners: −0.007, [−0.049, +0.035]. Contact and control mismatch alone give
intervals that include zero in every set as well.
**Reading.** A bounded null: in these species an association of contact mismatch with lifespan larger than
about a tenth of maximum lifespan per standard deviation is excluded, and nothing smaller is claimed. Like H2, H7
and H8, a lifespan question whose answer here is no detectable effect. The
mitoport plan's stop rule is reached for its contact-residue branch.
**Limits.** Complex IV only, with one human structure for every species; mismatch counted against the consensus
residue, not against each species' own nuclear partner (the partners exist for 176 of these species); no result
of the cybrid calibration that RM11 proposed first is recorded in this repository.

## H10. Aligning by structure changes the columns of the mitochondrial proteins (RM13)

**Hypothesis.** The 13 proteins are membrane proteins of a few families of folds, and sequence aligners may
misplace columns that structure would put elsewhere.
**Test.** `mitoage download pdb 5XTD 5Z62 8H9S 5XTE 5XTH 5XTI 6ZQN 7AJD`, `structure within`, `structure
fold` (GPU), `structure msa`; into `within_species.tsv`, `chains.tsv`, `plddt.tsv` and `structural_msa.tsv`.
First run on the cluster on 2026-11-14 and 15 by another agent; rerun from this code, `within` and `msa` (on
that run's structures) give the same tables, timings aside; the folding was not rerun, and its 806 structures
are published (`mitoage download prepared --topic structures`). First, the 13 proteins from solved human structures (ND1-6 and ND4L from
5XTD, COX1-3 from 5Z62, CYTB from 5XTE, ATP6 and ATP8 from 8H9S; identity to MitoAge's human sequences 0.979
to 1.000), superimposed all against all with US-align. Then predicted structures of the 13 proteins in 62
mammals (806 structures, ESMFold2, checkpoint `biohub/ESMFold2-Fast`, one sample), aligned by FoldMason and
compared with MAFFT L-INS-i on the same sequences: the share of residue pairs placed in one column by both,
over 500 species pairs per gene, overall and in the conserved core.
**Result, folds.** Of the 78 pairs of mitochondrial proteins three exceed a TM-score of 0.5 normalised by the
longer chain: ND4-ND5 0.69, ND2-ND4 0.67, ND2-ND5 0.52 (340 to 447 aligned residues, RMSD 2.5 to 2.8 Å), the
antiporter-like triad of complex I; the next is 0.41 (COX1-ND4). The rest are unrelated folds. Scores
normalised by the shorter chain run high for short subunits and are not used.
**Result, alignments.** Median agreement 0.9997 over the 13 genes and 1.000 in the core; the lowest are ND6
(0.970; core 0.996) and ND5 (0.982), and no core is below 0.993 (ATP6). Mean pLDDT 84.8: 79 to 97 for ten
proteins, 56 to 62 for ND3, ND4L and ATP8. Disagreement does not follow model confidence (r = −0.03 across
genes).
**Cost.** Aligning is trivial (18 s FoldMason, 9 s MAFFT for the 13 genes); folding took 11,949 s for 806
structures on one A100 (40 GB), 14.8 s per structure, which extrapolates to about 120 GPU-hours on that card
for the 2,290 species of the build with all 13 proteins (a rough figure: proteins of other phyla differ in
length).
**Reading.** At the divergence of mammals structure does not move the columns: the sequence alignments stand,
and H0's soft genes (ND6, ND5) are also where structure and sequence differ most, by little.
**Limits.** Mammals only (62 species); one prediction per sequence; agreement between two methods is not
correctness. Distant phyla, where sequence aligners are more likely to err, are not assessed; at the same rate,
20 species across phyla would take about one GPU-hour.

## H11. The composition signal that survives the controls follows the slow life history, not lifespan alone

**Hypothesis.** In the site's correlation table, base composition is the one family of mtDNA features whose
correlation with lifespan survives phylogenetic contrasts and body mass: GC content and a T-to-C shift at
third codon positions. Once the contrast test is calibrated it survives in mammals only (the site
repository's `docs/PIC_CALIBRATION.md`). Two traits that go with lifespan have been tied to the same
composition:
- generation time, through the A>G signature on the heavy strand, which rises with generation time across
  mammals (Mikhailova et al. 2022, NAR 50:10264, doi:10.1093/nar/gkac779);
- gestation, which with mtDNA GC content was among the two most influential traits in a model of mammalian
  lifespan (Kulaga et al. 2021, IJMS 22:1073, doi:10.3390/ijms22031073).

The entry first asked "lifespan or generation time". That was the wrong question. Gestation, maturity,
birth weight, body size and lifespan evolve together as one slow or fast strategy, so a feature can be tied
to all of them for one evolutionary reason. The question here is how the composition signal is shared
across the slow life history.
**Test.** `mitoage longevity generation-time`, into `generation_time.<maturity>.tsv` and
`.classes.tsv`. The traits are log10 lifespan, log10 adult mass, and three traits from AnAge, each log10:
age at female maturity (the build's proxy for generation time), gestation (incubation in birds) and birth
weight. Each of the 43 features is standardised and regressed on the dated tree with Pagel's lambda, and
with lambda fixed at 1 as a check. Every term, mass included, is reported.
- For each trait T, on the species that have it: lifespan alone, T alone, each with mass, and lifespan, T
  and mass together.
- On the species with every trait: mass alone, lifespan alone, all five together, and two slow-fast axes.
  The axis "with size" is the first phylogenetic principal component (Revell 2009, correlation matrix) of
  the four traits and mass. The axis "beyond size" is the component of the four after mass is regressed out
  on the tree, fitted together with mass. Both are oriented so that lifespan loads positively and scaled to
  unit SD.
- Models without mass give the total association. Models with mass ask what a trait adds beyond size,
  which removes real signal if size is on the same path: size is part of the slow life history. Both are
  reported.
- Models with several traits show how a signal is shared. They are never read as one trait explaining it
  "instead of" another: the trait measured with less error takes the shared part, and maximum lifespan is
  a noisy record.

The 43 features were fixed before any result was seen:
- GC content, GC skew and AT skew of the protein-coding genes and of the whole record;
- A, C, G and T at fourfold sites (`fourfold`, gene sense);
- `c3_share`, the share of pyrimidine-ending codons that end in C;
- the 32 codons ending in U or C, each as a percentage of all complete codons.

q is Benjamini-Hochberg over the 43 features, within class, trait set, model and term.

**Run.** The numbers below are from build 2 with the EURING birds: longevity-genie/mitoage revision f75bdc9,
2,344 species. They come from Slurm job 6489877 on haumea, the "after" run of `slurm/manuscript-rerun.sbatch`
from commit 2ee27a1, with tables in `results/euring/`. The build before EURING (c525adb) was run twice:
- as job 6487335 from 858974b (tables in `docs/analyses/h11/`). All 14 tests pass in the job, and a second
  run reproduces both maturity modes byte for byte;
- as the "before" run, job 6489876.
`results/euring/comparison.txt` compares the two builds. The mammalian tables are identical: EURING adds only
bird tips. In birds only the maturity set changes, from 281 to 282 species. In an earlier version of this
test, the site's TypeScript contrast code, an independent implementation, gave the same t statistics as the
lambda-1 fits, to within 0.017 over 152 fits.

**Result, mammals.** Every slow trait follows the composition features. The counts are features of 43 at
q < 0.05:

| Trait set (species) | Lifespan alone | Trait alone | Lifespan + mass | Trait + mass | Lifespan, trait and mass: lifespan / trait |
|---|--:|--:|--:|--:|--:|
| Maturity (527) | 20 | 12 | 11 | 1 | 4 / 1 |
| Gestation (557) | 9 | 21 | 6 | 17 | 1 / 15 |
| Birth weight (470) | 18 | 23 | 5 | 19 | 0 / 16 |

On the 439 species with every trait:
- **Slow-fast axis beyond size: 32 of 43**, the strongest of all. All four traits load positively
  (lifespan +0.51, maturity +0.44, gestation +0.56, birth weight +0.48), 39% of the variance beyond size.
- **Slow-fast axis with size: 27 of 43.** All five traits load positively, 55% of the variance.
- Lifespan alone: 24 of 43. Mass alone: 3 of 43.

Examples, in SD of the feature per SD of the axis:

| Feature | Axis beyond size | Axis with size | Lifespan alone | Mass alone |
|---|---|---|---|---|
| GC content, whole record | +0.28 (q = 8 x 10^-5) | +0.34 (q = 2 x 10^-4) | +0.61 per tenfold (q = 0.002) | +0.12 per tenfold (q = 0.06) |
| `c3_share` | +0.33 (q = 6 x 10^-5) | | | |
| T at fourfold sites | -0.30 (q = 6 x 10^-5) | | | |

The direction is the one the A>G signature on the heavy strand predicts: in gene sense, less T and more C
along the slow end. Body mass by itself barely follows composition, so correcting for it changes little.
Lambda is 0.94 to 1.00, and with lambda fixed at 1 the axes give 31 and 26.
**Result, birds.**
- Lifespan follows nothing in any model, except codon UAU in the 115 species with a birth weight.
- Gestation, birth weight, mass and both slow-fast axes follow nothing (107 to 175 species).
- Only age at maturity follows composition, in 282 species: 7 features alone, 8 with mass, and 8 with
  lifespan and mass. The direction is the opposite of the mammalian slow end: later-maturing birds carry
  more T and less C in gene sense. With mass, in SD per tenfold maturity:
  - AT skew of the protein-coding genes, -0.77 (q = 0.013);
  - codon ACC, -0.71;
  - codon GCU, +0.94;
  - AT skew of the whole record, -0.67;
  - GC skew of the protein-coding genes, +0.55;
  - T at fourfold sites, +0.60 (q = 0.032);
  - GC skew of the whole record, +0.51;
  - codon CUU, +0.63 (q = 0.047).
- On the build before EURING (281 species), the same signs held, with 7, 6 and 7 features; CUU and the
  whole-record GC skew were then just above the threshold.
- Gene sense does not depend on which strand RefSeq deposited, so the reversal is not a strand artefact.
- Beyond size, the four bird traits do not share one axis: lifespan loads against the other three.
  Including size they do (all positive, 53% of the variance), and that axis follows no feature.
- With the male age filling a missing female one (`--maturity any`), mammals are unchanged (no mammal
  gains a value). Birds have 7 features with mass and 7 with lifespan and mass (283 species; 5 each before
  EURING), with the same signs.

**Reading.**
- **Mammals:** mtDNA base composition follows the slow life history as a whole. Taken together, the slow
  traits carry the signal more strongly than lifespan alone (32 against 24 features on the same species).
  Separately, gestation and birth weight carry it about as strongly as lifespan does. What each trait adds
  once the others are in the model is small and goes to whichever is measured best. The effect is modest:
  about 0.3 SD of each composition feature per SD of the slow-fast axis.
- **Birds:** neither lifespan nor the slow-fast axis follows composition. Maturity does, with the opposite
  sign, so the composition signal is not one relationship with life history shared across these classes.
**Limits.**
- Age at maturity is not generation length as Mikhailova et al. used it (Pacifici et al. 2013), which
  includes part of the reproductive lifespan; generation length is not in the build.
- Bird traits are coarse and sparse. Of the 292 live birds with a female maturity age, 71.9% are whole
  years and 129 are exactly 365 days. Only 107 birds have all four traits. (An earlier count, 296, included
  retired records.)
- Missing traits cost species: in mammals, 527 of 604 have a maturity age, 557 a gestation, 470 a birth
  weight and 439 all three.
- The bird maturity reversal comes from a screen of 43 features, with no replication.

## H12. The structural findings in every class: complex IV contacts, hydropathy, motifs (additions to H8-H10)

**Question.** H8 to H10 found nothing about lifespan, but H9 asked it only of all species pooled and of
mammals, and the protein features the build publishes (GRAVY hydropathy) were never put to this repository's
model at all. Does any structural feature follow lifespan in mammals or in birds (and in ray-finned fishes and
reptiles where the data allow)? Kept apart from H11, which is about base composition and life history.
**Fixed before any result was seen** (written into the session's notes before the first run, 2026-11-15; of
the new tests, nothing had been run):
- Species, per class: H1's default dataset (`aa_pgls.select`): not retired, lifespan record of high or
  acceptable quality and not estimated, body mass, all 13 proteins, fourfold counts, a tip of the dated tree;
  the class's part of the tree (`aa_pgls.class_covariance`). Under 25 species: not assessed.
- Model: log10 maximum lifespan on the feature, standardised over the species of the fit (a slope is log10
  years per standard deviation), with log10 adult body mass and the fourfold composition (A4, C4, T4), Pagel's
  lambda by maximum likelihood for each fit (`aa_pgls.fit`). Also, for comparison with H9's original model,
  with mass only. q: Benjamini-Hochberg within class, family and model; mass with composition is the primary
  model.
- Family "contacts" (H9 taken into each class): per species, the share of the 298 contact columns of COX1-3
  (5Z62, `mitoport burial`) whose residue differs from the column's modal residue among that class's tested
  species, the same over the 298 burial-matched controls, and their difference, which is the test. Three
  variables per class.
- Family "hydropathy": GRAVY of each of the 13 proteins and of all 13 together, from
  `protein_gene_features.parquet`. Fourteen features per class.
- H8's motifs: `longevity motif-pgls` rerun unchanged (lambda estimated on its model without the feature) and
  read class by class.
**Test.** `mitoage longevity structure-pgls` (after `download pdb 5Z62`, `mitoport contacts`, `mitoport
burial`, `longevity controls`, `fourfold`, `dated-tree`), into `structure_pgls.tsv`; on haumea,
`slurm/structure-pgls.sbatch`, which also runs `motif-pgls`, the tests, a check that each test fails when the
code it guards is broken, and a rerun. The first round's build is longevity-genie/mitoage at c525adb; the numbers
below are from the rerun on f75bdc9 (see "Which build"), with the first round's in parentheses where they differ.
**Run.** Slurm job 6487327 on haumea, 2026-11-15, from commit 858974b; tables and log in `results/h12/`. The
26 tests of `test_structure_pgls.py` and `test_motif_pgls.py` pass, each of 11 deliberate breaks of
`structure_pgls.py` makes a test fail, and a second run of `structure-pgls` gives the same bytes. The contact map
(`sites.tsv`) is byte-identical to H5's cluster run. `motif_pgls.tsv` is byte-identical to the table H8 reports.
`composition.tsv` differs from the laptop's only in the last bit of log10 lifespan and mass (at most 9 x 10^-16).
**Species.** Mammals 604, birds 381 (379), ray-finned fishes 235 (most fishes lack a body mass), reptiles 97.
**Result, contact mismatch (H9 in every class): no class on the first round's alignment; reptiles, narrowly, on
the realigned build.** Slopes are log10 years per standard deviation of
mismatch, primary model (mass and composition):

| Class | Contact minus control | 95% interval | q | Contact alone | q |
|---|---|---|---|---|---|
| Mammals (604) | -0.007 | [-0.029, +0.015] | 0.83 | +0.006 | 0.83 |
| Birds (381) | -0.015 | [-0.037, +0.008] | 0.20 | -0.045 (p = 0.02) | 0.059 |
| Ray-finned fishes (235) | -0.015 | [-0.071, +0.041] | 0.84 | -0.006 | 0.84 |
| Reptiles (97) | -0.111 (p = 0.031) | [-0.211, -0.010] | **0.046** | -0.114 (p = 0.03) | **0.046** |

First round (c525adb): birds (379) -0.014 [-0.037, +0.009], q = 0.24, contact alone -0.044, q = 0.065; reptiles
-0.107 (p = 0.035) [-0.207, -0.008], q = 0.052, contact alone -0.112, q = 0.052. Mammals and fishes are identical.
**The reptile result crosses q = 0.05 with no reptile added**: the same 97 species, but the contact columns and the
class's most common residues come from the alignment of all species, and the nine birds moved it. A result that
passes or fails with the alignment of another class is not read as a finding.

The mammalian interval is H9's again. Birds bound the contact-specific effect to at most about 8% shorter or 2%
longer maximum lifespan per standard deviation. Contact mismatch alone leans negative in birds and reptiles,
but so does mismatch at the burial-matched controls (birds -0.032, reptiles -0.082). That points to overall
divergence from the class consensus, not the contacts (reptile controls -0.082 on both builds). The reptile class is
diverse (mean mismatch to its consensus 0.23 at contacts, against 0.06 to 0.09 elsewhere), and 97 species is
the smallest set. With mass only, nothing passes (smallest q 0.051, reptiles; 0.056).
**Result, hydropathy (GRAVY): mammals yes, two of 14 genes; birds no.** Mammals, primary model: GRAVY of ND5
-0.026 per SD ([-0.042, -0.011], p = 0.001, q = 0.013) and of ND6 +0.030 ([+0.010, +0.050], p = 0.003,
q = 0.024). That is about 6% shorter and 7% longer maximum lifespan per standard deviation; the other 12 genes
and the total have q >= 0.07 (total -0.029, p = 0.016). With mass only, five pass: ND5 (q = 0.003), the total
of the 13 (-0.038, q = 0.003), ND6, ND1 and ND2 (q 0.005 to 0.023). Fourfold composition takes about a quarter
of the total's slope and the significance of ND1, ND2 and the total. Birds, primary model: nothing (smallest
COX2 +0.028, p = 0.007, q = 0.10; +0.029, q = 0.079); with mass only COX2 passes alone (q = 0.037; 0.027), so it
is not counted. Fishes
and reptiles: nothing (smallest q 0.23 and 0.66).
**Result, H8 motifs class by class (rerun, unchanged): no class.** Smallest q: mammals 0.16
(G-quadruplexes on the reverse strand, +0.017 per SD, p = 0.005, 604 species), birds 0.41 (forward-strand
G-quadruplexes over the major arc, p = 0.023, 334 of the 381; 0.42 and 332 of 379), fishes 0.21, reptiles 0.16; 0 of 132 at q < 0.05.
**Reading.** No structural finding follows lifespan in birds. In reptiles, contact mismatch passes on one alignment and
not on the other, with the same species, so it is not read as a finding. In mammals, mismatch at complex IV's contacts
and the motif features show nothing. The two genes whose hydropathy passes go in opposite directions: ND5,
like the total, is less hydrophobic in longer-lived mammals, and ND6, the only protein read from the other
strand, is more so. The 12-gene fourfold covariate does not cover ND6, so a strand-wide mutational bias, which
would push a gene read from the other strand the other way, was the first explanation. The second round tested
it and it does not hold: with each gene's own fourfold composition in the model both slopes stay (below), so
codon composition does not explain them. It is two of 14 features in one class, with no replication, and it is
not shared with birds.
**Seen before the run.** While the tests were being written, a failing assertion printed three q values of the
bird contact family (0.136, 0.165, 0.360). The run shows they were the mass-only model's (an earlier version of
this entry wrongly named the primary model). The features and the model above were written down before that and
were not changed after it.
**Not assessed in the first round, and why.**
- H10's predicted structures exist for 62 mammals only, and H10 is a check of the alignments with no lifespan
  hypothesis attached. Birds would need about 379 x 13 foldings, some 20 GPU-hours on an A100 at H10's rate.
  Not run, in either round.
- The per-residue ESM-2 embeddings cover 184 mammals and no bird. Assessed in the second round, with new
  embeddings of every species tested.
- H9's partner-specific branch needs nuclear partners, which existed for 176 mammals only. Assessed in the second
  round, with vertebrate ortholog sets.
**Second round (asked for on 2026-11-15: run what was not assessed unless it takes many hours).** Folding birds
(about 20 GPU-hours, and no lifespan hypothesis) stays not assessed. Three tests were added, each fixed before any of
its results, with H12's species, model and correction:
- *Hydropathy with each gene's own composition.* A third model for the 14 GRAVY features, "mass+own composition":
  log10 mass and the shares of A, C and U at the fourfold sites of the same region (gene sense) in place of the
  12-gene composition, which leaves out ND6. Reading fixed in advance: if ND6's (and ND5's) GRAVY slope loses
  p < 0.05 under its own composition, its signal is read as composition; if q < 0.05 stays, it is not explained by
  its own fourfold composition.
- *Discordance of the mito-nuclear contact pairs* (`mitoport partners`). NCBI's vertebrate ortholog sets of the
  eleven nuclear subunits, one protein per species, MAFFT L-INS-i on one thread. A pair is a COX1-3 contact residue
  and a nuclear residue within 4.5 A of it in 5Z62; its control keeps the nuclear residue and takes the mito
  residue's burial-matched control. In a species a pair is discordant when exactly one side differs from its
  column's consensus among the class's species. Variables: discordance over the contact pairs, over the control
  pairs, and their difference (the test); a species needs half the contact pairs filled on both sides.
- *ESM-2 divergence at the contacts* (`structure embed` on a GPU, then `structure-pgls --family embedding`). ESM-2
  650M per residue on the ungapped COX1-3 of every species H12 tests; per species and gene the mean over the contact
  residues and over the control residues; the cosine distance of each from the class's mean vector, weighted over
  the three genes. Variables: contact, control, contact minus control (the test).
Jobs: `slurm/h12-partners.sbatch` (CPU) and `slurm/h12-embed.sbatch` (GPU, after the first).
**Second round, run.** On haumea, 2026-11-15: CPU job 6488097 and GPU job 6488098, from commit 2933535; tables and
logs in `results/h12/round2/`. Two earlier CPU runs stopped on their own checks. The first (6488055) refused to
rebuild the control pairs from `sites.tsv`, whose rounded accessibilities give a different COX1 matching, so
`mitoport burial` now records the pairing (`control_pairs.tsv`), with `sites.tsv` byte-identical. The second
(6488075) found one of 19 code breaks uncaught (species filling under half the contact pairs scored), and a test now
covers it. In the runs reported: 39 and 20 tests pass, 19 of 19 and 2 of 2 code breaks are caught, and every result
table is byte-identical on a second run. The rows the first round also fitted are identical to it. The ESM-2
checkpoint is facebook/esm2_t33_650M_UR50D at commit 08e4846 (the one snapshot in the cluster's cache; the job
failed to record it itself, which `structure embed` now does from the cache's file paths). The 3,945 sequences
(1,315 species x 3 genes) took 260 s on one A100.
**Second round, hydropathy with each gene's own composition: the mammalian result stands.** Mammals (604): ND5
-0.027 per SD ([-0.042, -0.012], p = 5 x 10^-4, q = 0.007) and ND6 +0.033 ([+0.013, +0.053], p = 0.002,
q = 0.011), almost exactly the slopes under the 12-gene composition (-0.026, +0.030) and with mass only (-0.028,
+0.034). Under the reading fixed in advance, neither signal is explained by its own gene's fourfold composition.
The total of the 13 (-0.029, q = 0.07) and ND1, ND2 and COX3 (q 0.076) stay below the threshold. Birds (379):
COX2 +0.030 ([+0.010, +0.050], q = 0.041). It passes under this model and with mass only (q = 0.027), but not
under the primary model (q = 0.079): one gene of 14, borderline, and in mammals COX2 shows nothing (p = 0.76). On
the build with the EURING birds (381, job 6489877) it passes under none but mass only: +0.029, q = 0.056 with its
own composition, q = 0.037 with mass only.
Fishes and reptiles: nothing (smallest q 0.32 and 0.57).
**Second round, discordance of the mito-nuclear contact pairs: no class.** 687 contact pairs (COX1-3 with all
eleven nuclear subunits). Species with a partner sequence and half the contact pairs filled: mammals 177, birds 78,
ray-finned fishes 41, reptiles 14 (not assessed). Mean discordance at contact and control pairs: mammals 0.135 and
0.135, birds 0.101 and 0.108, fishes 0.116 and 0.102. Contact minus control, primary model: mammals -0.021 per SD
([-0.056, +0.014], q = 0.34), birds -0.006 ([-0.056, +0.043], q = 0.80), fishes +0.026 ([-0.061, +0.114],
q = 0.55). Smallest q in any class and model 0.077 (fishes, contact alone, mass only).
**Second round, ESM-2 divergence at the contacts: no class.** Contact minus control, primary model: mammals (604)
+0.005 per SD ([-0.032, +0.042], q = 0.79), birds (379) -0.007 ([-0.026, +0.013], q = 0.50), fishes (235) -0.014
(q = 0.42), reptiles (97) -0.001 (q = 0.98); smallest q in any class and model 0.27. The distances themselves are
small (mean cosine distance to the class mean 0.0005 to 0.0026): mean ESM-2 vectors over a hundred residues differ
little between species of one class, and standardising the variable does not create variation the model never
sees.
**Second round, reading.** Mito-nuclear compatibility at complex IV's contacts, measured three ways (mito-side
mismatch, discordance of the two sides, embedding divergence), follows lifespan in no class, with intervals that
exclude effects larger than about 12% of maximum lifespan per standard deviation in mammals and birds. The contact pairs
are not more concordant than the control pairs in mammals, and only slightly so in birds (not tested). The one
structural feature that does follow lifespan is the hydropathy of ND5 and ND6 in mammals. It is not codon
composition, neither that of the 12 genes nor the gene's own. With opposite signs in the two genes and nothing
comparable in birds, what it means is open.

## H13. Audit of the build's records: strand, gene order, genetic code, duplicates and mislabels (an audit, not a hypothesis test)

**Why.** A referee of the manuscript asked for checks of the records themselves: mislabelled species, duplicate
genomes, misassigned strands. The site's roadmap item RM8 ("Detect possible duplicate species") is part of it. This
entry tests nothing about longevity. It describes every mitogenome record of the build and names candidates for a
person to check; no candidate is a verdict, and nothing is removed from the build because of it.
**Fixed before any result was seen** (written and committed with the code, 2026-11-15, before any audit step had run
on the build; the tests had run on a few records):
- *Records.* Every species of the build longevity-genie/mitoage at c525adb (H12's revision, 2,335 species), with
  its RefSeq record read from RefSeq release 237, the release build 2 read (`download refseq --from-snapshot`, md5
  checked against the pin; the CONTIG-only records' sequence from `download contig --from-snapshot`). Features are
  named with build 2's rules (`build2/genome.py`). A species whose record is not in the release (the retired ones) is
  a row in every table, marked unknown with the reason. Every input's SHA-256 is written to `inputs.tsv`.
  `mitoage audit records`.
- *Strand* (`audit strand`, `strand.tsv`, `strand_summary.tsv`). The vertebrate convention: 12 protein genes on one
  strand, ND6 on the other. Every annotated copy of the 13 is read. Conventional: every copy on its conventional
  strand with the 12 on the record's + strand; reversed: every copy on the other strand; mixed: copies both ways (the
  genes against the convention are named, and which way the majority goes); unknown: no record, outside Vertebrata
  (by the record's lineage; the convention is the vertebrates'), no protein gene annotated, or fewer than 7 of the 13
  annotated. Beside it, the record's GC and AT skew over the whole sequence as the build gives them
  (`region_stats`, total_mtDNA, + strand as deposited), and the build's own `orientation` column; the summary counts,
  per class and orientation, the records with a negative and a positive skew. The question it answers is whether the
  sign of the skew follows the orientation of the deposit.
- *Gene order* (`audit order`, `order.tsv`, `order_counts.tsv`, `order_arrangements.tsv`, `order_flags.tsv`). The 13
  protein genes, 2 rRNAs and tRNAs as named by build 2, every copy, signed by strand, in order of start on the record;
  the control region left out. Canonical by rotation: started at COX1 when the record has one COX1, else the first
  rotation in string order. A second form, canonical under rotation and reversal, tells a reversed deposit from a
  rearrangement. Two tiers: "full" (with tRNAs; unknown when no tRNA is annotated) and "genes" (13 genes and 2
  rRNAs). Arrangements counted per class, order and family (the build's taxonomy). A record is flagged when its
  family has at least 3 records with a known arrangement at that tier, its arrangement occurs once there, and all the
  other records of the family share one arrangement. Each flag says what differs, checked in this order: reversed
  deposit (the family's arrangement read from the other strand), gene content (genes absent or extra), strand only
  (same order, some genes on the other strand), tRNAs only (full tier: the 15 genes in the family's order), order.
  Whether a difference is a real rearrangement or an annotation error is not decided: every flag but a reversed
  deposit says "cannot tell from the record".
- *Genetic code* (`audit code`, `code.tsv`, `code_genes.tsv`). Declared: the /transl_table most of the record's
  protein genes state (the lower number on a tie; none stated is unknown). Expected: the mitochondrial genetic code
  NCBI Taxonomy (the dump of 2026-11-13, `download taxonomy --from-snapshot`) gives the record's taxon, made current
  through merged.dmp. Internal stops: each copy of the 13 genes in its own sense from codon_start, complete codons,
  every codon but the last complete one; codons with an ambiguity code and codons a /transl_except gives an amino
  acid are not stops. Under the declared and the expected code. A record is flagged ("neither clean") when both codes
  are known and both give at least one internal stop; it is not flagged when either is clean; otherwise unknown. A
  code the audit does not know (it knows NCBI tables 1, 2, 4, 5, 9, 13, 14, 21, 24 and 33) is unknown, never clean.
- *Near-identical pairs* (`audit duplicates`, `pairs.tsv`). Nucleotide identity over the 13 protein genes, from the
  build's codon alignments concatenated: matching bases over the columns where both records have A, C, G or T, at
  least 5,000 such columns. Threshold **0.99**, fixed now and why: mitochondrial divergence within an animal species
  is usually under 1 %, while sister species usually differ by more than 2 % at COX1 (the barcode gap; Hebert et al.
  2003, Proc R Soc B 270:313-321). Two species' records within 1 % over all 13 genes are where one species would be:
  one animal under two names, a misidentified specimen, a hybrid's or introgressed mitogenome, or a very young pair of
  species. Pairs at 0.999 or more are also marked near-identical; pairs whose whole sequences are equal are reported
  whatever the alignment says. Each pair carries the rank at which its taxonomy first agrees (same genus, same family,
  ... different class); pairs from different genera are read first.
- *Out-of-place records* (`audit duplicates`, `nearest.tsv`). Chosen in advance: amino-acid identity over the 13
  proteins concatenated (the build's protein alignments; the 20 amino acids, at least 1,500 shared residues), not
  COX1 alone: more sites, and less saturated than nucleotides at the depth of orders and classes. A record is out of
  its order when some record outside its order is strictly closer to it than every other record of its order; the
  same for class. A record whose order (class) has no other record compared with it is unknown. Candidates are ranked
  by the margin (identity to the closest outside minus identity to the closest inside). The build's classes include
  paraphyletic ones (Reptilia without birds, Actinopterygii beside Sarcopterygii), so some out-of-class candidates
  are expected from the scheme alone.
**Checks of the code.** `tests/test_audit.py` on real records of release 237: a record reverse-complemented by the
test swaps conventional and reversed and the sign of its GC skew; a record started elsewhere on its circle keeps its
arrangement; a bovid edited by the test (reverse-complemented, one gene on the other strand, ND6 removed, ND1 and ND2
swapped, two tRNAs swapped) is flagged in a family of three bovids with the kind the edit made, and nothing is flagged
without the edit; our translation under the declared code equals RefSeq's /translation; under the standard code every
internal TGA is a stop; a stop written into COX1 is counted once and survives reverse complement; an unknown code and
a record without code or genes are unknown; two identical records pair at identity 1 and records relabelled into
another order sit outside it. `slurm/job3-audit-mutations.py` breaks 23 guards one at a time and checks that a test
fails each time; `slurm/job3-audit.sbatch` runs every step twice and compares the bytes.
**First run, stopped.** Slurm job 6490037 on haumea, 2026-11-15, from commit 1d5f8aa. Every step ran and the 27 tests
passed, but the job stopped on its mutation check: 21 of 22 breaks were caught, and ignoring /codon_start was not (a
test now moves COX1's start one base back with codon_start 2 and asks for the same codons). Its tables were read
before the second run, and they showed one departure from the rules above: the build writes 44 missing orders and 6
missing families as an empty string, which the code read as an order (family) named "", so records with no order were
compared with each other as if they shared one. The rules above already say a record without an order is unknown; the
code now reads an empty label as missing, and a test checks it. Nothing else was changed after that run.
**Run.** Not yet run.

## H14. In complexes I, III and V too, the coupling belongs to whole proteins, not to the residues that touch (H5's open item)

**Question.** H5 found, for complex IV, that its nuclear subunits co-evolve with the mitochondrial COX genes as whole
proteins, while the residues that touch carry nothing extra: on the mitochondrial side contact residues couple with
the nuclear partners no more than burial-matched non-contact residues, on the nuclear side interface columns no more
than the rest, and (H12, second round) contact pairs are no more concordant than control pairs. Does that hold in
complex I (5XTD; ND1-6, ND4L), complex III (5XTE; CYTB) and complex V (8H9S; ATP6, ATP8)? This is a question about
mito-nuclear co-evolution in general, not about longevity.
**Fixed before any result was seen.** This entry was committed with the code, before the job that computes the
readouts ran. Complex IV goes through the same code and readouts again, so that it can be compared with the others.
- *Contact maps* (`mitoport contacts --complex I|III|V`, `mitoport burial --complex ...`): H5's definitions unchanged.
  A contact is a residue of a mitochondrial subunit with a heavy atom within 4.5 Å of a nuclear subunit. Controls are
  matched on accessibility in the isolated chain, within the gene, without replacement. Chains are named by their
  UniProt references (`mitoport/complexes.py`, checked by hand against each entry's entity names and chain ids). Each
  mitochondrial gene is read from one reference chain. Complex III is a dimer, so cytochrome b is read from chain J
  and the other copy (V) counts as a mitochondrial partner, not a nuclear one. A nuclear gene in several chains (the
  c-ring, alpha and beta, the dimer's copies) takes the union of its chains' interface residues.
- *What the structures give* (seen while the maps were built; they are the design, not results): complex I has 621
  contact residues (ND1 116, ND2 105, ND3 32, ND4 144, ND4L 17, ND5 176, ND6 31), and 27 of its 37 nuclear subunits
  touch a mitochondrial one. Complex III has 104 (CYTB), with 6 of 9 subunits touching. Complex V has 129 (ATP6 95,
  ATP8 34), with 6 of 15 touching. ATP8 runs out of candidate controls: 17 of its 34 contacts have none, and pair tests
  leave them out. All contact and control residues of the four complexes map to an alignment column. The subunits with
  at least 5 interface residues in the structure number 24 in complex I, 5 in complex III and 5 in complex V. Complex
  IV's contact map from the generalised code is byte-identical to H5's and H12's (`sites.tsv` SHA-256 6f219d3a...).
- *Seen before this entry was written, and therefore not blind*: `mitoport burial` prints the median identity to
  human of the contact, control and other residues gene by gene. For complexes I, III and V that printout was read.
  Readout R1a (conservation) is that comparison pooled, so it is reported as descriptive only and enters no decision.
  No coupling, interface or pair statistic of complexes I, III or V had been computed.
- *Sequences.* For R0 to R2 (H5's `iface` design): NCBI's mammalian ortholog sets (`download orthologs`) of H5's 62
  genes, plus the 37 nuclear genes of 5XTD and the 15 of 8H9S; one protein per mammal (`mitoport select`); MAFFT
  L-INS-i on one thread (`iface-align --threads 1`). The sequence distance only: no per-residue embeddings exist for
  the new genes, so every `esm_` column is not assessed (`--skip-esm`). For R3 (H12's `partners` design): NCBI's
  vertebrate ortholog sets of the nuclear subunits that touch (complex I 27, III 6, IV 11, V 6); one protein per
  species of mammals, birds, ray-finned fishes and reptiles; L-INS-i on one thread. The species of each class are
  those of H12 (H1's default dataset with a partner sequence, half the contact pairs filled on both sides).
- *Readouts, per complex.*
  - R0, the positive control (is there whole-protein coupling to localise?). `seq_control_other` of `mitoport iface`:
    over pairs of mammals, with patristic distance partialled out, the correlation of a subunit's non-interface
    columns with the mitochondrial burial-matched controls, so that no contact residue sits on either side. The
    complex's nuclear subunits are compared with the 30 cytosolic controls, Mann-Whitney, two-sided. H5's own group
    contrast (`seq_contact_iface`, subunits with an interface) is reported beside it.
  - R1a, mitochondrial side, conservation: identity to human at the contact against the control residues, the
    complex's genes pooled (Mann-Whitney, as H5), and over the matched pairs (Wilcoxon signed-rank). Descriptive.
  - R1b, mitochondrial side, coupling: `seq_excess_mt` (the interface's coupling with the mitochondrial contacts
    minus with their controls) of the subunits with at least 5 interface columns, against the cytosolic controls
    (random column sets of the usual interface size, seed 0). Mann-Whitney, two-sided. The interface hypothesis
    predicts partners above the controls.
  - R2, nuclear side: `seq_excess_nuc` (interface columns minus the other columns, both against the mitochondrial
    contacts), with the same sets and test; the hypothesis predicts partners above the controls.
  - R3, pair concordance: per species, discordance at the contact pairs minus discordance at the control pairs, as
    H12 defined them. The mean over species is estimated as the intercept of a phylogenetic regression with Pagel's
    lambda by maximum likelihood on the class's part of the dated tree, with a 95% interval. Mammals are primary;
    birds, ray-finned fishes and reptiles are reported where they have 25 species or more. The hypothesis predicts a
    negative intercept (contact pairs more concordant).
  - Correction: Holm over the 12 interface tests (R1b, R2 and R3 in mammals, for the four complexes), at 0.05.
  - `mitoport readouts` writes all of them into `contact_readouts.tsv`.
- *Reading, per complex.*
  - If R0 shows no whole-protein coupling (`seq_control_other` of the subunits not above the controls at p < 0.05,
    uncorrected), the complex has no coupling to localise. The interface question is then not assessable there, and
    the complex is not counted for or against H5's reading.
  - Otherwise H5's reading holds in the complex when none of R1b, R2 and R3 (mammals) passes Holm in the predicted
    direction. It fails when at least one does: the contacts then carry coupling beyond the whole protein in that
    complex, and the entry says on which side.
  - A readout that passes Holm in the opposite direction is reported as such and does not count against H5's reading.
  - A readout that cannot be computed (fewer than 3 proteins on a side, fewer than 25 species) is not assessed, and is
    said as such; the complex's reading then rests on the readouts that remain, and the entry says so.
  - Power is not equal: complexes III and V have about 5 subunits with an interface each, so a null there bounds
    little. The medians and, for R3, the interval are reported so that the size can be read, not only the p value.
- *Secondary, not read for this question*: `mitoport iface` also fits H9's lifespan model on each complex's contact
  mismatch, and `mitoport partners` fits H12's lifespan models on the pair discordance. Their tables are written
  as always and reported with numbers only, with no claim.
**Test.** `slurm/job1-contacts.sbatch` on haumea (the downloads are made beforehand on the login node). It runs every
step for the four complexes, the tests, `slurm/job1-contacts-mutations.py` (each deliberate break of the new code must
make a test fail) and H12's mutation check, and then the iface, partners and readouts steps twice, to check that the
rerun gives the same bytes. It also runs complex IV through the code before the generalisation on the same inputs and
checks that the contact map, `coupling.tsv`, `lifespan.tsv`, `partners.tsv` and `partners_species.tsv` are the same
bytes. The build is longevity-genie/mitoage at
c525adb, and the dated tree is the site's phylogeny.json at 6fe8e72, as in H12.
## H15. The published mtDNA-encoded peptides across the whole build, against nulls of the host gene (RM16)

**Question.** H4's pilot read four peptides with a window null. RM16 (site repository) asks, in order: is each
published frame there in each species, under the species' own code and the standard code; is it kept beyond what its
host gene's conservation gives for free, against windows of the same gene and then against the host's synonymous
codons shuffled; does keeping or losing it follow the tree, and lifespan. RM16's fourth question, a blind scan of
every frame of 20 codons or more, is not part of this entry.
**Fixed before any result was seen** (written with the code, 2026-11-15; no table of this entry had been computed,
and the only conservation numbers seen were H4's):
- *Catalogue* (`mitoage peptides catalogue`, `peptides/catalogue.py`). Each sequence and placement copied from the
  primary paper, each PMID checked by its title on PubMed: humanin and SHLP1 to SHLP6 (Cobb et al. 2016, PMID
  27070352, Fig. 1A; humanin discovered by Hashimoto et al. 2001, PMID 11371646), MOTS-c (Lee et al. 2015, PMID
  25738459, Fig. 1A), MTALTND4 (Kienzle et al. 2023, PMID 37198654, Table 1), MTALTCO1 (Robitaille et al. 2026,
  PMID 41923604; its residues from Kienzle et al. 2023, Table 1, OpenProt IP_306389), MENTSH (Yen et al. 2026, PMID
  42559346, Fig. 1F). The human row (NC_012920.1) must read each published sequence exactly once under its published
  reading, followed by a stop, and the printed nucleotides of MOTS-c and MENTSH exactly; otherwise the step refuses.
  It does: MTALTND4 and MTALTCO1 start at the published rCRS positions, MENTSH's ATT holds m.4977, the SHLPs sit 2
  and humanin 1 below the paper's numbers. Readings: standard code for humanin, the SHLPs and MOTS-c (cytosolic, as
  published); the mitochondrial code with AGA/AGG read as arginine for MTALTND4 and MTALTCO1, as their authors read
  them; the mitochondrial code as NCBI has it for MENTSH, which ends at AGG. SHLP1 to SHLP5 lie on the strand
  opposite the 16S rRNA. **Left out** (`left_out.tsv`): MASL (D-loop: the build aligns no control region), SHMOOSE
  (starts over the serine tRNA, which the build does not align; PMID 36127429), Gau (antisense of COX1, mostly
  computational), the other candidates of Kienzle et al.'s Table 1, SHLP2b (a prediction).
- *Frame in a species* (`mitoage peptides frames`). The species' first base at or after the human start column, read
  on the peptide's strand until a stop. Kept: a start codon of the reading there (standard: ATG; own: the NCBI
  start codons of the species' table) and no stop in the first ceil(0.9 n) codons. Lost: an early stop, or open
  without a start. Not assessed, a third answer never counted as either: under half the human window's columns hold
  a base, a base other than ACGT before the decision, or the gene ends first. Both readings for every species: its
  own NCBI table (with AGR as arginine for the two peptides read that way), and the standard code. Species: every
  row of the build that is not retired. A class is summarised with at least 20 assessed species.
- *Window null* (H4's, made exact). Every window of n codons in the same gene, same strand and, in a protein gene,
  the same frame, not overlapping the peptide in its own frame, that is an ORF in human under the same reading
  (statistic "kept"), or only stop-free in human (statistic "open", H4's). Statistic: the share of the class'
  assessed species where the window is kept (open). p = (1 + windows at or above) / (1 + windows); under 20 windows
  is not assessed.
- *Synonymous null* (the three frames inside protein genes). In each species, the host gene's codons are permuted
  among positions coding the same amino acid under the species' own code (protein and codon counts unchanged), 1,000
  times; the frame is read at the species' window in each shuffle. Statistic: the share of species kept; p = (1 +
  shuffles at or above) / 1,001.
- *Dinucleotide null* (the eight frames inside rRNA genes, which have no synonymous codons). Each species' window
  (3n bases from the start) shuffled 1,000 times keeping its dinucleotide counts and its first and last base
  (Altschul-Erickson, by a random Eulerian path). Same statistic and p. **Its weakness**: it holds composition, not
  the rRNA's structure, so a frame more open than its composition predicts may be kept by the rRNA's own
  constraint. Only the window null can say "beyond the host's conservation" for an rRNA frame; the dinucleotide
  null is reported as context and decides nothing.
- *Relatedness.* The shuffles draw each species independently while relatives inherit their frames, so at the level
  of species they overstate the evidence. Both shuffle nulls are also computed with one species per family (the
  family's assessed species with the smallest SHA-256 of its slug), and the family level is the one read.
- *Tree* (`mitoage peptides tree`), on the build's dated tree (the same as H12's: site repository at 6fe8e72).
  Per peptide, reading and class (and all species): the minimum number of changes between kept and lost (Sankoff,
  unit costs; polytomies exact; species not assessed are missing tips), losses and gains on one most parsimonious
  reconstruction (root kept on a tie, changes as close to the tips as allowed), and clustering: the changes against
  1,000 random reassignments of the same states to the same tips, p_clustered = (1 + at most as many) / 1,001.
- *Lifespan* (`mitoage peptides pgls`, secondary): H12's model unchanged (`aa_pgls.select`, `class_covariance`,
  `fit`) in mammals, birds, ray-finned fishes and reptiles: log10 lifespan on the frame kept (1) or lost (0), with
  log10 body mass and A4, C4, T4, Pagel's lambda; the published reading; not assessed species left out; at least 25
  species with 5 kept and 5 lost, or not assessed. The state is binary and not standardised (H12 standardised its
  features): the slope is log10 years kept minus lost.
- *Thresholds and reading.* q is Benjamini-Hochberg over peptides and classes within null, statistic, level and
  reading (nulls), over peptides within class (lifespan). A frame is called **kept beyond its host's conservation**
  in a class only if, under its published reading, the window null (statistic kept) gives q < 0.05 and, for a
  protein-gene frame, the synonymous null at the family level gives q < 0.05 too. Loss is called clustered at
  p_clustered < 0.05. Lifespan: q < 0.05. Seeds: 20261115 for every shuffle and permutation, each species' draws
  seeded from (seed, peptide, species) by SHA-256.
**Test.** `mitoage peptides catalogue`, `frames`, `nulls --jobs 8`, `tree`, `pgls` (after `longevity dated-tree`,
`controls`, `fourfold`), into `data/output/peptides/`; on haumea `slurm/job4-peptides.sbatch`, which also runs the
tests, a check that each test fails when the code it guards is broken (`slurm/job4-peptides-mutations.py`) and a
second run of every step. Build: longevity-genie/mitoage at c525adb, as H12.


## H16. The published cybrids, simulated: does mismatch at the mito-nuclear contacts predict which complex fails?

**Question.** H5, H9 and H14 ask whether the residues where mtDNA-encoded subunits touch nuclear ones carry
anything beyond the whole protein, using sequences alone. Xenomitochondrial cybrids measure it: a donor species'
mtDNA in cells of a host species, so every nuclear subunit is the host's and every mtDNA-encoded one the donor's.
Before any of our interface scores is read as incompatibility, it should order the published cybrids as their
measured complex activities do. This is RM11's calibration, first proposed and not run until now.
**Data.** `src/mitoage_pipelines/mitoport/cybrid_measurements.tsv`, transcribed by hand with the figure or table of
each value: human cells with chimpanzee, bonobo and gorilla mtDNA (doi 10.1074/jbc.273.23.14210); mouse cells with
*Mus spretus* and rat mtDNA (10.1074/jbc.M004070200), rat again from a second laboratory (10.1074/jbc.M004053200),
and *M. caroli*, *M. terricolor* (the papers' *Mus dunni* line), *M. pahari* and *Otomys irroratus*
(10.1093/molbev/msg132); donors that gave no respiring cybrid: orangutan in human cells (10.1073/pnas.94.17.9131;
the paper names no species, so both *Pongo* species are scored) and two hamsters in mouse cells
(10.1074/jbc.M004070200). In all, 30 measured (host, donor, complex) points: 9 reduced, 20 not, and one where the
two laboratories disagree (rat, complex III: 37% of control in one, not reduced in the other), which is left out.
Only 13 points have a stated value; the rest are given as not significantly different, without a number, and are
used only as not reduced. Complex III in the murid papers is a linked II+III assay, which shows only severe defects.
**Method.** `mitoport cybrid`, fixed in its docstring before the first run. Contact maps of the four human complexes
(I 5XTD, III 5XTE, IV 5Z62, V 8H9S) with H14's chain table and definition (4.5 A heavy atoms, the reference chain of
each gene, complex III's other cytochrome b a mitochondrial partner): 621, 104, 298 and 129 contact positions, the
same sets `mitoport contacts` finds (tested). Host and donor proteins from RefSeq at fixed versions, placed on human
positions by a global alignment; the cytochrome b the 2003 paper sequenced from its own cell lines replaces RefSeq's
for *M. spretus*, *M. terricolor* and *Otomys* (identical for the first two, 2 residues apart for *Otomys*). Scores per
(host, donor, complex): differing contact positions (`d_contact`), their share (`f_contact`), those with a negative
BLOSUM62 score, and the same counts at the resolved non-contact positions and over every position (`d_all`, `f_all`:
divergence, the baseline any interface score has to beat). Tests: T1 Spearman of the stated activities against each
score; T2 the area under the ROC curve for reduced against not reduced; T3 whether the highest-scoring complex of a
pair is its most reduced one; T4 whether donors that gave no cybrid score above every donor that did, in the same host.
No significance is claimed anywhere: the points are few and not independent (one host, one laboratory per series).
**Result** (`results/h16/summary.tsv`, `calibration.tsv`).
- T1, 13 stated values: the contact count orders them worse than divergence (Spearman −0.20 against −0.36 for
  `d_all`; as shares −0.46 against −0.59).
- T2, 9 reduced against 20 not: `d_contact` 0.80, `d_all` 0.74, non-contact positions 0.73; as shares 0.68, 0.61 and
  0.54. The edge of the contact count over divergence is about 10 of the 180 reduced-against-not orderings.
- T3, six pairs with a reduced complex: the raw counts, contact or not, always name complex I, the largest, and so are
  right for the three apes and the rat and wrong for *M. pahari* (complex IV reduced) and *Otomys* (complex III at 2%);
  shares do worse (2 of 6 for `f_contact`, 1 of 6 for `f_all`). No score here says which complex fails.
- T4: both orangutans and both hamsters score above every donor that respired, with the contact count and with
  divergence alike.
- Exploratory, added after the above had been read (T2b, T2 within one complex): complex IV 0.86 for the contact count
  against 0.72 for divergence (3 reduced, 6 not); complex III 0.86 against 0.71 (one reduced); complex I 0.60 against
  0.55. Complex I is where hosts differ: the apes' complex I is reduced with 42 to 53 differing contacts, while
  *M. pahari*'s is normal with 100.
**Reading.** The cybrids confirm what divergence alone predicts: defects grow with the distance between host and
donor, and the donors that fail are the most distant. A count of differing contact residues adds little to it:
slightly better at telling reduced from normal complexes, worse at ordering the measured activities, and unable to
say which complex fails in a given pair. That agrees with H5, H9 and H14 (the coupling belongs to whole proteins),
and it means the interface mismatch scores used so far cannot be read as a measure of incompatibility. It neither
confirms nor refutes an interface mechanism: a count weighs every substitution alike, while the murid paper points
to a few cytochrome b residues near the Rieske protein's path. What could test that is an energy at the threaded
interface (FoldX or Rosetta on the chimeric complexes), which is heavier and belongs on haumea.
**Limits.** RefSeq individuals stand in for the donor cell lines except for cytochrome b. The mouse host is scored
on human structures through its own sequences, not on a mouse structure. Activities are percentages of each paper's
control, by different assays and normalisations. Failure to form a cybrid can have causes other than the respiratory
chain (replication, translation). T2b was added after the first run.
