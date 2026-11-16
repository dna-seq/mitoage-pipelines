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
discordance, ESM-2 embeddings, hydropathy with each gene's own composition) was rerun on the build after that,
0f4e970, which the manuscript cites; see H12. H2's power calculation was not rerun and describes the build before.

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
calibration (its log warns of a non-positive length 210 times, once per branch). "On the dated tree" above
means its topology only. Sampling about a hundred ages jointly with the
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
**Two-step result** (haumea, job 6490192, run from d6581f1; `results/h3/two-step/`). There are 209
branches, of which 106 are terminal and last 12 to 84 Myr. The synonymous fraction is 0.260. Terminal dS has a
median of 0.66 and a maximum of 1.97, so no branch passes 2 and the sensitivity rows equal the main ones.
Terminal dN/dS has a median of 0.053 (0.029 to 0.121). α is 0 on 7 short internal branches and on no terminal
one.
- **The positive control fails.** On mass alone, the slope of log10 dN/dS is +0.008 per log10 g (SE 0.009,
  p 0.34, λ 0.29). By the reading fixed above, nothing about lifespan is read from this test.
- **The secondary control holds.** dS per Myr falls with mass: −0.037 (p 0.0005, λ 0.99).
- **Not read, because the control failed.** Given mass, the lifespan slope of log10 dN/dS is +0.092
  (p 0.10). Given mass, the lifespan slope of log10 dS per Myr is −0.129 (p 0.016, Holm 0.033). The dS
  slope has Coevol's sign, but the rule fixed in advance does not let it count.
- **Rerun.** A second run in the same job is not byte-identical: HyPhy's 16-thread optimiser stopped
  0.1 log-likelihood units apart, with α at 0 on 1 branch instead of 7. Terminal dS and dN/dS agree within
  0.15 %. Every slope and p above agrees to three significant figures, and every verdict is the same.
- **What it says.** With one species per family, a terminal branch averages dN/dS over 12 to 84 Myr of
  ancestors, while the trait is today's. Popadin et al. (2007) found their body-mass effect on young
  terminal branches. This design does not reproduce it, so it cannot speak for dN/dS. That explanation is
  ours and has not been tested. The question stays with the Coevol chains.
**Coevol, resumed (2026-11-15).** Four chains were submitted on haumea in 3-day legs, with a checkpoint every
6 h, towards 15,000 points each. They use Coevol's MPI build on 8 ranks, at 64 s per point on the test node and 72 to 75 s on the
compute nodes, against 90 to 96 s serial; `docs/COEVOL_SPEED.md` compares the options.
- Chains 1 and 2: the published chains, continued. Jobs 6490194 to 6490197 and 6490198 to 6490201.
- Chains 1 and 2 of `coevol-fixed`: new chains holding the dated ages with -fixtimes. Jobs 6490202 to
  6490205 and 6490206 to 6490209.
Nothing from them is read until tracecomp gives effective sample sizes above 100 and small relative
differences.
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
**Second round, published.** The vertebrate ortholog sets, selections and alignments of complex IV's nuclear
subunits are in longevity-genie/mitoage-nuclear-partners at 7deb2e0 (`download prepared --topic partners`). The
contact map with its control pairing, the residue coupling of H5 (`coupling.tsv`, the first cluster run's), H9's
`lifespan.tsv`, both rounds' H12 tables and the ESM-2 means are in longevity-genie/mitoage-pipelines at 126d237,
under `output/` (build c525adb), which `download prepared --topic results` puts in `data/reference`.
**Second round, rerun on 0f4e970** (the build manuscript v0.3 cites: f75bdc9's EURING birds and the African grey
parrot's ND6; dna-seq/mitoage 03c0949, 2,344 species). Slurm jobs 6490380 (CPU) and 6490381 (GPU) on haumea,
2026-11-15, from commit 7ab48d2, with the vertebrate ortholog sets as published (nuclear-partners 7deb2e0), so that
only the build changes; tables and logs in `results/h12/round2/0f4e970/`. The checks are the same: 39 and 20 tests
pass, 19 of 19 and 2 of 2 code breaks are caught, and every table is byte-identical on a second run. The ESM-2
checkpoint is recorded by the job this time (08e4846, as above).
- Unchanged in every row: all mammal and fish rows of the three tables, and the whole partner table (the same 177
  mammals, 78 birds and 41 fishes have a partner sequence and enough pairs).
- Birds, 382 (379): hydropathy of COX2 with its own composition +0.030, q = 0.045 (0.041); under the primary model
  q = 0.085 (0.079). Contact mismatch alone -0.045, q = 0.054 (0.065). Embedding divergence, contact minus control,
  -0.007, q = 0.46 (0.50).
- Reptiles, 97 (the same species, realigned genes): contact minus control -0.111 and contact mismatch alone -0.114,
  both q = 0.046 (0.052), as in the first round's rerun on f75bdc9. It is read as before: the same species pass on one alignment and not on the
  other, which is not a finding.
- Nothing else crosses q = 0.05 in either direction. The mammalian ND5 and ND6 hydropathy, and every null of the
  second round, stand on the current build.
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
**Run.** Slurm job 6490112 on haumea, 2026-11-15, from commit 49d52d5, build c525adb; tables and log in
`results/job3-audit/`. The 28 tests pass, 23 of 23 deliberate breaks are caught, and a second run of every step gives
the same bytes. 2,335 species: 2,317 records are in release 237 (the 2,316 live ones and the retired *Andrias
davidianus*), 18 retired records are not, and every table carries those 18 as unknown. Everything below is a
candidate for a person to check, never a verdict.
**Result, strand: 2,201 conventional, 3 reversed, 7 mixed, 124 unknown** (106 outside Vertebrata, 18 without a
record). The sign of the GC skew follows the orientation almost without exception: negative in 2,200 of 2,201
conventional records and positive in all 3 reversed ones; the one conventional record with a positive GC skew is the
bonefish *Albula vulpes* (NC_082992.1, +0.043), whose genes are laid out conventionally: a reversal of the
strand asymmetry itself, as has been reported in some fishes (not checked here), not a deposit of the other strand. The AT skew follows less (positive in 2,034 conventional records, negative in 167, most of them fishes).
Against the build's own orientation column the two agree on every record except where the build asks for all 13
genes (5 records with one gene not annotated are conventional here and "noncanonical" there), where its chordate
rule differs (lancelet, ascidian; outside Vertebrata here), and the retired *Andrias davidianus* (unknown in the build,
conventional here from its record).
- Reversed (the whole record is the other strand; GC skew +0.27 to +0.28 in each): *Dasycercus cristicauda*
  NC_057634.1, *Pseudocaranx dentex* NC_058961.1, *Pseudoscaphirhynchus hermanni* NC_058992.1.
- Mixed: five *Trematomus* records with ND1 alone on the other strand (*T. eulepidotus* NC_057674.1, *T. hansoni*
  NC_086619.1, *T. loennbergii* NC_048965.1, *T. newnesi* NC_068275.1, *T. nicolai* NC_065826.1; one genus, so a
  shared annotation source or a real inversion, which the record cannot tell apart); ND6 alone on the + strand in
  *Anguis fragilis* NC_012431.1 and *Phocoenoides dalli* NC_053753.1.
**Result, gene order.** With tRNAs: 244 arrangements among 2,316 records (242 when a reversed deposit is not counted
apart); the commonest holds 61 %. Mammals 42 arrangements in 624 records, birds 26 in 418, ray-finned fishes 70 in
893, reptiles 28 in 129. Without tRNAs: 74 arrangements among 2,317 records; mammals 4, birds 8, fishes 10. *Crocidura
russula* NC_056768.1 is the one record with no tRNA annotated (and no rRNA). 71 flags in all, 53 at the full tier and 18
at the genes tier:
- Reversed deposits: the three reversed records above and the bivalve *Prolasmidonta heterodon* NC_037431.1.
- One gene or rRNA missing, extra or moved where the family shares one arrangement: *Pluvialis fulva* NC_033966.1
  (no CYTB), *Phalacrocorax auritus* NC_052822.1 (no ND6), *Cacatua pastinator* NC_040142.1 and *Turdus philomelos*
  NC_029147.1 (a second ND6), *Dryobates minor* NC_061952.1 (a second 12S rRNA), *Protonibea diacanthus* NC_024573.1
  (two 12S, no 16S), *Ostrea edulis* NC_016180.1 (no 12S), *Hemilepidotus gilberti* NC_034653.1 (16S on the other
  strand), *Trachypenaeus curvirostris* NC_050695.1 (both rRNAs on the other strand), *Galaxiella nigrostriata*
  NC_008448.1, *Rana amurensis* NC_030042.1 and *Sepioteuthis lessoniana* NC_007894.1 (genes in another order).
- With tRNAs: 30 records where one to five tRNAs sit on the other strand of an otherwise shared arrangement (tRNA-Gln,
  -Glu, -Pro, -Ala most often), 14 with tRNAs absent, extra or unresolved (e.g. *Carcharodon carcharias* NC_022415.1,
  *Galeorhinus galeus* NC_084376.1 with three extra tRNA-Thr, *Paleosuchus palpebrosus* NC_009729.1 missing three), 6
  with tRNAs placed differently.
All but the reversed deposits are "cannot tell from the record"; a single tRNA on the other strand in a family of
identical arrangements looks like an annotation slip, but the audit does not decide it.
**Result, genetic code: no mismatch.** Every one of the 2,317 records declares the code NCBI Taxonomy assigns to its
lineage (2: 2,211; 5: 97; 9: 6; 4: 2; 13: 1), and every one states a transl_table. Internal stops: none in 2,316
records under either code. One record is flagged, neither code clean: *Pelomedusa subrufa* NC_001947.1 (code 2; 6
internal stops in ND3 and 8 in ND4L). A likely reading, not checked: turtles' ND3 carries an untranslated extra
nucleotide that an annotation without join() reads through, and ND4L is read out of frame.
**Result, near-identical pairs: 84 pairs at identity >= 0.99, 15 of them >= 0.999; no two records with the same
sequence.** 75 pairs are within a genus (congeneric pairs this close are expected among recently split species: e.g.
*Agapornis lilianae* and *A. pullarius* NC_045369.1/NC_045368.1, identical over 11,346 bases; *Lepisosteus oculatus*
and *L. platyrhincus* NC_004744.1/NC_029715.1, 1 difference; *Saimiri boliviensis* and *S. sciureus*, 3; *Eubalaena
australis* and *E. glacialis*, 8; *Bos javanicus* and *B. taurus*, 16). 9 pairs are from different genera, the
strongest mislabel candidates:
- *Echis coloratus* NC_060592.1 and *Vipera berus* NC_036956.1: identical over 11,265 bases (*Echis carinatus*
  NC_060591.1 at 3 differences from both). The *Vipera berus* record sits inside *Echis*.
- *Anas platyrhynchos* NC_009684.1 and *Tadorna tadorna* NC_024750.1: 5 differences in 11,373 bases (0.99956).
- *Carcharhinus acronotus* NC_024055.1 and *Negaprion brevirostris* NC_072260.1: 0.99833.
- *Arius maculatus* NC_045222.1 and *Netuma thalassina* NC_030323.1: 0.99815.
- *Didelphis marsupialis* NC_057518.1 and *Metachirus nudicaudatus* NC_006516.1: 0.99727.
- *Catostomus rimiculus* NC_066944.1, *Chasmistes brevirostris* NC_066948.1 and *Deltistes luxatus* NC_054258.1: three
  pairs at 0.991 to 0.993, the Klamath suckers, among which hybridisation is reported.
**Result, out-of-place records: 136 out of order, 4 out of class, 102 unknown** (43 records with no order in the build,
39 alone in their order or class, 19 in no alignment). 86 of the 136 have their nearest neighbour in their own family:
the build's orders are not one scheme (AnAge's and NCBI's names side by side: Perciformes against Scombriformes,
Lutjaniformes or Centrarchiformes; Artiodactyla against Cetacea; Soricomorpha against Eulipotyphla), so these say
the labels disagree, not that a record is misplaced. The other 50 have their nearest neighbour in another family;
most margins are below 0.02 in amino-acid identity, at the depth where the 13 proteins resolve fish and bird orders
poorly. The candidates that stand out:
- *Seriolella porosa* NC_025665.1 (Centrolophidae, Perciformes) and *Piaractus brachypomus* NC_025315.2 (a pacu,
  Characiformes) are each other's nearest neighbour at 0.975, far above *Seriolella*'s best within its order (0.869):
  the *Seriolella* record looks like a characiform's.
- *Neogobius melanostomus* NC_052757.1 (Gobiidae) is nearest to the bonito *Sarda sarda* at 0.919, against 0.845 to
  any goby-like record; *Psenopsis anomala* NC_021460.1 is in turn nearest to it.
- *Cookeolus japonicus* NC_082750.1, *Monotaxis grandoculis* NC_010957.1, *Coryphaena equiselis* NC_087992.1,
  *Pterois volitans* NC_025290.1, *Alepocephalus bairdii* NC_013567.1 (margins 0.06 to 0.03), then birds whose nearest
  record is a flamingo (*Todus mexicanus* NC_052799.1, *Fregata magnificens* NC_052797.1, *Crotophaga sulcirostris*
  NC_052800.1, margins about 0.02), where deep bird relationships are the likelier reading.
- Out of class: the coelacanth and the lungfish (nearest to sturgeons, a scheme effect: Sarcopterygii is a class of
  its own here) and two bivalves nearest to gastropods at identities below 0.5.
**Reading.** The records are clean where it matters most for the analyses: every record declares the expected code
and only one translates with internal stops, 3 of 2,211 vertebrate records are deposited from the other strand (and
their GC skew says so), and the 13 protein genes keep one arrangement within almost every family. The candidates a
person should check first are the identical *Vipera berus* / *Echis coloratus* records, *Tadorna tadorna* at mallard
identity, *Seriolella porosa* at pacu identity, *Neogobius melanostomus* near the bonito, the four reversed deposits
(which the analyses reading strand-specific composition should orient or leave out), and *Pelomedusa subrufa*'s
frame. The out-of-order count mostly measures the build's mixed naming of orders, which is worth fixing in the build
before any analysis groups species by order.
**Rerun on 0f4e970.** The manuscript's round v0.3 cites build 2 at Hugging Face revision 0f4e970 (dna-seq/mitoage at
03c0949: the nine EURING birds and the patched ND6 of *Psittacus erithacus*, 2,344 species). Slurm job 6490297 on
haumea, 2026-11-15, from commit fb1cf00, in the same clone; tables and log in `results/job3-audit/0f4e970/` (those of
c525adb stay in `results/job3-audit/`). The 29 tests pass, 25 of 25 breaks are caught, and a second run of every step
gives the same bytes. Two things were added to the code for this build, each with a test and a mutation: the audit
compares every record's annotation as RefSeq release 237 deposits it with the build's current features
(`annotation_vs_build.tsv`, and a `build_annotation` column in strand, order and code), and it keeps reading the
record, not the patch. A first try (job 6490246, cancelled) also compared the feature rows the build keeps as retired
history and so named *Crocidura russula* and *Lemur catta*; only current rows are compared now.
- *What differs between the record and the build:* one record, *Psittacus erithacus* NC_027847.1, whose ND6 (CDS
  17376..17891, - strand) is in the build only, as the build's changelog records (`annotation_patch`). Strand, order
  and code read it without ND6 (conventional, "not annotated: ND6"; not flagged; code clean). The near-identity and
  nearest-neighbour steps read the build's alignments, which carry the patched ND6. The retired *Andrias davidianus*
  has no current features in the build; its record is audited as RefSeq has it.
- *The nine EURING birds* (*Tadorna ferruginea*, *Phylloscopus canariensis*, *P. ibericus*, *Certhia brachydactyla*,
  *Cyanopica cyanus*, *Aythya nyroca*, *Tichodroma muraria*, *Picus canus*, *Montifringilla nivalis*) are all in release
  237, all conventional, all translate cleanly under code 2, none is in a near-identical pair, and none is out of its
  order or class. One is flagged: *Montifringilla nivalis* NC_025911.1, tRNA-Gln on the other strand of the arrangement
  its family shares. *Tadorna ferruginea*'s nearest record is *Aythya marila* (0.971), not *Tadorna tadorna*, which
  fits the *T. tadorna* record being a mallard's.
- *Counts* (c525adb in brackets): 2,344 species, 2,326 records in the release (2,317), 18 without. Strand: 2,210
  conventional (2,201), 3 reversed, 7 mixed, 124 unknown, the same records as before; the build now calls *Psittacus*
  canonical, so four records, not five, are conventional here and "noncanonical" there. Gene order: 245 arrangements
  with tRNAs (244), 74 without (74), 72 flags (71: the one added is *Montifringilla*). Code: no mismatch in 2,326; the one
  record with stops under both codes is still *Pelomedusa subrufa*. Pairs: the same 84, byte for byte the same table.
  Out of place: 136 out of order and 4 out of class, as before; with the realigned proteins *Kryptolebias marmoratus*
  (margin 0.0004) leaves the out-of-order list and *Lipotes vexillifer* (margin 0.0003, a family-level label case:
  Cetacea against Artiodactyla) enters it, and identities move in the third or fourth decimal. Every candidate named
  above stands unchanged.
**What it cannot tell.** Whether a flagged arrangement is a rearrangement or an annotation error; whether a
near-identical pair is one species, a misidentification or a young pair of species; which record of a pair is the
mislabelled one. Each needs the specimen or the original submission.

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
bytes. The build is longevity-genie/mitoage at 0f4e970 (the site repository's 03c0949, the build round v0.3 of the
manuscript cites), and the dated tree is the site's phylogeny.json at 03c0949.
**Changed after the first job, before any readout existed.** Job 6490155 (from bf691b3, build c525adb) passed its tests
and both mutation checks (18 of 18 and 19 of 19 caught), then stopped in `mitoport iface` for the first complex other
than IV: the species of each protein's coupling were still taken from the COX1 alignment. They are now the species in
every mitochondrial alignment of the complex (`iface.pair_species`); for complex IV that is the same set, and its
tables stay byte-identical to the earlier code's. The readouts it had written were complex IV's only (its contact map,
coupling and partner tables, which H5 and H12 already report); `contact_readouts.tsv` had not been written. With the
rerun the build moves from c525adb (H12's) to 0f4e970, for consistency with the manuscript. Complex IV's species change
slightly, so its tables are no longer byte-identical to H12's; the job still runs complex IV through the earlier code
on the same inputs and compares the bytes, which is what shows that the generalisation changed nothing. The readouts,
their sets and the reading above are unchanged.
**Run.** Slurm job 6490344 on haumea, 2026-11-15, from commit 7ab48d2, build longevity-genie/mitoage 0f4e970 and the
site's tree at 03c0949; tables and log in `results/job1-contacts/`. 55 tests pass (one skipped: the pinned digests of
complex IV are of H12's build); 20 of 20 deliberate breaks of the new code and 19 of 19 of H12's are caught; every
table is the same bytes on a second run. On the same inputs, complex IV's contact map, `coupling.tsv`, `lifespan.tsv`,
`partners.tsv` and `partners_species.tsv` from this code are byte-identical to those of the code before the
generalisation, and its coupling is the same bytes whether or not the alignments of complexes I and V sit in the
folder. All numbers below are from `contact_readouts.tsv`; p values are two-sided, "Holm" is over the 12 interface tests.
**Result, R0 (is there whole-protein coupling?): yes, in every complex.** `seq_control_other`, subunits against the 29
cytosolic controls with a value, medians: complex I +0.49 against +0.09 (37 subunits,
p = 7.0 x 10^-7), III +0.40 against +0.02 (9, p = 6.8 x 10^-5), IV +0.42 against -0.02 (11, p = 6.4 x 10^-5), V +0.21
against +0.04 (15, p = 9.1 x 10^-4). H5's own contrast, `seq_contact_iface` of the subunits with an interface: I +0.47
(23 subunits, p = 8.7 x 10^-7), III +0.39 (5, p = 4.3 x 10^-4), IV +0.58 (11, p = 6.4 x 10^-5), V +0.37 (5, p = 0.046).
**Result, R1a (conservation, descriptive and not blind).** Identity to human, contact against matched control residues:
I 0.41 against 0.55 (621 pairs; paired median difference -0.029, p = 0.016), III 0.87 against 0.90 (104; -0.005,
p = 0.91), IV 0.92 against 0.97 (298; -0.021, p = 0.010), V 0.44 against 0.35 (112 pairs of 129 contacts; +0.037,
p = 0.08). Contacts are less conserved than their controls in I and IV, as H5 found for IV; not in III; in V the
difference leans the other way.
**Result, R1b (mitochondrial side, coupling): nothing.** `seq_excess_mt`, medians of subunits with an interface against
the controls: I +0.024 against +0.029 (p = 0.45), III +0.048 against +0.024 (p = 0.41), IV +0.049 against +0.014
(p = 0.36), V +0.010 against -0.006 (p = 0.71); Holm 1 in every complex.
**Result, R2 (nuclear side): nothing.** `seq_excess_nuc`: I -0.017 against -0.087 (p = 0.38), III -0.149 against -0.035
(p = 0.34), IV -0.021 against +0.030 (p = 0.72), V -0.133 against -0.055 (p = 0.54); Holm 1 in every complex.
**Result, R3 (pair concordance), mammals (primary).** Phylogenetic mean of contact minus control discordance per species:

| Complex | Species | Contact / control (raw means) | Mean difference | 95% interval | p | Holm |
|---|---|---|---|---|---|---|
| I | 177 | 0.251 / 0.284 | -0.047 | [-0.081, -0.012] | 0.009 | 0.11 |
| III | 177 | 0.138 / 0.193 | -0.027 | [-0.093, +0.040] | 0.43 | 1 |
| IV | 177 | 0.135 / 0.135 | +0.028 | [-0.012, +0.068] | 0.17 | 1 |
| V | 176 | 0.197 / 0.231 | -0.023 | [-0.080, +0.034] | 0.42 | 1 |

Other classes, reported and outside the correction: complex I in birds -0.012 ([-0.018, -0.005], p = 0.001, 75 species)
and fishes -0.008 ([-0.044, +0.028], p = 0.66, 51); III in fishes -0.056 ([-0.110, -0.002], p = 0.042) and birds -0.013
(p = 0.42); IV in birds -0.016 (p = 0.11) and fishes +0.016 (p = 0.25); V in birds -0.046 ([-0.090, -0.002], p = 0.042)
and fishes -0.056 (p = 0.19). Reptiles: not assessed in any complex (no species fills half the contact pairs with
vertebrate partners).
**Reading, as fixed.** R0 is positive in all four complexes, so each has whole-protein coupling to localise. No
interface readout passes Holm in either direction (the smallest adjusted p is complex I's R3, 0.11). Under the rule
written before the run, **H5's reading holds in complexes I, III, IV and V**: the nuclear subunits co-evolve with the
mitochondrial genes of their complex as whole proteins, and the residues that touch carry no detectable extra coupling
on the mitochondrial side, on the nuclear side, or as concordant pairs. Beyond the rule, and not counted: complex I's
pairs lean the predicted way in mammals (p = 0.009 uncorrected) and in birds (p = 0.001), the one place where the
contacts look more concordant than their controls; it is the lead to test again, not a finding. In complex IV the
phylogenetic mean leans the other way (+0.028) while the raw means are equal; neither excludes zero.
**Limits.** Complexes III and V have 5 subunits with an interface each, and the readouts compare them with 29 proteins:
a null there bounds little. One human structure stands for every species. The sequence distance only: no `esm_` column
is assessed. R1a was seen gene by gene before the entry was written. The pair test needs vertebrate partners, which
cover 176 to 177 mammals, 75 to 78 birds and 41 to 52 ray-finned fishes.
**Secondary, not read (lifespan).** `mitoport iface`'s lifespan model on contact minus control mismatch, mammals (604):
I +0.002 ([-0.021, +0.025]), III +0.010 ([-0.006, +0.027]), V -0.001 ([-0.017, +0.015]) log10 years per standard
deviation. `mitoport partners`, discordance contact minus control, mass and composition, mammals: I +0.025
([-0.003, +0.053]), III +0.005 ([-0.023, +0.033]), V +0.006 ([-0.026, +0.039]). No row of the partners tables reaches
q < 0.05 in any class; in the iface tables the only p < 0.05 is complex V's control mismatch alone in mammals (p = 0.048),
not the test variable.

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
**Run.** Slurm job 6490085 on haumea, 2026-11-15, from commit 2be936d (job 6490061, from af54126, stopped at
`frames`: a numpy bool in its rows broke the summary; fixed with a test, nothing else changed); build
longevity-genie/mitoage c525adb, dated tree from the site repository at 6fe8e72 (dated.nwk byte-identical to H12's).
Tables and log in `results/job4-peptides/`. The 60 tests of `tests/test_peptides.py` pass, each of 19 deliberate
breaks of the code makes a test fail, and a second run of every step gives the same bytes for all seven tables. The
steps take 7 minutes on 8 cores.
**Catalogue.** Eleven peptides placed, as fixed above; every one read once from the human row and followed by a stop.
**Question 1, the frame in each species** (`kept.tsv`; published reading; "start" is ATG under the standard code and
the NCBI start codons under the species' own). Species with the region, open at 90 %, with the start, kept:

Window null, open: H4's statistic, p over stop-free windows of the gene. Shuffle null: synonymous for the three
protein-gene frames (decisive), dinucleotide for the eight rRNA frames (context only); p's floor is 1/1,001;
n.a.: under 20 families.

| Peptide (reading) | Class | With the region | Assessed | Open >= 90% | Start codon | Open and start (kept) | Window null, open: p | Shuffle null, kept, one per family: p (q) |
|---|---|---|---|---|---|---|---|---|
| humanin (standard) | mammals | 624 | 623 | 557 | 463 | 412 | 0.098 | <= 0.001 (0.002) |
| humanin (standard) | birds | 418 | 418 | 381 | 357 | 332 | 0.165 | <= 0.001 (0.002) |
| humanin (standard) | reptiles | 129 | 128 | 65 | 49 | 34 | 0.233 | <= 0.001 (0.002) |
| humanin (standard) | amphibians | 45 | 45 | 42 | 28 | 26 | 0.068 | n.a. |
| humanin (standard) | ray-finned fishes | 891 | 891 | 497 | 865 | 478 | 0.402 | <= 0.001 (0.002) |
| humanin (standard) | cartilaginous fishes | 87 | 87 | 87 | 1 | 1 | 0.154 | 0.385 (0.615) |
| SHLP1 (standard) | mammals | 624 | 621 | 403 | 53 | 48 | 0.154 | <= 0.001 (0.002) |
| SHLP1 (standard) | birds | 418 | 417 | 167 | 29 | 5 | 0.531 | 0.130 (0.216) |
| SHLP1 (standard) | reptiles | 129 | 129 | 43 | 2 | 0 | 0.496 | 1.000 (1.000) |
| SHLP1 (standard) | amphibians | 45 | 45 | 22 | 3 | 1 | 0.335 | n.a. |
| SHLP1 (standard) | ray-finned fishes | 891 | 891 | 123 | 7 | 1 | 0.819 | 0.634 (0.976) |
| SHLP1 (standard) | cartilaginous fishes | 87 | 87 | 67 | 2 | 0 | 0.085 | 1.000 (1.000) |
| SHLP2 (standard) | mammals | 624 | 621 | 83 | 575 | 77 | 0.817 | <= 0.001 (0.002) |
| SHLP2 (standard) | birds | 418 | 413 | 91 | 13 | 1 | 0.746 | 1.000 (1.000) |
| SHLP2 (standard) | reptiles | 129 | 128 | 20 | 0 | 0 | 0.732 | 1.000 (1.000) |
| SHLP2 (standard) | amphibians | 45 | 45 | 8 | 0 | 0 | 0.656 | n.a. |
| SHLP2 (standard) | ray-finned fishes | 891 | 890 | 267 | 0 | 0 | 0.607 | 1.000 (1.000) |
| SHLP2 (standard) | cartilaginous fishes | 87 | 87 | 10 | 0 | 0 | 0.371 | 1.000 (1.000) |
| SHLP3 (standard) | mammals | 624 | 623 | 72 | 45 | 8 | 0.258 | 0.009 (0.019) |
| SHLP3 (standard) | birds | 418 | 413 | 106 | 7 | 3 | 0.226 | 0.037 (0.070) |
| SHLP3 (standard) | reptiles | 129 | 129 | 17 | 4 | 1 | 0.371 | 0.046 (0.084) |
| SHLP3 (standard) | amphibians | 45 | 45 | 9 | 5 | 1 | 0.177 | n.a. |
| SHLP3 (standard) | ray-finned fishes | 891 | 889 | 245 | 51 | 7 | 0.226 | 0.055 (0.096) |
| SHLP3 (standard) | cartilaginous fishes | 87 | 87 | 3 | 0 | 0 | 0.274 | 1.000 (1.000) |
| SHLP4 (standard) | mammals | 624 | 623 | 108 | 601 | 98 | 0.626 | <= 0.001 (0.002) |
| SHLP4 (standard) | birds | 418 | 418 | 39 | 3 | 0 | 0.868 | 1.000 (1.000) |
| SHLP4 (standard) | reptiles | 129 | 129 | 45 | 5 | 0 | 0.416 | 1.000 (1.000) |
| SHLP4 (standard) | amphibians | 45 | 45 | 33 | 0 | 0 | 0.164 | n.a. |
| SHLP4 (standard) | ray-finned fishes | 891 | 890 | 488 | 27 | 22 | 0.283 | 0.002 (0.004) |
| SHLP4 (standard) | cartilaginous fishes | 87 | 87 | 33 | 0 | 0 | 0.247 | 1.000 (1.000) |
| SHLP5 (standard) | mammals | 624 | 623 | 288 | 39 | 34 | 0.321 | 0.016 (0.032) |
| SHLP5 (standard) | birds | 418 | 418 | 203 | 0 | 0 | 0.450 | 1.000 (1.000) |
| SHLP5 (standard) | reptiles | 129 | 129 | 44 | 3 | 2 | 0.506 | 1.000 (1.000) |
| SHLP5 (standard) | amphibians | 45 | 45 | 15 | 1 | 1 | 0.480 | n.a. |
| SHLP5 (standard) | ray-finned fishes | 891 | 891 | 748 | 0 | 0 | 0.070 | 1.000 (1.000) |
| SHLP5 (standard) | cartilaginous fishes | 87 | 87 | 15 | 0 | 0 | 0.362 | 1.000 (1.000) |
| SHLP6 (standard) | mammals | 624 | 620 | 547 | 620 | 547 | 0.135 | <= 0.001 (0.002) |
| SHLP6 (standard) | birds | 418 | 418 | 405 | 417 | 405 | 0.108 | <= 0.001 (0.002) |
| SHLP6 (standard) | reptiles | 129 | 129 | 104 | 129 | 104 | 0.126 | <= 0.001 (0.002) |
| SHLP6 (standard) | amphibians | 45 | 45 | 43 | 45 | 43 | 0.056 | n.a. |
| SHLP6 (standard) | ray-finned fishes | 891 | 890 | 864 | 889 | 863 | 0.102 | <= 0.001 (0.002) |
| SHLP6 (standard) | cartilaginous fishes | 87 | 87 | 74 | 87 | 74 | 0.211 | <= 0.001 (0.002) |
| MOTS-c (standard) | mammals | 624 | 620 | 331 | 566 | 296 | 0.535 | <= 0.001 (0.002) |
| MOTS-c (standard) | birds | 418 | 413 | 342 | 334 | 272 | 0.234 | <= 0.001 (0.002) |
| MOTS-c (standard) | reptiles | 129 | 129 | 62 | 54 | 26 | 0.581 | <= 0.001 (0.002) |
| MOTS-c (standard) | amphibians | 45 | 45 | 14 | 18 | 7 | 0.738 | n.a. |
| MOTS-c (standard) | ray-finned fishes | 892 | 892 | 472 | 669 | 327 | 0.662 | <= 0.001 (0.002) |
| MOTS-c (standard) | cartilaginous fishes | 87 | 87 | 63 | 84 | 62 | 0.353 | <= 0.001 (0.002) |
| MTALTND4 (own, AGR=R) | mammals | 625 | 622 | 131 | 573 | 125 | 1.000 | 0.844 (1.000) |
| MTALTND4 (own, AGR=R) | birds | 418 | 416 | 158 | 379 | 133 | 1.000 | 0.103 (0.772) |
| MTALTND4 (own, AGR=R) | reptiles | 129 | 129 | 32 | 102 | 26 | 0.609 | 0.758 (1.000) |
| MTALTND4 (own, AGR=R) | amphibians | 45 | 44 | 8 | 31 | 7 | 0.957 | n.a. |
| MTALTND4 (own, AGR=R) | ray-finned fishes | 893 | 891 | 243 | 378 | 97 | 0.203 | 0.942 (1.000) |
| MTALTND4 (own, AGR=R) | cartilaginous fishes | 87 | 87 | 2 | 0 | 0 | 1.000 | 1.000 (1.000) |
| MTALTCO1 (own, AGR=R) | mammals | 625 | 625 | 18 | 386 | 13 | n.a. | 0.338 (1.000) |
| MTALTCO1 (own, AGR=R) | birds | 418 | 414 | 57 | 19 | 3 | n.a. | 0.998 (1.000) |
| MTALTCO1 (own, AGR=R) | reptiles | 129 | 128 | 2 | 32 | 0 | n.a. | 1.000 (1.000) |
| MTALTCO1 (own, AGR=R) | amphibians | 45 | 45 | 0 | 17 | 0 | n.a. | n.a. |
| MTALTCO1 (own, AGR=R) | ray-finned fishes | 893 | 890 | 4 | 117 | 1 | n.a. | 1.000 (1.000) |
| MTALTCO1 (own, AGR=R) | cartilaginous fishes | 87 | 87 | 0 | 10 | 0 | n.a. | 1.000 (1.000) |
| MENTSH (own, AGR stop) | mammals | 625 | 622 | 441 | 34 | 32 | 0.291 | 0.964 (1.000) |
| MENTSH (own, AGR stop) | birds | 418 | 418 | 342 | 39 | 37 | 0.299 | <= 0.001 (0.015) |
| MENTSH (own, AGR stop) | reptiles | 129 | 129 | 98 | 12 | 12 | 0.239 | 0.559 (1.000) |
| MENTSH (own, AGR stop) | amphibians | 45 | 45 | 32 | 6 | 6 | 0.239 | n.a. |
| MENTSH (own, AGR stop) | ray-finned fishes | 893 | 892 | 592 | 70 | 49 | 0.410 | 0.672 (1.000) |
| MENTSH (own, AGR stop) | cartilaginous fishes | 87 | 87 | 49 | 7 | 7 | 0.333 | 0.540 (1.000) |

Under the published reading, only SHLP6 is kept in most of every vertebrate class (81 to 97 %), as H4 found;
humanin is kept in 66 % of mammals, 79 % of birds and 54 % of ray-finned fishes, but in 1 of 87 cartilaginous fishes,
whose frame is open in all 87 and lacks the ATG; MOTS-c in 48 % of mammals and 66 % of birds. SHLP1 to SHLP5 are kept
in at most 16 % of any class (SHLP4 in mammals), and mostly in mammals. Of the frames inside protein genes, MTALTND4
is kept in 20 % of mammals and 32 % of birds; MTALTCO1 in 18 species (13 mammals: human, three macaques and the
Allen's swamp monkey, five carnivores, a bovid, the treeshrew and the aardvark, plus three birds, the coelacanth and a
parrotfish; no other great ape, as its authors found), at 41 to 54 % residue identity with human; MENTSH in 5 to 13 % of each
class. The standard reading of every frame is in `kept.tsv`.
**Question 2, beyond the host gene.**
- *Window null, statistic kept (the decisive one): not assessed for every peptide and class.* The genes hold almost
  no other ORF of the same length that starts with a start codon in human: under the published readings 0 to 9 such
  windows per peptide (21 only for MOTS-c read with the mitochondrial code, which is not its published reading),
  under the 20 the entry requires. Under the reading fixed in advance, **no frame can be
  called kept beyond its host's conservation**; the decisive null could not be run as designed. Requiring a start
  codon in human was too strict for genes this short; it is not relaxed here after the fact.
- *Window null, statistic open (H4's), every peptide and class assessed (61 to 345 windows):* nothing passes (smallest
  p under the published readings 0.027, SHLP4 in Malacostraca; every q >= 0.64, either reading). SHLP6 is open in 88 to 97 % of vertebrate classes against a
  window mean of 42 to 56 % (p 0.06 to 0.14): H4's result again.
- *Synonymous null (protein genes), one species per family:* MENTSH in birds is kept in 12 of 102 families against
  4.0 % expected from the shuffled codons (p <= 0.001, q = 0.015); in no other class (mammals 4.7 % against 8.8 %,
  p = 0.96). MTALTND4: no class (birds 35 % against 29 %, p = 0.10, q = 0.77). MTALTCO1: no class (mammals 1.9 %
  against 1.1 %, p = 0.34). At the species level (not read) MENTSH in birds (q = 0.024) and MTALTCO1 in mammals
  (q = 0.036) pass.
- *Dinucleotide null (rRNA, context only):* humanin, SHLP6 and MOTS-c are kept far above what their windows'
  composition gives in every vertebrate class with 20 families (amphibians have 16), humanin except in cartilaginous
  fishes (mammals, one species per family: SHLP6 86 % against 5.8 %, humanin 58 % against 2.1 %, MOTS-c 48 % against
  2.4 %); SHLP1, SHLP2, SHLP3 and SHLP5 in mammals only, SHLP4 in mammals and ray-finned fishes (q <= 0.032). As fixed in advance, this says the frames are more open than
  the rRNA's composition makes them, not that a peptide keeps them: the rRNA's own constraint does that too.
**Question 3, the tree** (`tree.tsv`, published reading, all species in the tree). SHLP6: kept at the root, 67 losses
and 15 gains among 2,276 tips (82 changes against 225 for random tips; clustered, p <= 0.001); in mammals 25 losses
and 6 gains. Humanin: 206 losses and 51 gains (257 changes against 745); MOTS-c 104 losses and 154 gains. The frames
inside protein genes are lost at the root and arise again and again: MTALTND4 19 losses and 262 gains, MENTSH 4 and
113, MTALTCO1 0 and 14; SHLP1 to SHLP5 likewise (0 to 5 losses, 15 to 49 gains). Every peptide's changes are
clustered on the tree overall (p <= 0.001). Read with the window null's failure: what the tree shows for the protein-gene
frames is many independent openings, not one old frame kept.
**Question 3, lifespan** (`pgls.tsv`, secondary). Of 25 fits with enough kept and lost species, one passes q < 0.05:
humanin in reptiles, +0.194 log10 years for species that keep the frame (23 kept, 73 lost; [+0.053, +0.336],
p = 0.008, q = 0.039, lambda 0.73). In mammals humanin goes the other way (-0.016, p = 0.19), in birds and fishes it
shows nothing (p 0.55, 0.61). Next smallest q: 0.22 (MTALTND4 in mammals +0.027, p = 0.029; MTALTCO1 in mammals,
13 kept, +0.062, p = 0.04). One class, the smallest, with no replication: not read as a lifespan signal.
**Reading.** The published peptides that look kept across species are the three in the most conserved stretches of
the rRNA genes (SHLP6, humanin, MOTS-c), and for them the only null that can separate the peptide from the rRNA
could not be run; the open-window version shows them no wider kept than other open windows of their gene. The frames
inside protein genes (MTALTND4, MTALTCO1, MENTSH) are not kept beyond what the host protein and its codon usage give,
in any class but one: MENTSH in birds, kept in 12 % of families against 4 % expected. A kept frame is not a
translated peptide.
**Not done.** RM16's fourth question (a blind scan of every frame of 20 codons or more); MASL, SHMOOSE and the other
left-out peptides (`left_out.tsv`); the window null with a criterion that has enough windows (it would be a new
test, fixed before it is run).
**Rerun on 0f4e970** (the build round v0.3 of the manuscript cites). Slurm job 6490249 on haumea, 2026-11-15, from
commit 00e5d12 (the analysis code of 2be936d; 00e5d12 adds only results), build longevity-genie/mitoage 0f4e970
(dna-seq/mitoage 03c0949: the 9 EURING birds and the parrot ND6 patch; 2,344 species, 2,325 not retired), dated tree
from dna-seq/mitoage at 03c0949. The tree matches the build: its 2,314 tips are all species of the build, the 9 new
birds are among them, and the 30 live species without a tip are the same 30 as on c525adb (they are read in `frames`
and `nulls`, and are absent from `tree` and `pgls`, as before). Again 60 tests pass, 19 of 19 code breaks are caught
and every table is byte-identical on a second run; tables and log in `results/job4-peptides/0f4e970/`, the c525adb
ones kept. **What changed.** The build realigned every gene with the new birds (all five host alignments have new
bytes); the catalogue is byte-identical. Of the 2,316 species in both builds, under the published reading, the state
is unchanged for humanin, SHLP6, MTALTND4, MTALTCO1 and MENTSH and changes in one species for MOTS-c, but not for the
frames on the 16S's other strand: SHLP2 25 species (718 change open), SHLP3 14, SHLP4 46 (399 change open), SHLP5 37, SHLP1 1. Their
columns move with the alignment, so SHLP1-5's counts are alignment-dependent and are quoted with that caveat. The
reading fixed in advance gives the same answers: the window null with ORF windows is not assessed for any peptide
under its published reading, so no frame is called kept beyond its host's conservation; H4's open-window null passes
nowhere (smallest p 0.0099, q >= 0.61); the synonymous null passes for MENTSH in birds alone (12 of 102 families
kept against 4.0 % expected, p = 0.002, q = 0.030). Context only, the dinucleotide null now also passes for SHLP2 in
birds and cartilaginous fishes and for SHLP4 in cartilaginous fishes (SHLP4 in ray-finned fishes moves to q = 0.050).
Tree: the same losses and gains for every frame but MTALTND4 (263 gains, from 262) and SHLP4 (4 losses and 39 gains,
from 5 and 43); every frame still clustered. Lifespan: unchanged but for the birds' n; humanin in reptiles again
+0.194 (q = 0.039), not read as a signal.


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

## H17. The cybrids again, with scores that weigh substitutions: co-evolved pairs and ESM-2 (H16's open item)

**Question.** H16's count of differing contact residues did about as well as divergence and could not say which
complex fails. A count weighs every substitution alike. Do scores that weigh them order the cybrids better: one that
asks whether the chimera breaks residue pairs vertebrates keep together, and a protein language model's view of the
donor's residues in the host's proteins? Molecular dynamics and interface energies are the heavier alternative; these
two are cheap and open, and are tried first.
**Fixed before any score was computed.** This entry was committed with the code, before the jobs ran. The measured
table, the points and the readouts T1 to T4 (and T2b) are H16's, unchanged (`mitoport cybrid-readout`).
- *Co-evolution* (`mitoport cybrid-coev`; the definitions are in the docstring of `mitoport/cybrid_coev.py`). The contact
  pairs and their control pairs are H12's and H14's, for the four complexes. Over the vertebrates with both proteins,
  each pair of columns gets the pointwise mutual information of a mitochondrial with a nuclear residue (sequences
  weighted at 80% identity, pseudocount 0.05, host and donor left out). A pair where donor and host differ loses
  PMI(host mt, host nuclear) − PMI(donor mt, host nuclear). The score, `coev_excess`, is the sum over contact pairs
  minus the sum over control pairs, and `coev_excess_mean` is that per pair scored. Sensitivities that enter no
  reading: pseudocount 0.5, and the host's whole order left out of the counts.
- *ESM-2 650M* (`mitoport cybrid-esm`; `mitoport/cybrid_esm.py`). Masked marginals in the host's protein: `esm_contact`
  (and its mean per contact residue) is how unlikely the model finds the donor's residues at the contacts, alone;
  `esm_partner_excess` is what placing the host's nuclear partner beside the protein adds at the contact residues,
  minus the same at their control residues.
- *Baselines*, on the same points: H16's `d_contact`, `d_all` and `f_all`.
- *Reading.* A score adds to divergence if, against `d_all`, it does all three: a higher area in T2 by at least 0.10
  (`d_all` has 0.74), a more negative Spearman in T1 (`d_all` has −0.36), and the most reduced complex named in at
  least 5 of the 6 pairs of T3, which requires *M. pahari* (complex IV) or *Otomys* (complex III) to be right (the raw
  counts get 4 by always naming complex I). Two of three is reported as partial and claimed as nothing. If no score
  does, the answer is that these sequence-based scores do not find the cybrids' defects, and an interface energy or
  molecular dynamics is what remains.
- *Known weaknesses, said in advance.* The points are few and not independent. The pairs' statistics come from 267
  to 415 vertebrates per subunit, far fewer effective sequences, so pair frequencies are coarse. The contact map is
  the human one, used for the mouse host through the alignment. ESM-2 was not trained on complexes; a partner joined
  by a glycine linker is a known trick with a weak signal, and pairs over 1,022 residues are not assessed. `esm_contact`
  alone sees no nuclear subunit: if it works, it measures how unusual the donor's protein is, not a mismatch.
**Test.** `slurm/h17-cybrid.sbatch` (CPU: the contact maps, the co-evolution score, the tests and the breaks of each
guard, a rerun check) and `slurm/h17-cybrid-esm.sbatch` (GPU: the ESM-2 scores, then the readout), on haumea, on the
inputs of H14's job (build 0f4e970, its vertebrate partner alignments).
**Result** (haumea jobs 6493605 and 6493606, run from 329d9a6; `results/h17/`). No score meets the rule; none meets
two of its three parts.

| score | T1 Spearman (13) | T2 area (9 against 20) | T3 most reduced complex named (of 6) |
| --- | --- | --- | --- |
| `d_all` (divergence, the baseline) | −0.36 | 0.74 | 4 |
| `coev_excess` | −0.06 | 0.73 | 4 |
| `coev_excess_mean` | −0.03 | 0.67 | 3 |
| `esm_contact` | +0.33 | 0.33 | 2 |
| `esm_contact_mean` | +0.41 | 0.33 | 1 |
| `esm_partner_excess` | −0.12 | 0.58 | 4 |

- *Co-evolution.* `coev_excess` names *M. pahari*'s complex IV, which no count does, but names complex IV for the rat
  and for *Otomys* too, where complexes I and III are the most reduced. Its sensitivities fall to chance in T2 (0.51
  with the pseudocount 0.5, 0.54 with the host's order left out), so what it has rests on the host's relatives. The
  pairs' statistics are thin where the proteins are conserved: at 80% identity the 418 species of complex IV count as
  3.3 effective sequences, complex III's 420 as 16.8 (complex I 74.5, complex V 58.9).
- *ESM-2 alone* goes the wrong way: complexes whose donor residues the model finds less likely are the less reduced
  ones (T2 0.33, T1 positive).
- *ESM-2 with the partner.* `esm_partner_excess` names *Otomys*' complex III, the 2% of control that no other score
  finds, and the three apes' complex I; it misses the rat and *M. pahari*, and across all points separates reduced from
  normal complexes poorly (0.58).
- *T4*, donors that gave no cybrid: both hamsters score above every donor that respired with `coev_excess`,
  `esm_contact` and `esm_partner_excess`; of the two orangutans only *P. pygmaeus* does, and not with `esm_contact`.
  Divergence places all four.
**Reading.** By the rule fixed in advance, these sequence-based scores do not find the cybrids' defects: neither the
co-evolved pairs nor ESM-2 orders the measured activities or separates reduced from normal complexes better than
divergence. Two single hits are not claims and are named only as what a later test could ask again: the partner
reading of ESM-2 at *Otomys*' complex III, and the co-evolution score at *M. pahari*'s complex IV. What remains for
the interface question is an energy on the threaded chimeric complexes, or molecular dynamics at a few sites.
**Limits.** Those said in advance, and: the 80% identity weighting, standard for families of proteins, leaves almost
no effective sequences for the conserved COX proteins, so complex IV's co-evolution score is close to unweighted
noise; the library did not report the checkpoint's commit; one partner pair of complex I exceeds the model's length in
each host and was not assessed.

## H18. The cybrids by an interface energy on the chimeric complexes (H16's and H17's open item)

**Question.** Counts (H16), co-evolved pairs and a language model (H17) do not find the cybrids' defects better than
divergence. Does a physical score: is the interface between the mtDNA-encoded and the nuclear subunits weaker, by a
force field, in the chimeric complex than in the host's own?
**Fixed before any energy was computed.** This entry was committed with the code, before any job ran; the protocol's
numbers are in the docstring of `mitoport/cybrid_energy.py`. One trial was made while writing the code, on one subunit
(human COX2 with gorilla's residues, 100 iterations, one start), to see that the steps run; its number enters nothing.
- *Models.* The host's complex as solved: human 5XTD, 5XTE, 5Z62, 8H9S; mouse 8OM1 (I), 7O3H (III), 7O37 (IV). So
  the mouse host is scored on mouse structures, not on the human map as in H16 and H17. The system is the complex's
  mtDNA-encoded chains and every chain that touches them, protein only. The chimera takes the donor's residue at each
  resolved position where donor and host differ; the reference is the host's complex through the same steps.
- *Energy.* AMBER14 with the GBn2 implicit solvent in OpenMM, on CPU. Side chains within 8 Å of a substitution are
  free, everything else is held; three starts each for chimera and reference. `ddE` is the chimera's interaction
  energy between the mtDNA-encoded and the other chains minus the reference's (kJ/mol; positive is a weaker
  interface), with `ddE_sd` from the starts. `ddE_per_change` divides by the substitutions at the interface.
- *Readouts and rule.* H16's table, points and T1 to T4, and H17's rule, unchanged: `ddE` (and separately
  `ddE_per_change`) adds to divergence if, against `d_all`, it has a higher area in T2 by at least 0.10, a more
  negative Spearman in T1, and names the most reduced complex in at least 5 of the 6 pairs of T3. Two of three is
  partial and claimed as nothing.
- *Known weaknesses, said in advance.* No cofactors: complex IV's hemes and copper, complex III's hemes and the
  iron-sulfur clusters are removed, and they sit at or near the interfaces in question. Implicit water for membrane
  proteins, with no membrane. An interaction energy without entropy, on a minimised structure, with a fixed backbone:
  a substitution that would move the backbone or change assembly is not seen. Insertions, deletions and unresolved
  residues are not modelled (counted per task). The energy's own noise is measured by the starts and may exceed the
  differences sought. Mouse complex V is not scored.
**Test.** `slurm/h18-energy.sbatch`, an array over the 44 (host, donor, complex) tasks on haumea's CPU nodes (no GPU),
then `slurm/h18-energy-check.sbatch`: the collection, the readout, the tests and the breaks of each guard.
**Result.** Not yet run.
