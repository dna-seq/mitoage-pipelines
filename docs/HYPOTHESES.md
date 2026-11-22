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
**Coevol, intermediate (2026-11-18, the first legs still running; provisional, not to be quoted).** Read from
copies of the traces as they stood, burn-in 1,000, with tracecomp and, per sampled covariance matrix, the
correlation and the partial correlation given the other three variables (order dS, dN/dS, lifespan, mass).
- *Resumed chains (ages sampled), 4,329 and 4,281 points.* The covariances with the traits have effective sample
  sizes of 950 to 2,037 and relative differences between the chains below 0.08. Tree length is not there yet
  (effective size 46), and the variance of dS (137) and the covariance of dS with dN/dS (146) barely are. dS
  against lifespan −0.49 and −0.49 (posterior probability of a positive value 0.00); dN/dS against lifespan +0.44
  and +0.43 (0.99, 0.98); dS against mass −0.38 and −0.37 (0.00); dN/dS against mass +0.22 and +0.20 (0.87,
  0.84). As partial correlations: dS with lifespan −0.12 and −0.12 (0.15, 0.14); dN/dS with lifespan +0.29 and
  +0.28 (0.91, 0.92); dS with mass −0.16 and −0.16 (0.07); dN/dS with mass −0.22 and −0.23 (0.12, 0.11). These
  are the numbers of the first 1,100 cycles again, now from chains that agree: the raw fall of dS with lifespan
  goes with body mass, and the rise of dN/dS with lifespan weakens with mass held and stays short of a posterior
  probability of 0.975.
- *Chains with the dated ages held, 3,449 and 3,299 points.* Not mixed, against what was expected of them: the
  covariances with the traits have effective sample sizes of 15 to 141 and the chains differ by 0.30 to 0.59 in
  relative terms. Nothing is read from them. Seen, and said so that it is not found later as a surprise: both
  chains at this point have a partial correlation of dS with lifespan of −0.57 and −0.53 (0.00), which is not
  what the resumed chains say. With the ages held, dS is a rate per unit of the dated tree's time; with the ages
  sampled and no calibration, rate and time are not told apart. If the difference holds once these chains mix,
  it is the fixed-ages reading that answers the question as stated, and the sentence above about body mass
  would not stand. The two-step test's lifespan slope of dS per Myr given mass (−0.129, p 0.016), not read
  because its control failed, points the same way; neither counts.
- At about 47 points an hour the four legs end short of 15,000 points; a fifth is submitted when the fourth ends.
- *After the first leg* (it ended on 2026-11-18, all four jobs completed; 4,528 and 4,480 points resumed, 3,643
  and 3,510 with the ages held; burn-in 1,000). Resumed: effective sample sizes of the covariances with the
  traits 603 to 1,963, relative differences below 0.08; tree length 48, variance of dS 147; the correlations and
  partial correlations as above to the second decimal. Ages held: effective sample sizes of the covariances with
  the traits 20 to 86, relative differences 0.37 to 0.66: still not read.
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
**Corrected, 2026-11-18** (`docs/ANALYSIS_CORRECTIONS.md`; tables in `results/review-20261118/longevity/`, run on
build 0f4e970 from the corrected source, commit 4c3f9fc). Two things changed in the calculation: the two sisters'
distances to their outgroup now use one set of codon columns, resolved in all three species, and a pair with a
zero rate difference leaves the sign test. All 655 rows were computed again. The result stands: nothing passes
q < 0.05 on the concatenated genes (closest: mammals with body mass in the model, slope +0.112, p = 0.022,
q = 0.16, 115 pairs; birds per million years, p = 0.054, q = 0.25, 102 pairs with saturated ones, and p = 0.073,
q = 0.27, 96 without), and the same four single-gene tests pass: ND6 across all species (-0.196, q = 0.0086),
avian CYTB per million years (-0.022, q = 0.0089), mammalian ND5 with body mass (+0.219, q = 0.011) and
mammalian ND2 on dN (-0.017, q = 0.036). No regression changes side at q < 0.05; two rows change side at nominal
p < 0.05. The sign tests moved more: 32 rows change side at nominal p < 0.05
(`results/review-20261118/sign_test_changes.tsv`), mostly on dN, where many pairs have no difference at all and
used to be counted as if they had a sign. The numbers in the paragraph above are those of the runs before the
correction.
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
**Limit, and it is the important one.** **The control could not be run**, and since 2026-11-16 the reason
is counted rather than assumed. `mitoage download technology` fetches the RefSeq record of each of the
build's 2,344 species and the INSDC record each one names as its source; `mitoage longevity technology`
reads the `Sequencing Technology` field of both into `technology.tsv` (`results/h8/technology.tsv`). All
2,344 RefSeq records name a source and every source was read. 893 sources (38%) state a technology; no
RefSeq record states one that its source does not. By class of technology: **25 long-read** (PacBio or
Nanopore, hybrids included), 642 short-read, 226 Sanger, 1,451 unstated. The 25 long-read records are 10
ray-finned fishes, 5 mammals, 5 birds, 2 cartilaginous fishes and one each of amphibians, gastropods and
insects: no class reaches the 25 species the model needs, so `motif-pgls`, which now reads the table,
still writes no `long_read` rows. An unstated record is never counted as short-read, and nothing is
inferred from a record's date. The mitoVGP caveat therefore stays open: the repeat, motif and
contact-zone numbers above are computed on records of which at most 1% are known to come from long reads.
A comparison this table would allow, and which was not fixed in advance and has not been run: whether the
features themselves differ between short-read and Sanger records of one class.

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
**After the array, 2026-11-21, before the result was written.** The array (job 6500884, 44 tasks) ended; the check job
(6500885) collected and made the readout, then stopped at one failed test, so the breaks did not run. The first
readout's table was seen at that point (`ddE` below `d_all` in T2); what follows changes no readout and no rule.
- *The failed test, and what it was.* The test parts the mtDNA-encoded chain from the others by 30 nm and asks that
  the interaction energy be zero; it gave 102 kJ/mol. Three small jobs (6572017, 6573071, 6574125) on the same
  structure: in double precision (OpenMM's Reference platform) the interaction of the parted groups is 0.000 in every
  term, so the function is right. In single precision (the CPU platform, which the tasks use) it was -2.7, 102 or
  -3,608 kJ/mol from one build to the next: the structure as built is not relaxed, each build places the hydrogens
  anew, and its energy ranged from 1.3e5 to 1.6e8 kJ/mol, where single precision is off by up to thousands. The
  test now asks its question in double precision. On a structure of ordinary energy the two platforms agree to
  under 1 kJ/mol in the interaction energy (23,370.5 against 23,370.7).
- *What that showed about the noise.* Each task builds its chimera and its reference once, and its three starts
  share the build. So `ddE_sd` holds the starts and not the build. Hydrogens are free in the minimisation, which is
  capped at 2,000 iterations on up to 100,000 atoms; whether two builds end at the same energy was not measured.
- *Control added: the same task from a second build.* The ten human tasks of complexes IV and V (the shortest) run a
  second time into another folder; `mitoage mitoport cybrid-energy-repeat` writes both `ddE`, their difference and
  the spread the starts alone would give. Reading, fixed now: if the median |difference| is above twice the median
  spread by the starts, `ddE_sd` understates the noise and a `ddE` is said to be distinguished from zero only beyond
  the difference seen between builds; the readouts T1 to T4 are entered as they are either way, since they use
  `ddE`'s value and not its sd. Ten tasks of two complexes say nothing certain about complexes I and III, which
  are larger.
**Result (2026-11-21; energies job 6500884 from 313f3e3, readout, tests and breaks job 6575325 from d2c3646; tables in
`results/h18/`).** **The interface energy adds nothing to divergence: one of the rule's three conditions holds for
each score, and one of three is claimed as nothing.** All 44 tasks were computed; 38 tests pass and 16 of 16 breaks
are caught.

| readout | `ddE` | `ddE_per_change` | `d_all` (divergence) | `d_contact` |
| --- | --- | --- | --- | --- |
| T1, Spearman with activity (13 points; 12 for per change) | -0.371 | -0.546 | -0.355 | -0.195 |
| T2, area for reduced against not reduced (9 against 20; 18 for per change) | 0.617 | 0.475 | 0.744 | 0.800 |
| T3, names the most reduced complex, of 6 pairs | 0 | 1 | 4 | 4 |

- *Against the rule.* T2 had to exceed divergence's by 0.10: both scores are below it. T1 had to be more negative:
  it is for both, by 0.02 for `ddE`. T3 had to reach 5 of 6: 0 and 1, where plain divergence names 4.
- *T4, the donors that do not respire.* The two orangutans have a higher total `ddE` than any human donor that
  respired (610 and 152 against at most 100 kJ/mol). The two hamsters do not (423 and 435 against 1,053 for rat,
  whose complex IV chimera has one start that did not settle: `ddE` 807 with a spread of 1,320).
- *The energies themselves.* `ddE` is positive, a weaker interface in the chimera, in 33 of 44 tasks, and beyond
  twice its spread over the starts in 31 (25 positive, 6 negative). It rises with the number of substitutions
  (Spearman 0.48 over the 44 tasks; 0.89 in complex I, 0.68 in III, 0.15 in IV, 5 tasks in V), so in the large
  complexes it mostly repeats divergence. 70 positions over all tasks were
  not modelled (deletions, insertions, unresolved residues).
- *The control: ten tasks from a second build (array job 6575326 from d2c3646; `results/h18/h18_repeat_build.tsv`).*
  **`ddE` does not repeat between two builds of the same task.** The median |difference| is 60 kJ/mol and the largest
  162, against a median spread by the starts of 21; 4 of the 10 change sign. By the reading fixed above (more than
  twice the starts' spread), `ddE_sd` understates the noise. Complex V repeats better (median 19, 4 of 5 keep their
  sign) than complex IV (median 96, 2 of 5), where all five moved up and the reference's own interaction energy
  came out 28 to 70 kJ/mol lower in the second build: the two runs differ by more than chance there, and what in
  (the nodes were different ones; the energy code was the same) was not found.
- *What that does to the bullet above.* The median |`ddE`| over the 44 tasks is 67 kJ/mol, the size of the
  difference between two builds. 25 tasks exceed the median difference and 9 the largest (7 of them in complex I,
  which was not built twice). So "beyond twice its spread in 31 tasks" is withdrawn: most single energies are not
  distinguished from the noise of the build. The readouts' verdict stands as the rule reads it, from `ddE`'s value,
  and it now has a second reason: the score is too noisy at this protocol to rank 13 points.
- *What a better run would need.* Several builds per task with `ddE` averaged over them, or one build shared by
  chimera and reference (substitute after the hydrogens are placed), and a minimisation that converges.
- *What this does not say.* That the interfaces of cybrids are not weaker: the weaknesses named in advance stand
  (no cofactors, no membrane, a fixed backbone, no entropy), and a score that misses for those reasons looks the
  same as one that misses because there is nothing to find. With H16 and H17 it says that no score tried here, from
  counts to a force field, ranks the published defects better than the number of differing residues.

## H19. The coupling in birds and fishes, its fall with distance from the mtDNA-encoded core, and what its mitochondrial side is

**Question.** H5 and H14 found that the nuclear subunits of complexes I, III, IV and V follow the mitochondrial genes
as whole proteins, in mammals, and that the residues that touch add nothing. None of it follows lifespan (H9, H12).
Three things about the coupling itself, of interest beyond lifespan:
1. Does it hold outside mammals? The whole-protein result in mammals is known (Weaver, Rabinowitz, Thueson & Havird
   2022, Mol Biol Evol 39:msac233, doi:10.1093/molbev/msac233: 58 mammals of eight orders; its introduction says
   such rate correlations had not been run in a vertebrate before). Birds and ray-finned fishes have not been read.
2. Does it fall with a subunit's distance from the mtDNA-encoded subunits of its complex? H14 says the touching
   residues carry nothing extra; this asks whether the touching subunits do.
3. What is the mitochondrial side? Does a subunit follow the mitochondrial genes of its own complex more than those
   of the others, and does it follow the mitochondrial proteins' change or the neutral rate of the genome?

**Seen before this entry was written, and therefore not blind.** On H14's tables (mammals, `coupling.tsv` of the four
complexes) an exploratory look, 2026-11-16, gave: subunits that touch a mitochondrial subunit against those that do
not, median `seq_control_other` 0.51 against 0.20 in complex I (27 and 10, Mann-Whitney p = 0.0006), 0.45 against
0.25 in III (6 and 3), 0.37 against 0.20 in V (6 and 9); coupling against the distance defined below, ranked and
centred within complex, Spearman -0.43 over 72 subunits (p = 0.0002); and complex IV's and III's nuclear subunits
coupled with the mitochondrial genes of every complex about alike (IV: 0.49, 0.48, 0.42, 0.48 with I, III, IV, V;
III: 0.46, 0.40, 0.34, 0.41), complex V's own subunits with ATP6 and ATP8 at 0.21. So in mammals readout G below is
a second look at a pattern already seen, and S is partly seen. Also seen, as the design: how many species of each
class the 50 vertebrate alignments of H14's job hold (medians 182 mammals, 78 birds, 120 ray-finned fishes, 19
reptiles per gene), and H12's and H14's discordance means in birds and fishes, which are another statistic. No
coupling had been computed in any class but mammals, no control protein had a vertebrate alignment, and the fourfold
transversion distance had not been computed for any species.

**Fixed before any result was seen.** This entry was committed with the code, before the job ran
(`mitoport/classes.py`; the definitions are in its docstring).
- *Classes.* Mammals, birds, ray-finned fishes and reptiles of the build (longevity-genie/mitoage at 0f4e970, the
  site's tree at 03c0949). A species counts when it is in the 13 mitochondrial protein alignments, in the 12
  heavy-strand codon alignments, a tip of the dated tree and in at least one nuclear alignment. A protein is read in a
  class with 40 species or more, a correlation with 200 usable pairs or more. Reptiles are expected below that (not
  assessed, and said).
- *Proteins.* NCBI's vertebrate ortholog sets (`download orthologs --ortholog-taxon vertebrata`), one protein per
  species (`mitoport select`), MAFFT L-INS-i on one thread, of: the 72 nuclear subunits of the four solved human
  complexes (I 37, III 9, IV 11, V 15), cytochrome c, complex II (4), the 7 matrix enzymes and the 30 cytosolic
  controls of H5, and 15 proteins of the inner membrane that are in no respiratory complex (`genes.INNER_MEMBRANE`:
  carriers, the pyruvate carrier, import, MICOS, LETM1, the calcium uniporter, prohibitins), chosen here as the
  membrane's control. The 50 sets H14's job downloaded are reused as they are; the others are downloaded for this
  run.
- *Coupling.* H14's R0 statistic within one class: over pairs of the class's species, the Spearman correlation of
  the nuclear protein's distance (share of differing residues over the columns the human protein fills, the columns
  that touch a mitochondrial subunit left out) with the distance of a set of mitochondrial genes at their
  burial-matched non-contact control residues, the patristic distance of the dated tree and its logarithm partialled
  out. Sets: the genes of complex I, III, IV or V, and all 13. No contact residue is on either side.
- *Distance from the core.* Per nuclear subunit, in the human structure (5XTD, 5XTE, 5Z62, 8H9S): the median, over
  its heavy atoms, of the distance to the nearest heavy atom of any mtDNA-encoded chain; a gene in several chains
  takes the chain with the smallest median.
- *Readouts, per class.* Groups of proteins are compared by Mann-Whitney, two-sided, with 5 proteins or more on each
  side.
  - C1 (question 1), per complex: its subunits against the cytosolic controls, coupling with the complex's own
    mitochondrial genes. Predicted: subunits above.
  - G1 (question 2): over the subunits of the four complexes, each one's coupling with its own complex's genes against
    its distance from the core, both ranked within complex, the correlation pooled; p by 10,000 permutations within
    complex (seed fixed). Predicted: negative. G2: the same with the subunit's own divergence (its mean distance over
    the pairs) partialled out within complex, because the far subunits are the conserved ones and a protein that
    hardly changes has little to correlate. G3: subunits with 5 interface columns or more against subunits with
    none, ranked within complex, permutation p.
  - M1, M2 (the membrane): the inner-membrane proteins against the cytosolic controls, and the subunits that touch
    against the inner-membrane proteins; coupling with all 13 genes. O1 to O3, reported: complex II, the matrix
    enzymes and the subunits that touch nothing, each against the cytosolic controls.
  - S1 (question 3), per complex: the coupling with the complex's own genes, the distance of the other complexes'
    mitochondrial genes partialled out as well (as its rank over the pairs); subunits against cytosolic controls.
    Predicted under complex-specific coupling: subunits above. S2: the reverse, the coupling with the other
    complexes' genes given the complex's own.
  - T1 to T3 (question 3), all subunits against the cytosolic controls, all 13 genes: T1, the coupling with the
    neutral distance (share of transversions at the fourfold-degenerate third positions two species share in one
    codon family, 12 heavy-strand genes, at least 100 sites); T2, the coupling with the mitochondrial proteins given
    the neutral distance; T3, the coupling with the neutral distance given the proteins. The controls take out what
    both genomes share through the lineage's generation time.
  - `class_matrix.tsv`: the median coupling of every group with every set, descriptive.
- *Correction.* Holm within class over C1 (four), G1, S1 (four), T2 and T3.
- *Reading.*
  - Positive control: in mammals C1 must give the subunits above the controls at p < 0.05 in all four complexes, as
    H14's R0 did through the mammalian ortholog sets (medians +0.49, +0.40, +0.42, +0.21). If it does not, the
    vertebrate sets or this code differ from H14's in something that matters, and nothing below is read.
  - Question 1: the coupling holds in a class and complex when C1 passes Holm there with the subunits above. Classes
    are compared by their medians only: they differ in species, in depth and in rate, and no test between classes is
    made.
  - Question 2: the gradient holds in a class when G1 is negative and passes Holm. It counts as replicated when it
    holds in birds or in fishes; the mammalian G1 is reported and is not a test (seen above). If G1 holds and G2 has
    p >= 0.05, the gradient cannot be told from the far subunits being conserved. If the inner-membrane proteins are
    above the cytosolic controls (M1, p < 0.05) and the touching subunits are not above them (M2, p >= 0.05), touching
    cannot be told from sitting in the inner membrane.
  - Question 3, specificity: a complex's coupling is complex-specific in a class when S1 passes Holm there with the
    subunits above. Where C1 passes and S1 does not, the subunits follow the mitochondrial genome and not their own
    partners more than the rest; S2 is reported beside it.
  - Question 3, rate: T2 passes Holm and T3 does not: the subunits follow the change of the mitochondrial proteins
    beyond the genome's neutral rate. T3 passes and T2 does not: they follow the neutral rate. Both: both. Neither:
    the two are not separated in these data.
  - A readout that cannot be computed is "not assessed" with its reason, is left out of the correction, and is never
    counted as a negative answer.
- *Known weaknesses, said in advance.* Species pairs are not independent, which is why only groups of proteins are
  tested; proteins of one complex are not independent either. One human structure stands for every species.
  Ray-finned fishes carry duplicated genes, and `select` keeps one protein per species (the one nearest the human
  length), so a fish's protein may be the other copy. The distances are uncorrected shares; fourfold transversions
  saturate late, not never, and the fishes are deep. The mitochondrial sets differ in size (control residues: I 621,
  IV 298, V 112, III 104), so in S1 a small complex's own side is the noisier one and the comparison leans against
  complexes III and V. The inner-membrane set was chosen by hand. The mammals go through the vertebrate ortholog sets
  and one species set per protein, so their numbers will not equal H14's. The ortholog sets were downloaded on two
  dates (H14's job and this run).
**Test.** `slurm/h19-classes.sbatch` on haumea: the tests and the breaks of each guard first
(`slurm/h19-classes-mutations.py`), the contact maps, the selections and alignments (the 50 of H14's job are compared
with its own, byte for byte), `mitoport classes`, and a second run to check that the tables come out the same bytes.
**Result** (haumea job 6493878, run from c63661b on 2026-11-16 and read on 2026-11-17; tables and log in
`results/h19/`, bd23aca). The 26 tests passed and 33 of 33 breaks were caught; the 50 alignments of H14's job came
out identical; the rerun gave the same bytes for all four tables. Species with the 13 mitochondrial genes, a tip and
a nuclear sequence: 182 mammals, 82 birds, 137 ray-finned fishes and 17 reptiles, so every reptile readout is not
assessed. Of the 30 cytosolic controls, 29 have a value in mammals, 19 in birds and 14 in fishes; of the 15
inner-membrane proteins, 15, 13 and 12. Complex II has 4 proteins with a value in every class, so O1 is not assessed.
Each cell: the medians of the two groups, then p (Holm p for the readouts in the correction).

| Readout | Mammals | Birds | Ray-finned fishes |
| --- | --- | --- | --- |
| C1, complex I | 0.49 / 0.17, 2.4e-6 (2.2e-5) | 0.30 / 0.19, 0.026 (0.21) | 0.31 / 0.12, 0.0013 (0.013) |
| C1, complex III | 0.40 / 0.01, 1.8e-4 (0.0013) | 0.56 / 0.31, 0.0037 (0.041) | 0.28 / 0.19, 0.025 (0.18) |
| C1, complex IV | 0.45 / 0.00, 9.3e-5 (7.5e-4) | 0.57 / 0.34, 0.033 (0.23) | 0.24 / 0.10, 0.0016 (0.015) |
| C1, complex V | 0.21 / 0.01, 0.0014 (0.0070) | 0.13 / 0.01, 0.037 (0.23) | 0.21 / 0.17, 0.40 (0.67) |
| G1, rank correlation | -0.39, 9.0e-4 (0.0054), seen before | +0.009, 0.94 (1) | -0.28, 0.031 (0.18) |
| G2, divergence partialled out | -0.40, 9.0e-4 | +0.14, 0.26 | -0.23, 0.076 |
| G3, 5+ interface columns / none | 0.49 / 0.20, 1.0e-4 | 0.32 / 0.26, 0.86 | 0.40 / 0.16, 5.0e-4 |
| M1, inner membrane / cytosolic | 0.15 / 0.05, 0.49 | 0.37 / 0.27, 0.11 | 0.04 / 0.19, 0.19 |
| M2, touching / inner membrane | 0.49 / 0.15, 3.5e-7 | 0.45 / 0.37, 0.68 | 0.32 / 0.04, 6.1e-5 |
| O2, matrix enzymes | 0.28 / 0.05, 0.14 | 0.36 / 0.27, 0.33 | 0.11 / 0.19, 0.44 |
| O3, touching nothing | 0.17 / 0.05, 0.064 | 0.37 / 0.27, 0.13 | 0.17 / 0.19, 0.93 |
| S1, complex I | 0.36 / 0.10, 1.4e-6 (1.4e-5) | -0.005 / -0.04, 0.33 (1) | 0.26 / 0.07, 1.3e-4 (0.0015) |
| S1, complex III | 0.007 / -0.04, 0.20 (0.61) | 0.42 / 0.20, 0.0037 (0.041) | 0.19 / 0.09, 0.035 (0.18) |
| S1, complex IV | 0.11 / -0.03, 0.034 (0.14) | 0.48 / 0.20, 0.085 (0.43) | 0.16 / 0.005, 0.0023 (0.018) |
| S1, complex V | 0.004 / -0.02, 0.29 (0.61) | -0.08 / -0.09, 0.37 (1) | 0.12 / 0.07, 0.26 (0.67) |
| S2, complexes I, III, IV, V | p 0.83, 2.4e-4, 0.0018, 0.0065 | p 0.23, 0.027, 0.60, 0.65 | p 0.021 (subunits below), 0.14, 0.86, 0.21 |
| T1, neutral distance | 0.30 / 0.16, 9.0e-4 | 0.45 / 0.28, 0.0070 | 0.16 / 0.21, 0.33 |
| T2, proteins given neutral | 0.35 / -0.01, 1.2e-7 (1.3e-6) | 0.18 / 0.11, 0.086 (0.43) | 0.23 / 0.16, 0.027 (0.18) |
| T3, neutral given proteins | 0.15 / 0.14, 0.78 (0.78) | 0.29 / 0.19, 0.010 (0.092) | 0.11 / 0.18, 0.22 (0.67) |

- *Positive control.* In mammals C1 puts the subunits above the cytosolic controls in all four complexes (Holm p
  2.2e-5 to 0.0070; medians 0.49, 0.40, 0.45, 0.21 against H14's 0.49, 0.40, 0.42, 0.21). The rest is read.
- *Question 1, classes.* The coupling holds in birds in complex III (Holm p = 0.041) and in ray-finned fishes in
  complexes I and IV (0.013, 0.015). Elsewhere in those classes it is not shown: birds' complexes I, IV and V have
  their subunits above the controls at nominal p 0.026 to 0.037 and do not pass Holm; fishes' complex III does not
  pass Holm (0.18), and complex V's subunits are not above the controls (p = 0.40). Reptiles: not assessed.
- *Question 2, distance.* Not replicated. In birds there is no gradient (+0.009, p = 0.94); in fishes it has the
  predicted sign (-0.28) and does not pass Holm (0.18). The mammalian -0.39, which holds with each subunit's
  divergence partialled out (G2 -0.40), was seen before this entry and is not a test. In mammals and in fishes the
  touching subunits couple more than inner-membrane proteins without a mitochondrial partner (M2 p = 3.5e-7 and
  6.1e-5), and those proteins couple no more than the cytosolic controls (M1 p = 0.49 and 0.19): there, touching can
  be told from sitting in the inner membrane. In birds it cannot (M2 p = 0.68).
- *Question 3, specificity.* In mammals only complex I's subunits follow their own complex's mitochondrial genes
  beyond those of the others (S1, Holm p = 1.4e-5). Complexes III, IV and V do not (0.61, 0.14, 0.61), and they follow
  the other complexes' genes given their own (S2 p = 2.4e-4, 0.0018, 0.0065): their subunits follow the mitochondrial
  genome, not their own partners more than the rest. In fishes complexes I and IV are specific (0.0015, 0.018), in
  birds complex III (0.041), the complexes where C1 passes in each.
- *Question 3, rate.* In mammals T2 passes (Holm p = 1.3e-6) and T3 does not (0.78): the subunits follow the change
  of the mitochondrial proteins beyond the genome's neutral rate, and the neutral rate adds nothing once the proteins
  are in. In birds and fishes neither passes, so the two are not separated there.

**Reading.** Outside mammals the whole-protein coupling is found in some complexes and not shown in the others: one
of four in birds, two of four in ray-finned fishes, with no complex common to both. Not shown is not absent: birds'
subunits sit above the controls in every complex at nominal p < 0.05. That a subunit
couples more the nearer it sits to the mtDNA-encoded core was not replicated in either class. In mammals the coupling
is with the change of the mitochondrial proteins, not with the neutral rate of the mitochondrial genome, so a lineage's
shared mutation rate or generation time does not account for it; and only complex I's subunits follow their own
complex more than the rest of the mitochondrial genome.
**Limits.** Those said in advance, and: birds and fishes have fewer controls with a value (19 and 14 of 30), which
widens their tests; birds' cytosolic controls couple more than mammals' in complexes I, III and IV (medians 0.19 to
0.34 against 0.00 to 0.17), which narrows the margin their subunits have to clear. G3 and the M and O readouts are outside the correction and reported only.

## H20. Temperature, size and metabolic rate against the mitochondrial substitution rate and composition

**Question.** Does the mitochondrial genome evolve faster in warmer animals and in animals with a faster metabolism,
as metabolic theory predicts? Gillooly, Allen, West & Brown 2005 (PNAS 102:140, doi:10.1073/pnas.0407735101)
combined metabolic theory with the neutral theory: substitution rates should scale with mass-specific metabolic
rate, so rise with temperature as exp(-E/kT) and fall with body size. The build has body temperatures for 253
mammals (AnAge), FishBase's preferred water temperature for 256 ray-finned and 40 cartilaginous fishes, and metabolic
rates for 206 mammals and 75 birds; none had been tested. Fishes are the strong case: a fish's body is at the
water's temperature, which ranges here from -2 to 28 C, against 31 to 40 C for mammals. Beside it, a published
relationship retested on the tree: Lehmann, Segal, Muradian & Fraifeld 2008 (Rejuvenation Res 11:409,
doi:10.1089/rej.2008.0676), from MitoAge's first authors, found that mtDNA GC content and resting metabolic rate
together explain 77% of the variance of mammalian maximum lifespan over 140 species, without the tree; Lehmann,
Muradian & Fraifeld 2013 (Front Genet 4:111, doi:10.3389/fgene.2013.00111) took body temperature in as well.

**Seen before this entry was written, and therefore not blind.** While the code was written, on the laptop's copy
of build 0f4e970 (species.parquet 267a16311cca) and its dated tree, 2026-11-16:
- the matching of the 980 fishes to FishBase (955 by exact name, 25 by synonym, none ambiguous, none in conflict
  with the build's FishBase record), and how many species have each trait (above);
- a first version of the rate test on the plain cherries of the whole tree, whose log10 dS per Myr fits were printed
  in full: mammals, size -0.056 (p = 0.067, 152 cherries), temperature and size: 1/kT -0.53 (p = 0.071, 29
  cherries; E = 1.2 eV); ray-finned fishes, size -0.12 (p = 0.13, 198), temperature and size: 1/kT +0.029
  (p = 0.28, 40; E = -0.07 eV), the midpoint of the temperature range and size: 1/kT +0.12 (p = 0.15) and size
  -0.41 (p = 0.0006, 35); birds, size -0.099 (p = 0.026, 110); cartilaginous fishes, size +0.29 (p = 0.41, 27);
  metabolic rate had fewer than 20 cherries in mammals (19) and birds (11), and temperature in cartilaginous
  fishes (7). dN per Myr and dN/dS were computed and not looked at;
- then, before any further fit was looked at, the pairing below, chosen because plain cherries use 40 of 256 fishes
  with a temperature: its counts of pairs per model only.
- The composition screen and the Lehmann model were run for their species counts only (mammals: 250 with a
  temperature, 203 with a metabolic rate and mass, 185 with both and a lifespan; Lehmann's model 201 and 183); no
  slope or p was printed.
So the pairing was changed after the first fits were seen. The plain cherries are therefore kept and reported
beside every test, and a result that holds under one pairing and not the other is said to be so.

**Fixed before the job ran.** This entry was committed with the code, before the job ran (`longevity/thermal.py`; the
definitions are in its docstring).
- *Data.* Build 0f4e970 and the site's tree at 03c0949, the default dataset's filters not applied (lifespan is not
  part of the rate or composition tests), retired records out. FishBase v26.06 from source.coop, the `species`,
  `synonyms` and `stocks` tables, each md5-checked (`download fishbase-temperature`).
- *Traits.* 1/kT, with T the body temperature in kelvin: AnAge's for mammals (270-320 K kept); for fishes, the median
  over a species' stocks of FishBase's `TempPref50` (-3 to 40 C kept), and as a sensitivity the median of the
  midpoints of `TempMin` and `TempMax`. Size: log10 adult mass, or log10 FishBase maximum length for fishes. Metabolic
  rate: log10 AnAge W, always with log10 mass in the model. A fish's FishBase record is found by exact name or by
  synonym (status accepted, provisionally accepted or synonym); a name with two codes, or one that disagrees with the
  build's own FishBase record, stays unmatched.
- *Pairs.* Per class and model, the species with every trait of the model, a tip and all 13 protein alignments, and
  the largest set of pairs of them whose paths on the dated tree share no branch (Maddison 2000, J Theor Biol 202:195, doi:10.1006/jtbi.1999.1050): walking up from the
  tips, members left unpaired below a node are paired there two by two in order of name. Per pair, dS and dN between
  the sisters over the 13 concatenated genes (Nei-Gojobori, Jukes-Cantor, as H7), divided by twice the age of their
  common ancestor; traits are the mean of the two sisters', missing when either lacks one. Pairs with dS above 1 or
  with dS or dN of 0 are left out and counted. Regression over pairs by GLS, each pair at its first species' tip,
  Pagel's lambda by maximum likelihood; at least 20 pairs. The same fits on the plain cherries, reported beside.
- *Tests* (`readout.tsv`), all on log10 dS per Myr under the primary pairing; E = -ln(10) x the slope on 1/kT, in eV:
  - controls: P1 mammals and P2 ray-finned fishes, the slope on size alone, predicted negative (rates fall with body
    mass across mammals: Nabholz, Glemin & Galtier 2008, MBE 25:120, doi:10.1093/molbev/msm248);
  - primary: T1 ray-finned fishes, T2 cartilaginous fishes, T3 mammals, the slope on 1/kT with size in the model,
    predicted negative (E > 0); M1 mammals and M2 birds, the slope on log10 metabolic rate with log10 mass, predicted
    positive;
  - retests of Lehmann et al.: L1, mammals of the default dataset, log10 lifespan on the GC content of the whole
    record and log10 metabolic rate, the GC slope on the tree, predicted positive; L2, the same with log10 mass and
    1/kT added. The OLS R^2 of the published kind is reported beside.
- *Correction.* Holm over the primary tests that could be assessed. The controls and L1, L2 are read at p < 0.05.
  Everything else (dN per Myr, dN/dS, the temperature range, the models with several traits, the composition screen)
  is reported with Benjamini-Hochberg q and read as exploratory.
- *Reading.*
  - A class's control must pass (negative, p < 0.05), or its primary tests are "not read": the rates of that class
    are then too noisy to show the known effect of size.
  - A primary test with Holm p < 0.05 in the predicted direction is "as predicted"; in the other direction,
    "opposite to the prediction"; else "not detected", which says no more than that this many pairs did not show
    it, with E's 95% interval saying how large an effect is left open. A test with fewer than 20 pairs is "not
    assessed".
  - Metabolic theory's value of E is not a test here: E is reported with its interval, and compared with the
    published value in the text only.
  - L1 "as published" when the GC slope is positive with p < 0.05 on the tree; the OLS R^2 is reported for comparison
    with the published 77%, and the difference between the two is what the tree takes out.
  - The composition screen counts features at q < 0.05 per class, model and term; no direction was predicted, and no
    single feature is read as a finding.
- *Known weaknesses, said in advance.* FishBase's `TempPref50` comes from AquaMaps, modelled from where a species is
  recorded, not measured in the animal; maximum length mixes total and standard length. AnAge's body temperatures and
  metabolic rates come from many studies of uneven quality and measure different states. Pairs with deep common
  ancestors saturate, which costs most in fishes (46 of 126 pairs in the temperature model); dN, which saturates
  later, is reported beside. Mammalian body temperatures span less than 10 K, so E is weakly bounded there. Nei-
  Gojobori distances are strained at the depth of many fish pairs, as H7 says.
**Test.** `slurm/h20-thermal.sbatch` on haumea: the tests and the breaks of each guard first
(`slurm/h20-thermal-mutations.py`), `longevity thermal`, and a second run to check that the tables come out the
same bytes.
**Run.** Slurm job 6494078 on haumea from 78f8734, build 0f4e970 (species.parquet 267a16311cca), the site's tree at
03c0949, FishBase v26.06 (all three tables at their pinned md5). The 21 tests pass, all 18 breaks of the guards are
caught, and a second run gives the same bytes for every table. Tables and the log are in `results/h20-thermal/`.
The pair counts differ from those seen while writing this entry, which used the laptop's older copy of the tree.
**Result.** log10 dS per Myr, pairs with no shared branch; the plain cherries beside:

| Test | Class | Pairs | Slope | p | Holm p | E, eV [95%] | Cherries: pairs, slope, p | Reading |
|---|---|--:|--:|--:|--:|---|---|---|
| P1 | mammals, size | 232 | -0.049 | 0.16 | | | 152, -0.056, 0.067 | control not passed |
| P2 | ray-finned fishes, size | 288 | -0.148 | 0.036 | | | 198, -0.122, 0.13 | control passed |
| T1 | ray-finned fishes, 1/kT with size | 74 | +0.096 | 0.10 | 0.30 | -0.22 [-0.49, +0.04] | 40, +0.029, 0.28 | not detected |
| T2 | cartilaginous fishes | 17 | | | | | 7 | not assessed |
| T3 | mammals, 1/kT with size | 71 | -0.435 | 0.048 | 0.19 | +1.00 [+0.01, +2.00] | 29, -0.53, 0.071 | not read (P1) |
| M1 | mammals, metabolic rate with mass | 47 | +0.023 | 0.91 | 1 | | 19, not assessed | not read (P1) |
| M2 | birds, metabolic rate with mass | 33 | +0.171 | 0.66 | 1 | | 10, not assessed | not detected |
| L1 | mammals, GC with metabolic rate (201 species) | | +0.018 per % GC | 0.0072 | | | | as published |
| L2 | mammals, GC with metabolic rate, mass, 1/kT (183) | | +0.017 per % GC | 0.0084 | | | | as published |

- **No primary test passes.** In mammals the size control did not pass (the mtDNA rate falls with mass, -0.049, but
  p = 0.16 on 232 pairs), so the mammalian temperature and metabolic-rate tests are not read; T3's nominal p = 0.048
  would not pass Holm either (0.19).
- **Fishes, where the test is strongest, lean the other way.** The size control passes, and the slope on 1/kT is
  positive, colder water with the faster rate, in both pairings. E is -0.22 eV, and its interval excludes every
  positive E above 0.04 eV: within these 74 pairs, a temperature effect of the sign metabolic theory predicts and of
  more than a small size is not compatible with the data. Gillooly et al. 2005 take E = 0.65 eV from the kinetics of
  metabolism, and their type I regressions over rates compiled across taxa without a tree averaged 0.40 eV (n = 4);
  the interval excludes both. (Read from PMC544068 on 2026-11-16, after this entry's result was first written.)
- Exploratory, not tests:
  - In ray-finned fishes dN per Myr follows 1/kT in the same direction, +0.095 per unit 1/kT (q = 0.0045, 74 pairs;
    E = -0.22 [-0.36, -0.08] eV), while dN/dS does not move (+0.003, p = 0.94). Colder-water fishes change their
    mitochondrial proteins faster per million years, in proportion to their synonymous rate. Size is in the model;
    depth, generation time and the pair ages are not, and FishBase's temperatures are modelled from occurrences.
  - dN per Myr falls with fish length (-0.18, q = 0.015, 288 pairs), and with the temperature range's midpoint in
    the model more steeply (-0.30, q = 0.0006).
  - dN/dS rises with body mass in mammals (+0.042, q = 0.042, 232 pairs), the nearly-neutral expectation, which H3
    used as its positive control.
  - Birds: dN/dS falls with mass-specific metabolic rate (-1.0, q = 0.042) on 33 pairs; too few to read.
  - Composition: no feature follows 1/kT at q < 0.05 in mammals (250 species) or ray-finned fishes (255); 6 of 43
    do in cartilaginous fishes (40 species, unreplicated). Metabolic rate follows none in mammals or birds. Lifespan
    with mass follows 11 features in mammals (614 species) and keeps 8 with body temperature and metabolic rate in
    the model (185): the composition signal of H11 is not carried by either.
- **Lehmann et al. reproduces in direction, not in size.** On the tree, GC content keeps a positive slope on lifespan
  with metabolic rate in the model (L1) and with mass and body temperature added (L2). The OLS R^2 of their model is
  0.49 on 201 species, not 0.77 on 140. Their metabolic rate is a stand-in for mass: with log10 mass in the model its
  slope is -0.02 (p = 0.70), while mass gives +0.15 (p = 0.0005). Body temperature, not a test here, carries
  lifespan beyond mass and metabolic rate on the tree: +0.34 per unit 1/kT (p = 1.7 x 10^-5, 183 species), cooler
  mammals living longer.

**Reading.** Metabolic theory's prediction for the mitochondrial clock is not supported in these data. In ray-finned
fishes, the class with the widest temperature range and a passing size control, the rate does not rise with water
temperature; dN per Myr rises in colder water, which is exploratory and opposite to the prediction. The mammalian
tests cannot be read, because the known size effect did not show on these pairs. Lehmann et al.'s association of
mtDNA GC with lifespan survives the tree and body temperature, at a smaller share of the variance, and their
metabolic rate adds nothing once mass is in the model.

## H21. Diadromous fishes: the mitochondrial genes, and their nuclear partners, on the branches where diadromy arose

**Question.** `docs/MITO_LOAD_ADAPTATIONS.md` ranks the adaptations that should load mitochondria by what the build can
test; the maintainer chose its shortlist on 2026-11-16 ("Fish ecology + mammals"): H21, H22, H23. Diadromy (moving
between fresh and salt water) couples long sustained swimming with an osmoregulatory switch, and it arose many times.
Across 401 teleost mitogenomes, migratory fishes had lower Ka/Ks in the mitochondrial protein genes than
non-migratory ones (Sun, Shen, Irwin & Zhang 2011, Mol Biol Evol 28:39, doi:10.1093/molbev/msq256), and sedentary
fishes showed relaxed purifying selection in all 13 genes (Strohm, Gwiazdowski & Hanner 2015, Gene 572:27,
doi:10.1016/j.gene.2015.06.074). Neither counted independent origins. Does the build show stronger constraint on the
mitochondrial genes where diadromy arose, counting each origin once, and do the nuclear subunits of the respiratory
complexes shift rate on the same branches?
**Seen before this entry was written, and therefore not blind.** The design tally of `docs/analyses/mito-load/`
(2026-11-16): how many species carry each candidate trait, how many are unknown, the trait clades and the fewest
changes on the dated tree, and how many are in H12's partner set. For diadromy: 122 yes, 368 no, 403 unknown among 893
ray-finned fishes; 69 trait clades, 59 changes; 22 yes among the 61 fishes with complex IV partners. No rate,
selection, composition or coupling statistic was computed for this or any other candidate trait. H11, H14 and H19's
design were read; H20's fits were not.
**Fixed before any code was written.** This entry is committed before the code. The code must implement these
definitions; anything changed before the run is listed under "Changed before the run", with its reason.
The design below is shared by H21, H22 and H23, which refer to it as "H21's design".
- *Data.* Build 0f4e970 and the site's tree at 03c0949, retired records out. Species of the class with a tip and all
  12 heavy-strand codon alignments (as H19); ND6 is reported apart, exploratory.
- *Trait.* FishBase v26.06 (`download fishbase-temperature`; matched by H20's rule, `thermal.fishbase_codes`),
  `species.AnaCat`: anadromous, catadromous, amphidromous or diadromous = 1; potamodromous, oceanodromous,
  non-migratory or oceano-estuarine = 0; empty, or no FishBase match = unknown. Unknown species are left out of every
  test and counted, never read as 0.
- *Trait clades.* On the dated tree pruned to the species with a known state: the maximal clades whose tips are all 1.
- *R1, RELAX (primary).* Tips: every trait species, and for each trait clade the three state-0 species nearest its
  root on the tree. Test branches: the stem and every branch inside each trait clade; reference: all others. HyPhy
  RELAX (2.5.104) on the 12 genes concatenated, each species' genetic code. *Control:* 20 draws, each replacing every
  trait clade by a state-0 clade of the same class matched on species count (within a factor of 2) and on stem age
  (nearest), without replacement (seed 0), tips and branches built the same way. Readout: K, the likelihood ratio and
  its p.
- *R2, one pair per origin (primary).* For each trait clade: the trait species nearest the clade's root (ties by name),
  its nearest state-0 relative outside the clade, and an outgroup, the known species nearest the root of the smallest
  clade holding both, outside it. Pairs are accepted in order of increasing sister distance, and one whose path shares
  a branch with a pair already accepted is dropped and counted. dN and dS of each sister to the outgroup over the 12
  genes (Nei-Gojobori with the Jukes-Cantor correction, each species' code, as H7). Readout:
  Δω = log10(dN/dS of the trait sister) − log10(dN/dS of the other), Wilcoxon signed-rank over pairs, two-sided; the
  sign test and Δ of dN and of dS alone reported beside it. Pairs with dS above 1 to the outgroup are left out and
  counted, and reported with them as a sensitivity. With size: the intercept of Δω regressed on the pair's difference
  in log10 size (FishBase maximum length for fishes, AnAge adult mass for mammals), reported beside.
- *R3, composition (exploratory).* H11's 43 features against the trait (0/1) with log10 size, phylogenetic regression
  with Pagel's lambda by maximum likelihood, Benjamini-Hochberg q over the 43; the count at q < 0.05. No direction is
  predicted and no single feature is read.
- *R4, aBSREL (exploratory).* HyPhy aBSREL on the stem branch of each trait clade, Holm over stems, and on the stems of
  the first control draw of R1; the two counts reported.
- *C, the coevolution side.* Species: the class's species in the tree with the 13 mitochondrial proteins and every
  protein of H19's vertebrate alignments of the nuclear subunits of complexes I, III, IV and V and of the 30 cytosolic
  controls (`genes.CYTOSOLIC_CONTROLS`). If fewer than 5 trait species remain, the proteins kept are those present in
  at least 90% of the species, and the species those that have all of them; the entry says which. Branch lengths per
  protein on the dated topology pruned to them (`mitoport branch-lengths`, H5's models). Per branch b: each protein's
  length scaled by its mean over branches; the nuclear background the mean over the cytosolic controls (a control's
  own length left out of its background); a nuclear protein's relative rate r = its scaled length / background; the
  mitochondrial rate m = the mean scaled length of the 13 / background. Foreground: the stems and inner branches of the
  trait clades present.
  - C1: per protein, s = median over foreground of log r − median over the other branches of log r. The subunits of the
    four complexes against the cytosolic controls, Mann-Whitney, two-sided.
  - C2: median over foreground of log m − median over the other branches, against 1,000 draws of matched state-0
    clade sets (R1's matching), two-sided empirical p.
  - At least 3 trait clades and 5 trait species in the species set, or C is not assessed.
  - C is read only if H23's positive control (bats) passes; otherwise C is "not read" in H21, H22 and H23.
- *Minimums.* R2 needs 8 accepted pairs, R1 3 trait clades; below that a readout is "not assessed".
- *Correction.* Holm over the primary readouts that could be assessed: R1, R2, C1, C2.
- *Reading.*
  - R1 "as predicted" when K > 1 (selection intensified on diadromous branches) with Holm p < 0.05 and the likelihood
    ratio above at least 19 of the 20 control draws'; K < 1 under the same two conditions is "opposite to the
    prediction"; p < 0.05 without the control condition is "not separated from the control clades", said as such;
    else "not detected".
  - R2 "as predicted" when Δω is negative with Holm p < 0.05; positive, "opposite"; else "not detected". The count of
    pairs and the median Δω with its interval are given, so the bound can be read.
  - C: no direction predicted. C1 and C2 passing Holm with the same sign is a "joint shift" of the mitochondrial genes
    and their nuclear partners on diadromous branches; one alone is said as such.
- *Secondary, reported and outside the correction.* Any migration (`AnaCat` anadromous, catadromous, amphidromous,
  diadromous, potamodromous or oceanodromous) against non-migratory (86 species), which is closer to Sun et al.'s
  contrast. The trait and its pair's differences in size and in H20's preferred temperature are tabulated.
- *Known weaknesses, said in advance.* `AnaCat` is empty for 45% of the fishes, and the empty ones are not a random
  draw. FishBase's labels are per species, not per population; a partly diadromous species is labelled by its record.
  Twelve genes are one linked molecule, so per-gene results are not independent. Fish pairs are deep and saturate,
  which costs dS first. Ray-finned fishes carry duplicated nuclear genes, and `mitoport select` keeps one copy (H19).
  The positive control is a mammal branch, so it shows that the method can see a known shift, not that it can in fishes.
**Code, and what it decided where the text above leaves room** (2026-11-16, written before any job ran). The code is
`src/mitoage_pipelines/adaptations/` (`traits`, `clades`, `pairs`, `selection`, `coupling`, `fits`, `run`), one command
per step under `mitoage adapt`; its docstrings hold every definition. Where this entry left a choice open, the code took:
- *Control clades.* "Matched on species count (within a factor of 2) and on stem age (nearest), without replacement"
  cannot give 20 different draws if "nearest" means one clade. Read as: each trait clade's candidates are the all-0
  clades within a factor of 2 in size; its pool is the 10 of them nearest in stem age; a draw takes one pool member at
  random for each trait clade (largest first) that shares no tip with a clade already in the draw; draw k uses a
  generator seeded with (0, k). A trait clade whose pool runs out in a draw is left out of that draw and counted. C2's
  1,000 draws use the same rule.
- *RELAX.* `--models Minimal` (the test alone); the codon columns of each gene filled in at least half of the tips, stop,
  partial and ambiguous codons written as gaps, as H3's inputs were made; every tip in the vertebrate mitochondrial code.
- *Which definitions get which readouts.* R1 (RELAX) runs on each hypothesis's primary definition only; the secondary
  and sensitivity definitions (H21's any migration, H22's depth) get R2 and R3.
- *The positive control's foreground.* In the species set of C, a branch is foreground when its tips of known state
  are not empty and all lie in one trait clade; a tip of unknown state neither joins nor breaks a clade.
- *Readout.* `adapt readout` writes every readout with its status, Holm within the hypothesis over the primary readouts
  that were assessed, and the reading fixed above (`run.reading`).
**Seen before the run, while the code was debugged** (2026-11-16, on the laptop's copy of build 0f4e970, species.parquet
267a16311cca; the dated tree of `longevity dated-tree`). One debugging run of `adapt pairs` printed its test table, and its
last 30 lines were read: for H22, deep_habitat's d_dN with saturated pairs and d_omega given size, and every row of
deep_depth; for H23, the marine rows (1 unsaturated pair and 4 with saturated ones, so not assessed) and hibernation's
rows (3 unsaturated pairs of 15, so not assessed; with saturated pairs, d_omega median -0.011, p 0.14; d_dS +0.021,
p 0.15); for H24's L3 the same pairs (d_dS with saturated pairs +0.015, p 0.15). H21's rows and H22's primary rows were
above those lines and were not seen. After it, local runs printed statuses only: `adapt fits` runs on every model, H25's
B2 has 11 bat pairs (fewer than its 12, so not assessed). The outputs were deleted unread. No definition above was
changed after these were seen; the counts (3 of 15 hibernation pairs unsaturated, 11 bat pairs) mean that H23's R2,
H24's L3 and H25's B2 will be "not assessed" by the rules as fixed, and they stay as fixed.
**Test.** On haumea, from a clone of main (each job's header gives the login-node steps):
- `slurm/adapt-main.sbatch`: the tests (`tests/test_adapt.py`, `tests/test_hibernation.py`), the breaks of each guard
  (`slurm/adapt-mutations.py`), `adapt traits`, `adapt pairs`, `adapt fits`, `adapt selection-prep` for the five traits
  of R1, all again into a second folder and compared byte for byte;
- `slurm/adapt-selection.sbatch`: RELAX on 5 traits x 21 draws and aBSREL on draws 0 and 1, one array task each;
- `slurm/adapt-coupling.sbatch`: H19's vertebrate alignments (made there if H19's job has not), `adapt coupling-prep`,
  `mitoport branch-lengths`, `adapt coupling` for mammals and ray-finned fishes;
- then `mitoage adapt readout`.
**Result.** Not yet run.

## H22. Deep-sea fishes: relaxed or adapted mitochondria where the deep sea was entered

**Question.** Two published readings point in opposite directions. Oxygen consumption falls with depth in
benthopelagic fishes and, less strongly, in benthic ones (Drazen & Seibel 2007, Limnol Oceanogr 52:2306,
doi:10.4319/lo.2007.52.5.2306), which predicts relaxed constraint on the mitochondrial genes. Yet 11 mitochondrial
genes carry convergent or parallel amino-acid changes across nine independent deep-sea lineages, read as repeated
adaptation (Shen, Pu, Chen, Murphy & Shen 2019, Front Genet 10:925, doi:10.3389/fgene.2019.00925). On the branches
where fishes entered the deep sea, counting each entry once, is the constraint on the mitochondrial genes relaxed,
intensified, or neither?
**Seen before this entry was written, and therefore not blind.** As H21: the design tally only. By habitat, 38 deep-sea
species of 893 ray-finned fishes, 32 trait clades, 32 changes, 2 in the partner set; by depth (≥ 1000 m), 52 yes,
519 no, 322 unknown, 41 clades, 40 changes, 8 in the partner set.
**Fixed before any code was written.** H21's design, with:
- *Trait (primary).* FishBase v26.06 `species.DemersPelag` bathydemersal or bathypelagic = 1, any other value = 0; no
  FishBase match = unknown.
- *Trait (sensitivity).* `species.DepthRangeDeep` ≥ 1000 m = 1, < 200 m = 0, 200 to 999 m left out, empty = unknown.
- *Readouts.* R1 and R2 primary, Holm over the two. R3 and R4 exploratory; R4 is where Shen et al.'s positive
  selection would show. C is expected to be not assessed with the partner set of 2026-11-16 (2 trait species); it is
  computed if H19's fish alignments hold at least 5.
- *Reading.* No direction is predicted. R1 with K < 1 under H21's two conditions is "relaxed: as the low-metabolism
  reading predicts"; K > 1, "intensified". R2: Δω positive with Holm p < 0.05 goes with relaxation, negative with
  intensification. R1 and R2 disagreeing is said as such. A sensitivity definition that disagrees with the primary is
  reported and not used to choose.
- *Known weaknesses, said in advance.* FishBase's depth habitat is a coarse label, and a species' recorded depth range
  is its deepest record, not where it mostly lives. Deep-sea fishes are cold, so a relaxation here cannot be told
  from H20's temperature without H20's result beside it. The deep-sea lineages of the build are the ones with
  mitogenomes in RefSeq, a fraction of those Shen et al. sampled.
**Test.** As H21.
**Result.** Not yet run.

## H23. Hibernating and marine mammals, with bats as the positive control: the mitochondrial genes and their nuclear partners

**Question.** The mammals are where the nuclear partners are (177 species in H12's set). Bats are one origin and so no
test of flight, but both the mitochondrial and the nuclear OXPHOS genes show adaptive evolution on the ancestral bat
branch (Shen et al. 2010, PNAS 107:8666, PMC2889356): a method that claims to see a joint mitochondrial and nuclear
shift must see it there first. Hibernators suppress metabolism to about 6% of basal rate during torpor bouts, daily
heterotherms to about 35% (Ruf & Geiser 2015, Biol Rev 90:891, doi:10.1111/brv.12137). Marine mammals dive, and
entered the sea four times among the build's mammals; across three of those origins, convergent substitutions were
as common in terrestrial sister taxa (Foote et al. 2015, Nat Genet 47:272, doi:10.1038/ng.3198). On the branches
where hibernation or marine life arose, do the mitochondrial genes change constraint, and do their nuclear partners
shift rate with them?
**Seen before this entry was written, and therefore not blind.** As H21: the design tally only. Bats: 44 species, 1
clade, 13 in the partner set. Marine: 50 species (cetacean families whether filed under Cetacea or Artiodactyla,
Otariidae, Phocidae, Odobenidae, Sirenia, *Enhydra*), 4 clades, 4 changes, 26 in the partner set. Hibernation was not
tallied: its source is not yet downloaded. A rough hand list of hibernating genera was run in the brainstorm and is not
used.
**Fixed before any code was written.** H21's design, with:
- *Positive control (P).* Order Chiroptera. C1 on the bat stem branch alone: the subunits of the four complexes against
  the cytosolic controls, predicted positive (the subunits faster on the bat stem). P passes with p < 0.05 in that
  direction, read alone, outside the correction. If P fails, C is "not read" in H21, H22 and H23. R4 (aBSREL) on the bat
  stem is reported beside it.
- *Hibernation (trait A).* The species table of Ruf & Geiser 2015, downloaded from the publisher's supplementary file
  and pinned by SHA-256 at the run step (a new source; the maintainer's call, said in `docs/MITO_LOAD_ADAPTATIONS.md`).
  Hibernator = 1. Reference = mammals not in that table as hibernator or daily heterotherm, called "not recorded as
  heterothermic", never "homeothermic": the review lists the species studied. Daily heterotherms are left out of both
  sides and reported apart. Names matched by exact name or by NCBI synonym; an ambiguous match is unknown. If the table
  cannot be obtained in machine-readable form, trait A is not assessed, and the entry says so.
- *Marine (trait B).* The rule above. With 4 trait clades, R2 (8 pairs) is not assessed; R1 and C are.
- *Readouts.* For each trait: R1, R2 and C (C1, C2) primary, R3 and R4 exploratory. Holm over the primary readouts of
  both traits that could be assessed.
- *Reading.* No direction is predicted for either trait. Each readout is read as in H21, two-sided: K < 1 "relaxed",
  K > 1 "intensified", under H21's two conditions; C1 and C2 with the same sign, a "joint shift".
- *Known weaknesses, said in advance.* The reference side of trait A holds heterotherms nobody has studied. Torpor is
  spread widely across mammals, so gains and losses are hard to tell apart, and the trait clades may count losses as
  separate origins. Four marine origins bound little whatever R1 says. One positive-control branch can pass or fail by
  chance; it is a check that the method can see a published shift, not a calibration of its error rate. The coupling
  of H14 belongs to whole proteins and probably follows the lineage's overall rate, so a joint shift on trait branches
  says the partners changed speed together there, not why.
**Changed before the run** (2026-11-16, before any analysis code of H21-H23 existed and before any statistic of these traits was
computed). Trait A's source. The species table of Ruf & Geiser 2015 could not be fetched by a script: the publisher's
page, PMC's copy and the authors' repository each answer with a bot challenge, which this project does not get
around. In its place, and better suited to a two-state trait because it codes non-hibernators as well:
- the supplementary table of Heldstab, Isler & van Schaik 2018 ("Hibernation constrains brain size evolution in
  mammals", J Evol Biol, doi:10.1111/jeb.13353; 1,104 mammals; Dryad doi:10.5061/dryad.753d06g, CC0), which Dryad
  serves only to a browser, so it is downloaded by hand once and accepted only with the md5 of Dryad's file listing
  (`download hibernation`, `adaptations/hibernation.py`). Its column and values are fixed when it is first read, before
  any statistic, and recorded here; until then trait A is not assessed;
- for bats, the trait table of Wilkinson & Adams 2019 ("Recurrent evolution of extreme longevity in bats", Biol Lett
  15:20180860, doi:10.1098/rsbl.2018.0860; DRUM hdl:1903/21501, SHA-256 pinned), `Hibernation` Yes = 1, No = 0; it
  matches all 44 bats of the build (42 by exact name, 2 by NCBI synonym): 17 hibernate, 27 do not, 4 changes on the tree.
Where both sources code a species and disagree, it is unknown, and the count is reported. A species in neither is
unknown, as before. Daily heterotherms are no longer a separate class (neither source has one).
The mammal table, saved by the maintainer from Dryad on 2026-11-16 and accepted by its md5
(c5e58d7f194bc830bc7ba854ad54aabe, as Dryad lists it), was then read for its layout only, before any statistic: one
sheet, "Appendix S1", 1,104 rows; the column `hibernation` holds 1 (80 species) or 0 (1,024), nothing else, and is read
with `Genus species` (`hibernation.read_mammals`). Its other columns (brain and body mass, diurnality, substrate, diet,
mid-latitude) are not read. Combined with the bat table (`hibernation.combine`) over the build's 625 mammals: 41
hibernate, 394 do not, 190 are unknown (in neither source: cetaceans, sirenians and others absent from a brain-size
dataset); 19 bats are coded by both sources, all alike, so no species is unknown by disagreement. 16 families of
hibernators, 25 trait clades and 20 changes on the dated tree; 20 hibernators are in H12's partner set of 177 mammals.
**Test.** As H21.
**Result.** Not yet run.

## H24. Hibernation and longevity in mammals: is part of lifespan's mitochondrial signal hibernation's?

**Question.** The maintainer, 2026-11-16: "check for hibernation idea as it is connected with longevity". It is, in the
literature. Hibernating mammals have about 15% higher annual survival than non-hibernators of similar size, mature
later and have longer generation times, yet do not live longer than non-hibernators with similar survival (Turbill,
Bieber & Ruf 2011, Proc R Soc B 278:3355, doi:10.1098/rspb.2011.0190): the link runs through lower mortality. Among
bats, hibernation is among the strongest predictors of longevity, hibernators living about 6 years longer on average
(Wilkinson & South 2002, Aging Cell 1:124, doi:10.1046/j.1474-9728.2002.00020.x), and body mass and hibernation predict
bat longevity again in a larger sample (Wilkinson & Adams 2019, Biol Lett 15:20180860, doi:10.1098/rsbl.2018.0860).
Ageing itself slows during hibernation: epigenetic age essentially stalls through the hibernation of yellow-bellied
marmots (Pinho et al. 2022, Nat Ecol Evol 6:418, doi:10.1038/s41559-022-01679-1), and big brown bats age more slowly by
a multi-species bat clock during hibernation (Sullivan et al. 2022, Proc R Soc B 289, doi:10.1098/rspb.2022.0635).
H11 found that mtDNA composition follows the slow life history of mammals as a whole, more than lifespan alone.
Hibernation is a slow life history reached many times over. Two questions: does part of what reads as lifespan's
composition signal belong to hibernation? And do hibernators, with months at a few percent of basal metabolism and
longer generations, accumulate mtDNA substitutions more slowly per year?
**Seen before this entry was written, and therefore not blind.** The design tally (`docs/analyses/mito-load/`),
including the bat row: 17 hibernating and 27 non-hibernating bats, 7 trait clades, 4 changes, 6 hibernators with
partners. The bat table's file was opened to find its columns; it also carries longevity, body mass and latitude, and
its rows were printed while the columns were found, but no relationship in it was computed. After this entry was
first committed, the mammal-wide table was placed and read for its layout and its counts (recorded under H23's change:
41 hibernators, 394 not, 190 unknown, 20 changes, 20 with partners); it carries no lifespan, and no statistic was
computed on it. No lifespan, rate or composition statistic of hibernators was computed in the build. H11's tables were
read when H11 was written.
**Fixed before any analysis code was written.**
- *Trait.* H23's trait A as changed before the run: Heldstab et al. 2018 for mammals, Wilkinson & Adams 2019 for bats,
  disagreement and absence unknown (both tables are in hand as of 2026-11-16; see H23's change).
- *Species.* Mammals of H1's default dataset (lifespan record of high or acceptable quality, a tip of the dated tree),
  with an adult mass, all 12 heavy-strand codon alignments and a known hibernation state. Build 0f4e970, the site's
  tree at 03c0949.
- *P, the positive control.* Phylogenetic regression (Pagel's lambda by maximum likelihood) of log10 lifespan on
  hibernation (0/1) and log10 mass. Predicted: the hibernation term positive, p < 0.05. P failing does not stop the
  readouts below, but it is said first, because the published link is then not reproduced on these records.
- *L1, composition follows hibernation (descriptive).* H11's 43 features, each standardised, on hibernation and log10
  mass, as H11 fits them; Benjamini-Hochberg q over the 43; the count at q < 0.05 and each slope.
- *L2, lifespan's composition signal with hibernation in the model (primary reading, no p).* For each of the 43
  features, the lifespan slope of feature ~ log10 lifespan + log10 mass, fitted with and without hibernation added, on
  the same species. Readouts: the count of features whose lifespan term has q < 0.05 in each model, and the median over
  those passing in the model without hibernation of |slope with| / |slope without|. Reading: a median ratio below 0.5
  with the count at least halved is "largely hibernation's"; a median ratio above 0.8 is "not hibernation's"; between,
  "partly". Fewer than 5 features passing without hibernation: "not assessed" (nothing to divide).
- *L3, the substitution rate (primary).* H21's design R2, one pair per origin of hibernation, with the readout on dS per
  million years: Δ = log10 dS of the hibernating sister to the outgroup − log10 dS of the other; Wilcoxon signed-rank,
  two-sided. Predicted negative (slower in hibernators). At least 8 pairs, or not assessed. Δ of dN and Δω reported
  beside it.
- *Correction.* L3 is the one p-valued primary test; it is read at p < 0.05. L1's q values are read as H11's are.
- *Reading.* L3 "as predicted" when Δ is negative with p < 0.05; positive, "opposite"; else "not detected", with the
  median Δ and its interval. L2 as defined above. Neither is read as hibernation causing anything: hibernators differ
  from their relatives in latitude, size and ecology too.
- *Known weaknesses, said in advance.* Heldstab et al.'s coding, its definition of hibernation and its coverage of the
  build are unknown until the file is read. Turbill et al. found the longevity of hibernators matched their survival,
  so a hibernation effect here may be one of lower mortality, which composition would follow only through generation
  time. In bats, Wilkinson & South found every hibernator's longevity record came from the field, so record type and
  hibernation go together. Latitude predicts longevity among hibernating bats (Wilkinson & Adams 2019) and is not in
  the build. Torpor's gains and losses are hard to tell apart, so L3's pairs may count a loss as an origin.
**Test.** As H21 (`adapt fits` writes P, L1 and L2; `adapt pairs` writes L3).
**Result.** Not yet run.

## H25. Inside bats: lifespan, hibernation and the mitochondrial genome where flight is held constant

**Question.** The maintainer: "bats look interesting". Bats outlive non-flying mammals of their size, and Austad & Fischer
(1991, J Gerontol 46:B47, doi:10.1093/geronj/46.2.b47) read such long lives as reduced vulnerability to the environment
rather than slow metabolism. As later reviews summarise it, the heart mitochondria of the little brown bat produced
about half to a third of the reactive oxygen species of a shrew's and a white-footed mouse's (Brunet-Rossinni 2004,
Mech Ageing Dev 125, doi:10.1016/j.mad.2003.09.003). Across mammals, our tests found no lifespan signal in the
substitution rate (H7) and a composition signal that follows the slow life history (H11), with body mass the confound
in every model. Within bats lifespan is far less tied to mass, and flight, the bats' one origin, is shared by every
species, so it drops out. Does lifespan follow the mitochondrial genome within bats, and does hibernation, which arose
about four times among the build's bats, change the constraint on it?
**Seen before this entry was written, and therefore not blind.** As H24: the bat row of the design tally and the
family-by-family counts of bats with a lifespan record of high or acceptable quality (43 of 44) and with partners
(13). No lifespan, rate or composition statistic of bats was computed.
**Fixed before any analysis code was written.**
- *Species.* The bats of the build with a tip of the dated tree, all 12 heavy-strand codon alignments and an adult
  mass; for the lifespan readouts, a lifespan record of high or acceptable quality as well. Build 0f4e970, the site's
  tree at 03c0949. Hibernation from Wilkinson & Adams 2019 as fixed under H23's change.
- *B1, the positive control.* log10 lifespan on hibernation and log10 mass, phylogenetic regression with Pagel's lambda.
  Predicted positive, p < 0.05. Most of the build's bat lifespans come from AnAge, as Wilkinson & Adams's do, so this
  checks the tree and the records, not the published result independently.
- *B2, the rate (primary).* H20's sister pairs (`thermal.disjoint_pairs`: the largest set of bat pairs whose paths share
  no branch), dS and dN per million years over the 12 genes; regression over pairs of Δ log10 dS per Myr on Δ log10
  lifespan and Δ log10 mass, by GLS with Pagel's lambda, as H20 fits its pairs. The lifespan slope, two-sided. At least
  12 pairs, or not assessed. dN per Myr and dN/dS reported beside it.
- *B3, composition (descriptive).* H11's 43 features on log10 lifespan and log10 mass, and on hibernation and log10
  mass, phylogenetic regression with Pagel's lambda; Benjamini-Hochberg q over the 43; counts at q < 0.05 and each
  slope, set beside H11's mammalian counts as numbers only (different species sets, no test between them).
- *B4, hibernation and constraint (primary).* H21's design R1 (RELAX) with the bat hibernation clades as test branches
  and the control draws among non-hibernating bat clades; it needs 3 trait clades and has 7. R2 (pairs) needs 8 and is
  not assessed.
- *B5, the coevolution side.* 13 bats have partners, 6 of them hibernators: within bats the coupling is not assessed.
  The bat stem as a whole is H23's positive control.
- *Correction.* Holm over B2 and B4.
- *Reading.* B2: a lifespan slope with Holm p < 0.05 is reported with its sign; else "not detected", with the slope's
  95% interval, which is what bounds the answer: about 20 pairs leave large effects open. B4 as H23's R1, two-sided,
  under H21's two conditions. B3 as H11 reads its counts.
- *Known weaknesses, said in advance.* About 43 species and four origins of hibernation: B4 bounds little, and B2's
  interval will be wide. Bat lifespans are maxima of banding records of uneven effort, and hibernators' come mostly from
  the field. Bats' small size puts them at one end of the mammalian mass range, so nothing here is compared with other
  orders as a test. The bat mitogenomes in RefSeq are not a random draw of bats (Vespertilionidae, Pteropodidae and
  Phyllostomidae are 38 of the 44).
**Changed before the run** (2026-11-16, while the code was written, before any job ran). B2 as written above cannot be
computed: in H20's design a pair's dS is measured between its two sisters, so a pair has one rate and no "Δ log10 dS"
between them. B2 is therefore H20's design as H20 fits it: log10 dS per million years of each pair on the pair's mean
log10 lifespan and mean log10 mass, GLS with Pagel's lambda, each pair at its first species' tip (`fits.bat_pairs`). Its
minimum of 12 pairs, its two-sided reading and its place in the correction are unchanged.
**Test.** As H21 (`adapt fits` writes B1, B2 and B3; B4 is RELAX on `bat_hibernation` in `slurm/adapt-selection.sbatch`).
**Result.** Not yet run.

## H26. Long-lived for their body mass: the mitochondrial genes where long life arose, counting each origin once

**Question.** The maintainer, 2026-11-13, asked "to check for convergent evolution of unrelated species selected for
longevity", and on 2026-11-16, asked whether that step should be its own framework or one more trait in H21's design:
"wait the code to land and then add longevity". Long life for a given size arose many times in mammals and in birds. The
longevity hypothesis says selection lowers the mitochondrial mutation rate in long-lived species (Nabholz, Glémin &
Galtier 2008, Mol Biol Evol 25:120, doi:10.1093/molbev/msm248; Galtier, Jobson, Nabholz, Glémin & Blier 2009, Biol Lett
5:413, doi:10.1098/rsbl.2008.0662), and mitochondrial proteins evolve more slowly in long-lived mammals and birds
(Galtier, Blier & Nabholz 2009, Mitochondrion 9:51, doi:10.1016/j.mito.2008.11.006). On the branches where species
became long-lived for their mass, is selection on the mitochondrial protein genes stronger than on matched branches
where they did not, each origin counted once?
**What this entry does not test.** Convergence at single sites (the same amino acid reached independently): H21's design
has no such test. CSUBST, PCOC and RERconverge, named in the first plan, are not run here.
**Seen before this entry was written, and therefore not blind.** H1 to H3, H7 and H11 and their results: no single
residue robustly associated with lifespan; H3's RELAX on long-lived against short-lived tips, where ND3 and cytochrome b
came out relaxed in one draw of species (K = 0.71 and 0.51) and not in a second; H3's Coevol chains, unconverged, with
dN/dS rising and dS falling with lifespan; H7, where lifespan does not predict the substitution rate on 339 sister
pairs; H11, where composition follows the slow life history. Those lean against the prediction below, and it is kept
as the published one. H21 to H25 are written and had not run. For the trait defined below nothing was computed: no
count of long-lived species, of trait clades or of pairs, and no residual.
**Fixed before any code was written.** This entry is committed before the code. "H21's design" is the text of H21 above
and its code under `src/mitoage_pipelines/adaptations/`, which this hypothesis calls and does not change.
- *Data.* Build 0f4e970 and the site's tree at 03c0949, retired records out. Two classes, each its own analysis and its
  own correction: Mammalia and Aves.
- *Species with a residual.* Species of the class in H1's default dataset (`aa_pgls.default_dataset`: lifespan record
  of high or acceptable quality, not estimated) with an adult mass above 0, a tip of the dated tree and all 12
  heavy-strand codon alignments.
- *Residual.* Ordinary least squares of log10 maximum lifespan on log10 adult mass over those species, one fit per
  class. The fit is not phylogenetic on purpose: the residual says how long a species lives for its mass, and the tree
  enters where origins are counted.
- *Trait `long_lived` (primary).* Residual at or above the class's 75th percentile = 1; below the class's median = 0;
  between the two = unknown ("intermediate, left out"). A species of the class without a residual is unknown. Unknown
  species are left out of every test and counted, never read as 0.
- *Trait `long_lived_in_order` (sensitivity).* The same cuts on the residual of a fit with one intercept per order and
  a common slope, over the orders with at least 5 species with a residual; species of smaller orders unknown. It asks
  the same of a species against its own order, so that bats or parrots as a whole do not carry the result.
- *R1, RELAX (primary), R2, one pair per origin (primary), R4, aBSREL (exploratory).* As H21's design, with its trait
  clades, control draws, pair rule, saturation rule, minimums and code choices. R1 and R4 on `long_lived`; R2 on both
  definitions. R2's size term is the pair's difference in log10 adult mass.
- *C, the coevolution side.* Mammals only, on the species set and branch lengths H23 prepares
  (`out/coupling/Mammalia`), with `long_lived` as the foreground; C1 and C2 as H21's design, and read only if H23's
  positive control passes. In birds C is not assessed: H21 to H25 prepare no bird set.
- *R3 is not run.* Composition against lifespan and mass is H1's and H11's, on the continuous trait.
- *Predicted.* R1: K > 1, selection intensified on the branches of long-lived clades. R2: Δω negative. ΔdS is reported
  beside R2 and predicted negative by the longevity hypothesis, outside the correction. C: no direction.
- *Correction.* Holm within each class over its primary readouts that could be assessed: R1, R2, and in mammals C1, C2.
- *Reading.* As H21's: R1 "as predicted" when K > 1 with Holm p < 0.05 and the likelihood ratio above at least 19 of
  the 20 control draws'; K < 1 under both conditions "opposite to the prediction"; p < 0.05 without the control
  condition "not separated from the control clades"; else "not detected". R2 "as predicted" when Δω is negative with
  Holm p < 0.05; positive, "opposite"; else "not detected", with the count of pairs and the median with its interval.
  A result in one class and not the other is said as such and is not called convergent across classes.
- *Known weaknesses, said in advance.* The cuts (upper quarter against lower half) are arbitrary and make a 0/1 trait of
  a continuous one; H3's Coevol asks the continuous question. Maximum lifespan grows with the number of animals
  recorded, so well-studied species look long-lived. A species at 1 beside relatives at 0 may be a record, not an
  origin. In mammals the long-lived clades will hold many bats, which are also H23's positive control for C; the
  sensitivity definition is the answer to that for R2, and there is none for R1. Mass is one adult value per species.
  Twelve genes are one linked molecule.
**Code, and what it decided where the text above leaves room** (2026-11-17, written before any job ran). The code is
`src/mitoage_pipelines/adaptations/longevity.py`, with its commands under `mitoage adapt longevity-*`; it calls H21's
`clades`, `pairs`, `selection`, `coupling` and `run.reading` and changes none of them.
- *The percentiles.* The 75th percentile and the median are taken over the species of the class that have that
  residual (pandas' linear interpolation); a residual equal to the upper cut is 1, one equal to the median is unknown.
- *The fit within orders.* Lifespan and mass are centred within each order of at least 5 species and one slope is
  fitted through all of them; its cuts are taken over those species of the class.
- *One class at a time.* The trait table holds both classes; pairs and selection inputs are made per class, from the
  class's species alone, so a control clade or an outgroup is never of the other class.
- *ΔdS in the readout.* A row "R2 dS", outside the correction, with the reading by its own p.
- *C's positive control.* Read from H21 to H25's `readout.tsv`; while that file is absent, C is "not read".
**Seen before the run, while the code was checked** (2026-11-17, build 0f4e970, species.parquet 267a16311cca, the dated
tree of `longevity dated-tree`). The counts of the trait: `long_lived` 151 yes, 302 no, 172 unknown in mammals and 96,
191, 140 in birds; `long_lived_in_order` 145, 289, 191 and 90, 178, 159. One check of the pair and readout steps was
meant to print statuses only and also printed the readings: mammals have 6 unsaturated pairs on the primary definition
(26 on the sensitivity one; 58 with saturated pairs there), so the mammals' primary R2 will be "not assessed" by the
rule as fixed; birds have 41 unsaturated pairs (49 on the sensitivity definition). The readings seen: "not detected" for
the birds' R2 on both definitions, for the birds' ΔdS and for the mammals' R2 on the sensitivity definition. No
estimate, median or p was printed, and the tables were deleted unread. No definition was changed after this; the
mammals' minimum stays as fixed. R1 had not run.
**Test.** On the cluster, from a clone of main: `slurm/adapt-longevity.sbatch` (the tests `tests/test_adapt_longevity.py`,
the 15 breaks of `slurm/adapt-longevity-mutations.py`, `adapt longevity-traits`, `longevity-pairs`,
`longevity-selection-prep` for both classes, all again into a second folder and compared byte for byte); then
`slurm/adapt-longevity-selection.sbatch` (RELAX on 2 classes x 21 draws, aBSREL on draws 0 and 1); `adapt
longevity-coupling` on H23's prepared mammal set; `adapt longevity-readout`.
**Result, R2 only** (2026-11-17, Slurm job 6497484 from f9a914a; 9 tests passed, 15 of 15 breaks caught, every table
and selection input the same bytes on the rerun; tables in `results/h26/`). R1, R4 and C have not run.
- *Birds, `long_lived`:* 72 origins, 41 unsaturated pairs. Median Δω +0.009 (95% interval −0.014 to +0.018), Wilcoxon
  p = 0.67; 24 of 41 positive. **Not detected**; the interval excludes a shift of dN/dS larger than about 4% either
  way. ΔdS median −0.0005 (−0.011 to +0.003), p = 0.49. Given the pair's difference in mass, the intercept is +0.003
  (p = 0.71). On `long_lived_in_order`: 49 pairs, median +0.006, p = 0.61.
- *Mammals, `long_lived`:* 39 origins, 6 unsaturated pairs: **not assessed**, as the counts said in advance. With the
  19 saturated pairs added (a sensitivity, outside the reading): median Δω −0.025 (−0.044 to +0.016), p = 0.22. On
  `long_lived_in_order`: 26 unsaturated pairs, median +0.004 (−0.022 to +0.017), p = 0.94.
- Holm waits for R1: R2 is one of each class's primary readouts, and no p here is below 0.05 before any correction.
So far: where birds became long-lived for their mass, the ratio dN/dS of the 12 protein genes did not move against the
nearest relative that did not, and neither did dS; in mammals the pairs are too saturated to say on the main
definition, and show nothing on the one within orders.
**R2 corrected, 2026-11-18** (`docs/ANALYSIS_CORRECTIONS.md`; tables in `results/review-20261118/adapt-longevity/`,
from the corrected source, commit 4c3f9fc). The two sisters' distances to their outgroup now use one set of codon
columns, resolved in all three species; the pairs and their counts are the same. The readings stand.
- *Birds, `long_lived`:* 41 unsaturated pairs, median Δω +0.009 (−0.014 to +0.018), Wilcoxon p = 0.66; 24 of 41
  positive: **not detected**. ΔdS median −0.0005 (−0.011 to +0.003), p = 0.50. Given the difference in mass, the
  intercept is +0.003 (p = 0.68). On `long_lived_in_order`: 49 pairs, median +0.004 (−0.018 to +0.021), p = 0.58.
- *Mammals, `long_lived`:* 6 unsaturated pairs: **not assessed**. With the 19 saturated pairs added: median Δω
  −0.023 (−0.046 to +0.016), p = 0.21. On `long_lived_in_order`: 26 unsaturated pairs, median +0.004 (−0.024 to
  +0.018), p = 0.88.
The numbers of 2026-11-17 above and the tables in `results/h26/` are those of the run before the correction. R1
and R4 (RELAX, aBSREL) take their species from the pair table and do not use these distances, so the running array
is not affected; the readout and C run from a clone at the corrected commit or later.

**Result, R1 and R4, and the readout** (2026-11-20; RELAX and aBSREL: Slurm array 6500930 from f9a914a, 46 of 46
tasks; readout: job 6553573 from d20393f, which made the trait and pairs again, found them the bytes of the corrected
tables, prepared the selection inputs again and found every one the bytes the runs read; tables in
`results/h26-readout/`). A first readout job (6553170) stopped in H21's `selection.summarise`, which took `draws.tsv`
for a draw folder; fixed in d20393f with a test that fails on the code before it. No rule changed.
- *Mammals, R1:* K = 0.74, likelihood ratio 563 (p below machine precision), above all 20 control draws (their
  largest 74.7; their K from 0.89 to 1.10). **Relaxed: opposite to the prediction.** 39 trait clades, 202 tips, 262
  test branches.
- *Birds, R1:* K = 0.94, likelihood ratio 18.7 (p = 1.5e-5, Holm 3.1e-5), above 19 of 20 control draws (one control
  draw reaches 19.9 with K = 0.94). **Relaxed: opposite to the prediction**, by the rule as fixed, at its edge: one
  more control above it and the reading is "not separated from the control clades".
- *R2 after Holm:* birds p = 0.66, **not detected**; mammals **not assessed** (6 pairs).
- *R4, aBSREL (exploratory):* no stem of a long-lived clade with episodic positive selection in either class (0 of
  39 in mammals, 0 of 72 in birds, before and after correction); in the first control draw 1 of 39 and 2 of 72.
- *C:* **not assessed**: H23's mammal set does not exist yet.
**Reading.** Where mammals are long-lived for their mass, selection on the 12 protein genes is weaker, not stronger,
than on the rest of the tree, and far more so than on matched control clades; birds show the same direction, a small
effect barely apart from the controls. No sign of positive selection at the origins. This is the direction H3's Coevol
chains lean to (dN/dS rising with lifespan, mass given) and the opposite of the published prediction kept here.
Two things limit it, both said in advance: of the 151 long-lived mammals 87 are primates and 25 bats, so R1 in mammals
is largely primates against other mammals, a few deep origins and not 39 independent ones, and the definition within
orders, which answers that, exists for R2 only (where it shows nothing); and RELAX's K < 1 does not separate relaxed
purifying selection (as a smaller population gives) from weaker selection of any other cause. R1 on
`long_lived_in_order` was not fixed in this entry and would be a new hypothesis.

## H27. Do the results hold on records of high confidence, and on well-sampled lifespans?

**Question.** The maintainer, 2026-11-18, after the first author of the 2015 database reported mistakes in its
curation: "we should give all, but give confidence level (like classiciation) and for fast to compute things we have
to provide extra recomputed with high quality data only"; and, the same day, on maximum lifespan: "for non-popular
wild species we probably have dicreases MLS as we know only very small number of animals while for domestic and
popular zoo animals we know wayy more long-livers". The site's build now gives every species a confidence level with
its reasons (`data_confidence` in the species table of `dna-seq/mitoage`: high, caveat, low, unreliable). Do the
results of the earlier hypotheses change when the species set is restricted by that level, or by how well the
lifespan is sampled?
**Seen before this entry was written, and therefore not blind.** The counts of the levels (high 1,896, caveat 275,
low 111, unreliable 62 of 2,344) and the two audit reports behind them (`docs/analyses/curation-audit-2026-11-18/`
in `dna-seq/mitoage`). One look at the bias, unadjusted for phylogeny, on the default dataset: the residual of log10
lifespan on log10 mass within a class has a median 1.26 times the expectation in birds with a large AnAge sample
and 0.93 with a small one, 1.24 in captive and 0.95 in wild bird records; in mammals 0.86 (small), 1.03 (medium),
0.95 (large). No result of H1 to H26 was recomputed on any restricted set before this entry.
**Fixed before any code was written.**
- *Species sets.* For each analysis below, the species it used, restricted four ways:
  `not_unreliable`: the confidence level is not unreliable (the site's default dataset from 2026-11-19);
  `high`: the level is high;
  `well_sampled`: `not_unreliable`, and the lifespan is from AnAge with a sample size of medium, large or huge;
  `captive` and `wild`: `not_unreliable`, and the lifespan record is of that kind (mammals and birds only, each
  only where at least 50 species remain).
- *Analyses, the fast ones only.* H1 (amino-acid content against lifespan, phylogenetic regression with mass and
  with mass and composition, mammals and birds); H11 (composition against the life-history traits, mammals); H7
  (substitution rate against lifespan on sister pairs, the concatenated genes and the four genes that passed); H8
  (repeats and motifs, mammals); H20 (temperature, size and metabolic rate against the substitution rate). Each with
  the command, model, covariates and correction of its own entry, changed in nothing but the species.
- *Not rerun.* RELAX, aBSREL, Coevol and the co-evolution tables (H3, H5, H6, H9, H12 to H19, H21 to H26): for each,
  the species it used that are not of the level high are counted and named; a rerun is proposed only where a
  species of the level unreliable is in a foreground set or is one of fewer than 30 species of a test.
- *Sample size as a covariate.* For H1 and H11 in mammals and in birds, one more model on `not_unreliable`: the
  entry's model with AnAge's sample size added as an ordered term (tiny 0, small 1, medium 2, large 3, huge 4), on
  the species that have one. Reported: the lifespan slope with and without it.
- *Reading, per result that passed its entry's threshold on the original set.* On each restricted set: **holds**
  when the estimate has the original sign and passes the entry's threshold there; **holds, less precise** when it
  has the original sign, fails the threshold, and its 95% interval contains the original estimate (fewer species,
  the same effect); **moves** when its 95% interval excludes the original estimate, whichever way; **reverses**
  when the sign differs and the entry's threshold is passed. A result that did not pass on the original set and
  passes on a restricted one is reported as **new on the restricted set** and is not claimed: these sets were
  chosen after the results were known.
- *Predicted.* Most results hold or hold less precisely on `not_unreliable` and `high`, since 62 and 448 of 2,344
  species leave. The lifespan slopes of birds are the ones expected to move on `well_sampled`, `captive` and `wild`,
  by the bias above.
- *Known weaknesses, said in advance.* A level of high says no rule found a reason, not that the record was checked
  in every respect: 2,265 records were not compared with independent sequences. The level's lifespan reasons and
  the `well_sampled` set overlap. AnAge's sample size is five coarse classes, and FishBase, SeaLifeBase and the
  ringing schemes give none, so `well_sampled` drops every species whose lifespan is theirs. Captive and wild
  records differ by clade as well as by kind.
**Added on 2026-11-19, before any run: which source gives the lifespan.** The maintainer, on the bias above: "maybe
for 3 we switch to more reliable sources than anage". The build takes AnAge's value whenever AnAge has one. Seen
before this addition, from the curation audit: another source reports an observed value at least 1.5 times higher
for 87 species; of 254 birds with both an AnAge value and a ringing record the ringing record is more than 10%
higher in 68, and of 496 fishes with both FishBase's value is more than 10% higher in 129; one such higher value is
itself wrong (FishBase's 54 years under *Sparus aurata*). Fixed now:
- *A second lifespan, `highest_observed`.* Per species, the largest value among the build's lifespan observations
  (`lifespan_observations`) that are not estimated and not rated low or questionable, from any source; the chosen
  value when no other qualifies. A maximum lifespan is a record, so a documented higher value is the better
  estimate unless it is an error; errors are the reason this is a sensitivity run and not a new default.
- *Run.* The analyses above on `not_unreliable`, once more with `highest_observed` in place of the build's lifespan,
  and read by the same four words. Reported beside it: how many species change value, by class and by source, and
  the ratio's distribution.
- *What it decides.* Nothing by itself. If results move, the choice of rule for the build (AnAge first, as now;
  the highest observed; a source per class) goes to the maintainer with these numbers.
**Code for `highest_observed`, and what it decided** (2026-11-20, before any analysis ran on it).
`longevity/highest_observed.py` writes a copy of the build with the replaced lifespan, so an analysis run on the copy
changes in nothing else; `tests/test_highest_observed.py`, 3 tests, each of 14 deliberate breaks caught. Decided where
the text above leaves room: a qualifying observation is of record type captive, wild or observed and of quality high
or acceptable, so the 2015 database's own compiled values (an older AnAge, no quality) never replace a value; a tie
keeps the build's value; AnAge's sample size is dropped where the new value is another source's. Seen when the
command was tried on the build, before any analysis: 289 of 2,344 species get a higher value (125 ray-finned fishes
and 27 sharks and rays from FishBase; 126 birds from the ringing schemes; 6 others from SeaLifeBase), the ratio's
median 1.24 and largest 5.25 (*Lutjanus erythropterus*, 8 against 42 years). The audit read one of the top values as
an error (*Sparus aurata*, 11 against 54), and it is in this set: the run is a sensitivity, as said.
**Code** (2026-11-19, before any run). `src/mitoage_pipelines/longevity/species_sets.py` defines the five sets from the
build's species table. A run is restricted by the environment variable `MITOAGE_SPECIES_SET`: `aa_pgls.default_dataset`,
which every analysis of lifespan calls for its species, then keeps only the species of the named set, so an analysis
changes in nothing but its species; unset, every analysis is that of its own entry. Each restricted run writes into a
data root of its own (`MITOAGE_PIPELINES_DATA`), so no earlier table is overwritten. A species table without
`data_confidence` is refused, a species with no level is in no set, and an analysis that selects its species without
their slug is refused. `tests/test_species_sets.py`: 5 tests, each of 13 deliberate breaks of the code caught. The
sample-size covariate and the readout that sets each result beside its original are not written yet.
**Result** (2026-11-20; Slurm jobs 6538985, tasks 0 to 5, from 88b06f7, and 6540805 from 2650b4b, on the build at
Hugging Face revision 446dad6; every step of every task ended without an error; read by `longevity set-readout` at
b44a66c, 4 tests, each of 12 deliberate breaks caught; tables in `results/h27/`). The sample-size covariate is not
written and not run.
- *Not assessed, said first.* The set did not reach two of the analyses, whose tables are the same bytes on every
  set: H2's column screen (`longevity pgls`, which takes its species from the dated tree's quality list) and H20's
  sister-pair fits (`longevity thermal`, which calls the default dataset for one comparison only). Their 57 passing rows in
  `readout.tsv` say "holds" and mean nothing, and are left out of every count below; H2 and H20 are **not assessed** on any set. `captive` and `wild` were
  fixed for mammals and birds with at least 50 species: captive fishes (8 species) and wild reptiles (4) are in the
  tables and are not read; the three "reverses" of the readout are all captive fishes.
- *Species.* Mammals 604 on all, 599 not unreliable, 527 high, 400 well sampled, 545 captive, 53 wild. Birds 382,
  369, 332, 291, 104, 241. Ray-finned fishes 235, 232, 184, 226. Reptiles 97, 95, 85, 65.
- *`not_unreliable`.* Of the 1,368 results of H1, H7 and H11 that passed on all species, 1,309 hold and 59
  hold less precisely; none moves. All four gene results of H7 hold.
- *`high`.* 1,031 hold, 337 hold less precisely, none moves. H1 in mammals: 15 of 20 passing rows hold, 5 less
  precisely; the isoleucine slope is -0.242 against -0.257 on all. H7: mammalian ND5 and ND2 hold; ND6 across all
  species (-0.087 against -0.196, q 0.54) and avian CYTB (q 0.64) hold less precisely.
- *`highest_observed`.* 1,302 hold, 66 hold less precisely, none moves: taking the highest observed lifespan for
  289 species changes no reading that passed.
- *`well_sampled`.* 686 hold, 681 hold less precisely, 1 moves (alanine in the mammals' mixed model with
  composition, -0.041 against -0.174). On a third fewer species the estimates stay within their intervals of the
  originals and half lose their threshold.
- *`captive` (mammals and birds).* H1 in mammals: 19 of 20 hold. H7's mammalian ND2 moves (-0.004 against -0.017).
- *`wild` (mammals and birds).* In mammals, 53 species: none of H1's 20 passing rows holds, 8 move and 4 are
  unclear; the phylogenetic slopes of isoleucine (+0.160 against -0.257) and threonine (-0.073 against +0.183) have
  the other sign without passing. H11: of 1,327 passing rows about 90 hold, 1,099 hold less precisely, 112 move.
  The entry predicted that the birds' lifespan slopes would move on these sets; in birds nothing passed on all
  species, so that prediction could not be tested, and it is the mammals' results that move on wild records.
- *New on a restricted set, not claimed:* 2 amino acids in reptiles on `high`; 7 to 99 rows of H11 per set; 1 to 5
  rows of H7.
**Reading.** The results of H1, H7 and H11 do not depend on the records the audit found wrong, on the records with
any caveat, or on which source gives the lifespan: nothing that passed moves on `not_unreliable`, `high` or
`highest_observed`. The one set on which they do not stand is mammals with wild lifespan records, 53 species
against 545 captive. That is a different sample of mammals as much as a different kind of record, and the two
cannot be told apart here; the entry said so in advance.
**Extension: the two analyses the set did not reach** (2026-11-20, fixed before the code; no table of H2 or H20 on
a restricted set exists yet). Seen before: everything in the result above.
- *H2's column screen.* `longevity dated-tree` writes the list of species with a usable lifespan record, which
  `longevity pgls` (and the Coevol inputs) take their species from. With a set active, that list holds only the
  species of the set as well. Read on every set, mammals and birds, with H27's reading rule.
- *H20's sister-pair fits.* With a set active, `longevity thermal` forms its pairs from the species of the set only.
  H20 does not use lifespan, so only the two sets defined by the level are read for it, `not_unreliable` and
  `high`; on the sets defined by the lifespan record it is "not applicable", not "not assessed".
- *Not rerun.* The other analyses' tables on each set must come out the same bytes as in the result above, since
  their species do not change; the job says so or fails.
**Extension's result** (2026-11-20; Slurm array 6555962 from 1e65117, seven tasks into run folders of their own,
none with a failed step; check job 6559666: the pair-fit test passes on the build and fails on a copy with the
restriction removed; tables in `results/h27-extension/`).
- *Nothing else moved.* The unrestricted run is the same bytes as the earlier one in every table; on every set only
  `pgls.<class>.tsv` and the three pair tables of `longevity thermal` differ, as said in advance.
- *H2's column screen.* Two columns passed in mammals on all species (COX1 499 isoleucine, +0.146, q = 0.049; ND5
  316 glutamine, -0.213, q = 0.005). `not_unreliable`: ND5 holds (-0.219, q = 0.003), COX1 holds less precisely
  (+0.146, q = 0.050). `high`: both hold less precisely (+0.123, q = 0.21; -0.156, q = 0.45). `highest_observed`:
  ND5 holds, COX1 less precisely (q = 0.057). `well_sampled`: COX1 less precisely, ND5 not assessed; `captive`: both
  hold; `wild`: not assessed. None moves. In birds no column passes on all species or on any set. Four columns pass
  on a restricted set only (one on `well_sampled`, three on `captive`); found on subsets after the fact, they are
  listed and not read.
- *H20's pair fits.* Of the 17 rows that passed on all species: `not_unreliable` 12 hold and 5 hold less
  precisely; `high` 12 hold, 4 hold less precisely, 1 is not assessed; none moves. The other sets are defined by
  the lifespan record, which H20 does not use: not applicable.
- *Counts with the two analyses in.* Of 1,425 passing results: `not_unreliable` 1,360 hold, 65 less precisely, none
  moves; `high` 1,081, 343, none moves, 1 not assessed; `highest_observed` 1,353, 72, none moves.
So the reading above stands for H2 and H20 as well: no result changes direction without the unreliable records or
on high-confidence records; H2's two columns, already at the threshold, lose it on the 527 high-confidence mammals.

**Sample size as a covariate: what the entry left open, fixed on 2026-11-21 before the run.** `mitoage longevity
sample-size-covariate` (`longevity/sample_size.py`). Each model is fitted twice on the same species, those of
`not_unreliable` whose lifespan is AnAge's with a sample size: without the term and with it, so the two slopes differ
by the term alone. H1: both models (mass; mass and fourfold composition), by phylogenetic regression and by the mixed
model. H11: every model of the entry that holds a lifespan term (lifespan; lifespan and mass; lifespan, a competitor
and mass; all five), on Pagel's lambda and with lambda fixed at 1; the two slow-fast axes are not refitted, since
lifespan is inside the axis and not a term beside it. q is taken within the entries' own groups, separately without
and with the term. A result that passes without the term is read by the rule above, the fit without the term as the
original: holds, holds less precisely, moves, reverses. Also reported: the median of slope with over slope without
among the passing results, and how often the term itself has p < 0.05. A result that passes only with the term is
counted and not claimed. Not predicted in the entry; said here before the run: the slopes are expected to shrink
little in mammals, and more in birds, where the excess of lifespan over mass with sample size was seen.

**Sample size as a covariate: result** (2026-11-21; Slurm job 6564454 from 7dac172, on the data root of
`not_unreliable` of the extension's run; tables in `results/h27-sample-size/`). Species with an AnAge sample size:
mammals 599 of 599 in H1's test (tiny 43, small 156, medium 269, large 130, huge 1), birds 361 of 369 (4, 66, 218,
73, 0); H11's models 437 to 554 mammals and 103 to 265 birds. Of 2,052 fits, 206 pass their threshold without the
term: with it 187 hold and 19 hold less precisely (each with q between 0.050 and 0.091); none moves and none
reverses.

| Analysis | Class | Pass without | Hold | Less precise | Moves | Slope with / without, median (range) |
|---|---|---|---|---|---|---|
| H1 | mammals | 19 | 16 | 3 | 0 | 1.00 (0.84 to 1.19) |
| H1 | birds | 0 | | | | |
| H11 | mammals | 176 | 165 | 11 | 0 | 1.06 (0.88 to 1.43) |
| H11 | birds | 11 | 6 | 5 | 0 | 1.05 (0.95 to 1.16) |

The term itself has p < 0.05 in 87 of the 1,972 fits that report it, against 99 expected at that level by chance,
with both signs (37 positive and 50 negative: mammals 29 and 27, birds 8 and 23). 53 results pass only with the
term (52 in H11's mammals, 1 in H1's birds); they are counted, not claimed.
*Reading.* The lifespan slopes of H1 and H11 do not carry AnAge's sample size: adding it leaves them where they
were, and where they change they grow slightly more often than they shrink. What was said before the run, that the
slopes would shrink more in birds, is not what came out: in birds few results pass at all (none in H1, 11 in H11,
all with lambda fixed at 1), and those keep their size. *Limits.* The sample size is five coarse classes, entered as
equally spaced; almost no species is "huge" and few birds are "tiny", so the term is estimated from the middle
classes. It says nothing of lifespans from other sources, which have no sample size, and nothing of fishes or
reptiles, which the entry did not include.

## H28. The cybrids by a learned binding-energy model: StaB-ddG on the hosts' complexes (code written, not run; H16's to H18's open item)

**State: proposed on 2026-11-21 with ProteinMPNN; the maintainer agreed on 2026-11-22 and asked for a check of the
2026 literature for a better model first. The check changed the model to StaB-ddG, on 2026-11-22, before any code
and any score. The code followed the same day (`mitoport cybrid-stab`); not run.** The protocol below is fixed with
the code, before the first score; any change after this date is said as a change.
*Made exact with the code, 2026-11-22, before any score:* the system is H18's (not every chain of the entry); the
contacts and controls are H17's (the human structures' sites, not a new map on the host's structure); a break in a
chain is kept as masked places (chains are not renumbered); the score is the mean over the three seeds; the whole
run, not seed 0 alone, is repeated for the byte comparison.
*Added on 2026-11-22, after two tasks had started (their first lines of output were seen, no readout was computed)
and restarted from the beginning:* the maintainer asked for a verdict per ported protein, not a statement about
contacts ("porting this protein will have no effect but this one will break efficiency"). Two outputs are added and
the two scores read by the rule stay as fixed. (1) Per protein: for each donor and each mtDNA-encoded protein ported
alone, `stab_fold`, the model's predicted change in the stability of the whole system (substitutions anywhere in the
protein), and `stab_bind`, its change of binding to the rest (stab_genes.tsv). No rule is fixed for them: the
measurements are per complex, so only cytochrome b, one protein and one complex, is tested directly; for the other
twelve the numbers are predictions the cybrids can neither confirm nor refute one by one. (2) Per complex,
`stab_fold_all`, the same stability term for all substitutions, reported beside the scores and claimed as nothing.
What ground a per-protein number has is H29's question.
**Question.** Counts (H16), co-evolved pairs and ESM-2 (H17) and a force-field energy (H18) do not find the cybrids'
defects better than divergence. ESM-2 never saw the partner in space, and the energy does not repeat between builds
of one structure. Does a model trained to predict how substitutions change the binding energy of two partners, given
their structure, find the donor's mitochondrial subunits binding the host's nuclear subunits worse where the
measured complexes fail?
- *Why this model (literature read on 2026-11-22, a search and not a review).* StaB-ddG (Deng et al., ICML 2025,
  arXiv 2507.05502) is ProteinMPNN fine-tuned first on folding stabilities and then on measured binding energy
  changes; it predicts the change of binding energy between two groups of chains for one substitution or many at
  once, and its authors report it level with FoldX and Flex ddG on held-out interfaces at a fraction of the cost.
  Against plain ProteinMPNN it answers the question asked here (binding of one partner to the other, not how well a
  residue fits), takes all of a donor's substitutions together, and its repository carries the plain ProteinMPNN
  weights as well, so the first proposal becomes a control arm of the same code. Looked at and not taken: BioEmu
  (ensembles of single chains; its 2026 follow-up says complexes are not modelled); a 2026 preprint on inverse
  folding from surface and backbone (reports sequence recovery, no mutation effects, not peer reviewed); supervised
  predictors of binding energy changes from 2025 and 2026 (USP-ddG, Prot-ΔΔG, SimBinder-IF, ProtMutMap), which are
  trained on the same measured set of soluble complexes and were not shown to be better outside it. No work was found
  that scores cybrids or mitonuclear mismatch with any such model. What a wider search would find is not known.
- *Model.* StaB-ddG (github.com/LDeng0205/StaB-ddG at commit b19a5d7, MIT licence), its final weights `stabddg.pt`,
  the defaults of its prediction script: backbone noise 0.1, 20 Monte Carlo samples. The prediction is the model's
  folding-energy change of the whole complex minus those of the two partners taken alone, in kcal/mol, positive for
  weaker binding. Each set of substitutions is scored on its own, with seeds 0, 1 and 2; a score is the mean over
  the three.
- *Structures.* H18's: human 5XTD (I), 5XTE (III), 5Z62 (IV), 8H9S (V); mouse 8OM1 (I), 7O3H (III), 7O37 (IV). The
  backbone as deposited (N, CA, C, O); the system is H18's: the complex's mtDNA-encoded chains, every copy, and each
  protein chain with a heavy atom within 4.5 Å of one of them, without the mtDNA-encoded chains of another complex.
  Partner one is the mtDNA-encoded chains, partner two the other chains of the system. A break in a chain's
  numbering is kept as that many masked places, at most 33 (the model does not tell longer distances apart).
  Nothing is built or relaxed, so the step that made H18's energies differ between builds is absent.
- *Scores,* per (host, donor, complex), named as in H17. `stab_all`: the predicted change of binding energy for all
  of the donor's substitutions at resolved positions of the mtDNA-encoded chains, given together. `stab_contact` and
  `stab_control`: the same for the substitutions at the contact residues only and at their burial-matched controls
  only; contacts and controls are H17's, the sites of the human structures carried to the host's protein by
  alignment. Deletions, insertions and unresolved positions are counted and not modelled, as in H18.
- *Readouts and rule: H16's and H17's, unchanged* (`mitoport cybrid-readout`). The two scores read are `stab_all`
  and `stab_contact`. A score counts if it beats divergence (`d_all`) on all three: T1 a more negative Spearman than
  −0.36 over the 13 stated activities; T2 an area above 0.74 for 9 reduced against 20 not; T3 the most reduced
  complex named in at least 5 of the 6 pairs. Two of three is reported as partial and claimed as nothing. No
  significance is claimed: the points are few and not independent.
- *Controls.* `stab_control` beside the scores (a score that only repeats it is divergence again). The same run with
  the repository's two earlier weights: `proteinmpnn.pt`, plain ProteinMPNN (the first proposal), and
  `stability_finetuned.pt`. The model draws random noise: seeds 0, 1 and 2, and the spread over them reported beside
  every score; a readout whose verdict changes with the seed is said as not settled. The whole run twice, the two
  sets of tables compared byte for byte. An empty set of substitutions must give exactly 0, and a table where it
  does not is refused.
- *Known weaknesses, said in advance.* The model learned from small soluble domains and from measured soluble
  complexes, mostly with one to a few substitutions; here it scores membrane complexes of up to 45 chains with
  tens to hundreds of substitutions at once, far outside what it was tested on. Its "partners alone" are the same
  backbones taken apart, never relaxed. Its parent was trained on the Protein Data Bank, which holds these complexes
  or close relatives: it may prefer the host's residue everywhere, and then the score is a count of substitutions.
  It sees no cofactor, lipid or water, and it predicts a binding energy, not a rate of electron transfer.
- *If it fails:* with H16 to H18, four kinds of score would then agree that the published cybrids cannot be told from
  divergence with 30 points, and what remains is dynamics at a few named sites (cytochrome b near the Rieske protein),
  not another score over all contacts.

## H29. Does the model tell harmful from harmless substitutions within a mitochondrial protein? Disease variants against common ones (code written, not run; H28's check)

**State: proposed on 2026-11-22; the code followed the same day (`mitoport variant-stab`), not run.** The maintainer
asked whether published work can validate the method; told of this check, he had not answered when it was written.
It is written before any score, and runs after H28 unless he says no.
*Made exact with the code, before any score:* the control weights run with seed 0 only; the nuclear partner alone
holds no substitution, so its term is exactly 0 and is not computed; variants scored in one call share their
draws, as in the model's own script (four at a time on CPU; one at a time on the 11 GB cards it ran on, said with
the results); the site's variability is counted among the build's Mammalia.
Counted from the sources as fetched on 2026-11-22, with no score involved: 41 harmful and 331 harmless
variant-and-gene rows; none harmful in COX1, COX3, CYTB or ND2, one in COX2, so the test speaks for complexes I
and V and says nothing for III and almost nothing for IV. Only ATP6 and ND6 have 10 harmful variants.
**Question.** H28 gives a number per ported protein, and the cybrids cannot test it protein by protein. The one data
set that is per protein and not about contacts is human variation: does StaB-ddG's stability term separate
confirmed disease variants of the 13 mtDNA-encoded proteins from variants common in healthy people?
- *Variants.* Single-nucleotide missense variants of the 13 protein genes, the amino-acid change worked out here
  from the human reference (NC_012920.1) and not taken from either source. Harmful: ClinVar's mitochondrial records
  classified pathogenic or likely pathogenic and reviewed by its expert panel. Harmless: gnomAD v3.1's variants that
  pass its filters and are homoplasmic in at least 1 of 1,000 people (`AF_hom` of 0.001 or more). A variant in both
  lists, or pathogenic in ClinVar at any review level, is left out of the harmless ones and counted. The fetched
  files' SHA-256 go into the manifest. (Changed on 2026-11-22, before any code: MITOMAP's confirmed list was the
  first choice, but MITOMAP serves its pages to a browser only. ClinVar holds 93 such single-nucleotide records
  over the whole genome, tRNA genes included; how many are missense is counted by the code and may be few.)
- *Score.* Each variant alone on the human structures of H28 (5XTD, 5XTE, 5Z62, 8H9S), every copy of the chain, the
  weights `stabddg`, seeds 0, 1 and 2: `fold`, the predicted change in the stability of the whole system, and beside
  it the change of binding. A variant at an unresolved position is not assessed and counted.
- *Baselines the model must beat.* The BLOSUM62 score of the substitution, and the variability of the site among the
  build's mammals (how many different amino acids its alignment column holds). A model that only says "conserved
  sites matter" adds nothing to these.
- *Rule.* Over all proteins together, the area under the ROC curve for harmful against harmless is at least 0.70,
  and above each baseline's on the same variants, the difference's 95% bootstrap interval (resampling variants)
  lying above 0. Both, or the model does not count for per-protein verdicts. Per protein, the area is reported where
  the protein has at least 10 variants of each kind, and claimed as nothing.
- *Controls.* The weights `proteinmpnn` and `stability_finetuned`; the split by contact and not contact, reported;
  the run twice, compared byte for byte.
- *Known weaknesses, said in advance.* Disease variants were found partly because they sit at conserved sites, which
  favours every score, baselines included. Common variants are a human's own, so the structure's sequence is close
  to them. The model's parent saw these structures in training and may simply prefer the deposited residue. And the
  test is of one substitution at a time, while a donor carries tens to hundreds: passing here does not show the
  model adds them up rightly.
- *What follows.* Pass: H28's per-protein numbers may be shown as predictions with this evidence named. Fail: no
  per-protein verdict is printed from this model, whatever H28's complex-level result.

## H30. Two healthy mitochondrial genomes in one cell: are complexes assembled from a mix predicted to be less stable? (agreed, after H28; no code)

**State: agreed on 2026-11-22 ("Let's test second one after we get cybrid results"), no code, not run.** Written
before any score.
**Background.** Mice heteroplasmic for two normal mouse mtDNAs, each healthy alone on the same nucleus, are not
healthy as a mixture: reduced activity and impaired cognition for NZB with 129 (Sharpley et al., Cell 2012), heart
disease, frailty and a shorter life for NZB with C57 (Lechuga-Vieco et al., Circulation 2022), in the tissues that
cannot get rid of one type. Three explanations are on offer and none is shown: the cell cannot tune itself to two
genomes at once; complexes assembled from subunits of both fit themselves worse; the cell pays for trying to sort
them. Only the second is a statement about structure, and only it can be put to this model.
**Question.** For a complex with two or more mtDNA-encoded subunits, does StaB-ddG predict a mixed set (some
subunits from one genome, the others from the other) to be less stable than both pure sets?
- *Genomes.* C57BL/6 (the mouse reference) and NZB, and 129 as the second pair; the accessions are fixed with the
  code, not looked up yet. The substitutions are those between the two genomes' proteins; if a complex's subunits do
  not differ between the two, the complex is left out and said.
- *Structures and model.* H28's mouse structures (8OM1 for I, 7O37 for IV; complex III has one mtDNA-encoded
  subunit and nothing to mix; no mouse structure of V is in H28, so V is not assessed unless one is added before the
  run), the system, the weights `stabddg` and seeds 0, 1 and 2 of H28. The host sequence is the structure's
  (C57-like); a set is the list of subunits that take the other genome's residues.
- *Score.* `fold` of each set: the pure other genome (all subunits), and every mixed set (complex I: the 126 proper
  subsets of its seven subunits that are neither empty nor all). `mix_excess` of a mixed set is its `fold` minus the
  sum of the `fold` of its subunits taken one at a time: what the model says the combination costs beyond its parts.
  Mixing that is harmless by the model gives an excess near 0.
- *Rule.* Explanation 2 is supported, as a prediction, only if the pure other genome is not itself predicted worse
  than the largest mixed set, and the median `mix_excess` over the mixed sets is above 0 by more than twice the
  spread over the seeds, in the same direction for both genome pairs. Otherwise it is not supported by this model.
- *Known weaknesses, said in advance.* The genomes differ by about a dozen amino acids, so every number will be
  small and near the seeds' noise. The model's "excess" is non-additivity in its own likelihood, which was never
  tested against measured non-additivity. There is no measurement of a mixed complex anywhere: the animals'
  phenotype is predicted alike by all three explanations, so even a clear result makes explanation 2 plausible or
  less plausible, never shown. And whether H29 passes decides how much any of this is worth.
