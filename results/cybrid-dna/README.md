# Nucleotide baselines beside the cybrid scores (exploratory)

From the cluster's login node (counts and alignments of 15 mitochondrial genomes, seconds), `mitoage mitoport
cybrid-dna` run from a fresh clone at commit 495c30c on the records of `download cybrid` and H16's `scores.tsv`. The
tests (21 passed) and the breaks of each guard (15 of 15 caught) ran on the same commit. These scores were chosen
after H28's readout had been read: nothing here passes or fails a rule. What they are and why: `docs/HYPOTHESES.md`,
H28's last item, and the module's own text.

| file | one row is |
| --- | --- |
| `dna_scores.tsv` | one (host, donor, complex): the nucleotide differences in the complex's genes, all, silent and amino-acid-changing, and the pair's rRNA, tRNA, humanin, MOTS-c, control-region and conserved-block differences |
| `dna_calibration.tsv` | one measured (host, donor, complex) with its status, its activity and the scores |
| `dna_summary.tsv` | one readout (T1 to T4 of H16) of one score |
| `dna_correlations.tsv` | the rank correlations between the scores across the 13 donors |

What it shows: nucleotide scores follow the measured activities at least as well as the amino acid counts (T2 0.73
to 0.87 against 0.74 for `d_all`); silent differences do as well as amino-acid-changing ones (0.82 against 0.81);
the coding, rRNA and tRNA scores correlate 0.84 to 0.98 with each other. The points mostly carry the donor's
distance, and they do not rank one region above another.
