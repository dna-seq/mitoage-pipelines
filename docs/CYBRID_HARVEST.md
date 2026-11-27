# The cybrid harvest: rules fixed before anything is collected

The maintainer chose this on 2026-11-26 as the first of four directions, and it is reading, not computing.
Every structural and model-based attempt to say which nuclear gene limits cybrid compatibility (H16, H28,
H29, H34, H35 part 1, H36 steps 1 to 3) has ended in the same place, and the limit named each time was the
data: **29 cybrids, two usable hosts, and a handful of complexes with a stated outcome.** No better model
repairs that. More measured (host, donor, complex) points might.

This document fixes what the harvest counts as a point, how each one is extracted, and what the harvest has
to deliver before any of it is used, all before a single paper is read. It is written the same way a
hypothesis entry is, and for the same reason: a search whose inclusion rule is settled while the results are
in view finds what it hoped for.

## What the harvest is looking for

A **point** is one (host species, donor species, respiratory complex) with a stated outcome, measured in
cells or in an animal whose mitochondrial DNA comes from one species and whose nuclear genome comes from
another. Three designs qualify:

1. **Xenomitochondrial cybrids**: a host cell line emptied of its own mtDNA, repopulated with a donor
   species' mitochondria. This is what the existing 30 points are.
2. **Transmitochondrial mice and other animals**: a conplastic or xenomitochondrial line carried in a whole
   animal, where the mitochondrial genome is from another species or a divergent subspecies and the
   respiratory complexes were assayed.
3. **Conplastic lines within a species** (one species' nucleus with another population's mtDNA) **only when
   the two mitochondrial genomes differ by at least one amino acid in a complex's mitochondrial subunits,**
   which is checked against the build, not taken from the paper. Below that, the line tests something else.

Cell hybrids where both genomes come from the same species and the mtDNA carries a known pathogenic
variant are **out**: that is a disease model, and the variants are H28's and H29's subject.

## The outcome, and what is not one

- **reduced**: the paper states an activity below its own control with the difference called significant, or
  states no respiring cybrid could be recovered for that donor.
- **not reduced**: the paper states an activity and calls it not different from its control.
- **not assessed**: that complex was not assayed, or was assayed only in a linked assay that cannot see the
  complex alone (the murid II+III assay in the existing table is such a case and is marked).

A number without a stated comparison is **not** an outcome. A figure read off by eye is **not** a stated
value: it is recorded as the direction only, with `pct_low` and `pct_high` left empty and the note saying it
was read from a figure. Where two laboratories disagree on the same (host, donor, complex), both rows are
kept with their sources and the point is excluded from testing, as the rat complex III point already is.

## Extraction

Rows go into `src/mitoage_pipelines/mitoport/cybrid_measurements.tsv`, whose columns already exist: `host`,
`donor`, `complex`, `assay`, `outcome`, `pct_low`, `pct_high`, `significant`, `source_doi`, `location`,
`note`. Three rules that the existing rows follow and the new ones must:

- **`location` names the figure or table the value came from**, never the paper as a whole. A row whose
  location cannot be given is not entered.
- **A percentage this repository computed** from the paper's means, rather than one the paper states, says
  so in `note`. The existing murid rows do.
- **Species are named as the build names them**, and a paper's strain or line name is resolved to a species
  in the note, never silently. Where a paper names only a common name with more than one species behind it
  (the 1997 orangutan paper), every candidate species is entered, as those rows already are.

## The decision log

Every candidate paper the search returns is recorded in `results/h16/harvest/candidates.tsv` with its
identifier, what the search was, and one of: `included` with the number of points, or `excluded` with the
rule it fails, in the words of this document. A paper that cannot be read (paywalled with no accessible
full text) is `unread` with where it stopped, never `excluded`: unknown is not a negative. No gate that is
there on purpose is worked around, and nothing is mailed to any author.

## What the harvest must deliver before it changes anything

Fixed here, before the first search, so that a thin harvest is reported as thin:

- The harvest **counts as useful** only if it adds at least 10 new points with a stated outcome, from at
  least two source papers this repository does not already use.
- It **changes H16's reading** only if it adds reduced points outside complex I. Every one of H16's failures
  to name the failing complex comes from the raw counts always naming complex I, the largest; more complex I
  points cannot fix that, and more non-I points might.
- It **opens a new host** only if a third host species arrives with at least one reduced and one not-reduced
  point, which is what the two-host limit means in practice.
- Anything short of these is written up as what it is: the harvest ran, here is what exists, and the data
  for this question remains as thin as it was. **A null harvest is a result and is reported.**

No test is rerun on the enlarged table until the harvest is finished and its table committed, so that the
inclusion of a point is never decided with its effect on a test in view.
