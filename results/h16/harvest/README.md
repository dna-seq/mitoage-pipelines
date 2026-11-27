# The cybrid harvest (2026-11-27)

The rules are in `docs/CYBRID_HARVEST.md` and were committed before any paper was read. This folder is the
decision log that document requires, and it is **not finished**.

## Searched

PubMed, through its API, on 2026-11-27. First pass: `xenomitochondrial` (20 records, the whole of that term's
literature), `cybrid interspecific mitochondria` (7), `transmitochondrial mice xenomitochondrial` (4),
`xenomitochondrial cybrid respiratory chain activity` (3). Second pass, the transmitochondrial and conplastic
animal literature by its own terms: `transmitochondrial mice respiratory` (21), `conplastic mice mitochondrial
respiration` (3). Queries are recorded per row in `candidates.tsv`.

Searches still to run: forward citations of H16's five source papers, and the reviews' reference lists.

## Files

- `candidates.tsv`: every record the searches returned, with one of `already_used`, `included`, `excluded`
  (with the rule it fails), `unread` (with where the reading stopped) or `unscreened`.
- `new_points.tsv`: the points extracted so far, in the columns of
  `src/mitoage_pipelines/mitoport/cybrid_measurements.tsv`.

## A departure from the harvest's own document, and why

The document says extracted rows go into `cybrid_measurements.tsv`. They are being held in
`new_points.tsv` instead until the harvest closes, because that table is the input of H16 and of every
cybrid test: adding rows to it while the harvest runs would change those inputs silently, in pieces, and the
document's own rule is that no test is rerun until the harvest is finished and its table committed. The
rows are in the target's columns and move across unchanged when it closes. Nothing else in the document is
departed from.

## Two places where the rules, fixed in advance, decided against a point

Both are recorded rather than argued around, because the document's reason for fixing the rules first was
exactly that a rule bent with the result in view finds what it hoped for. Neither was changed here. Both are
for the maintainer to rule on after the harvest closes, not during it.

1. **Activity, not abundance.** The outcome rule speaks of an *activity* stated against its own control. The
   four-donor cybrid panel of 10.1371/journal.pone.0239804 states significant per-complex differences in
   **protein abundance** by BN-PAGE (complex I down in the *Mus caroli* and *Mus pahari* cybrids, complex IV
   down in *Mus pahari*, complex III down in *Mus caroli*), and no significant change in the respiration it
   states in text, which is basal and maximal and so cannot see one complex alone. Widening "activity" to
   "abundance" now would add four points, three of them outside complex I, which is precisely what the bar
   below names as the thing that would change H16's reading. That is the move the document forbids, so the
   paper is `excluded` and the question is left standing.
2. **Rule 3 cannot be executed as written for a within-species conplastic line.** It says the amino-acid
   difference is "checked against the build, not taken from the paper", but the build carries one
   mitochondrial genome per species, so it cannot adjudicate two strains of *Mus musculus*. The independent
   check for those is the strain mtDNA in GenBank. Until that check is run, a conplastic candidate stays
   `unread` or `unscreened` and no point is taken from it. No conplastic paper has yielded a point so far, so
   nothing turns on this yet.

## What the second pass found

**One new point, and it contradicts two old ones.** The xenomitochondrial XM mouse (C57BL/6J nuclear genome,
homoplasmic *Mus terricolor* mtDNA) is stated in 10.3390/antiox13020187 to have complex I-linked ADP-stimulated
respiration of **cardiac** mitochondria significantly decreased by about 40%. The same host, donor and complex
already stand in the table as **not** reduced twice over: at 83% of control in the cell cybrid, and at 90% in
skeletal-muscle mitochondria of the in-vivo line (10.1016/j.mito.2010.07.003, this harvest's first point).

The harvest's rule for a disagreement is that both rows are kept with their sources and the point is excluded
from testing. That rule is applied. But the disagreement here does not look like two laboratories differing:
it tracks the **tissue** — cells and skeletal muscle on one side, heart on the other. The table has no tissue
column, and H16 treats (host, donor, complex) as one point. If complex I in this pairing is reduced in heart
and not in muscle, then that point is not well defined without a tissue, and neither is any other row in the
table. This is the first thing the harvest has turned up that bears on H16's premise rather than on its
sample size, and it is reported as a finding in its own right.

Two more candidates were set aside, as the first pass's five were, for want of a readable copy: the reciprocal
NOD/ALR conplastic lines (10.1074/jbc.M609367200), whose abstract states that electron transport enzymatic
activities did not differ between the strains and which would therefore be several not-reduced points in both
directions; and the B6-mt(NOD) conplastic line (10.1016/j.freeradbiomed.2016.11.035). Mouse cells carrying rat
mtDNA at 40% respiratory activity (10.1093/hmg/ddl005) are `unread` for the same reason, and would in any case
need per-complex assays.

Seventeen records of the transmitochondrial-mouse literature are excluded as disease models: they are
single-species hybrids carrying a pathogenic or tumour-derived mtDNA variant, which the document's rules put
out from the start as H28's and H29's subject.

## Where it stands against the bar

The bar fixed in advance: at least 10 new points from at least two new papers to count as useful, reduced
points outside complex I to change H16's reading, a third host with one reduced and one not-reduced point to
open a new host.

After two passes: **2 points, from 2 papers**, both of them complex I in a host and donor the table already
has, and the two disagreeing with each other. Against the bar this is nothing, and it is recorded as nothing.

The searchable literature is close to exhausted. The whole `xenomitochondrial` term is 20 records, five of
them already H16's sources; the transmitochondrial literature is almost entirely disease models; the
conplastic literature that passes the design rule is three records. The candidates that could still carry
points are named in `candidates.tsv`, and the ones that matter are out of reach:

- **The hare cybrids (10.1016/j.yexcr.2025.114466, 2025)** are reciprocal brown-hare and mountain-hare
  cybrids with respirometry: a third and a fourth host species, and the only candidate that could open a new
  host by itself. The publisher refuses automated access (403 at the article page, and the DOI redirects to
  the same publisher). That gate is there on purpose and was not worked around. It is recorded as `unread`,
  which is not a negative.
- Seven others are `unread` for want of a readable copy (10.1073/pnas.0303184101,
  10.1023/B:JOBB.0000041778.84464.16, 10.1074/jbc.274.23.16188, 10.1080/10495398.2013.841709,
  10.1074/jbc.M609367200, 10.1016/j.freeradbiomed.2016.11.035, 10.1093/hmg/ddl005).

So the direction is still not answerable either way, and the reason is access rather than absence. The
judgement of whether to ask the maintainer for copies of the unread papers belongs to him; it is desk
question #97.
