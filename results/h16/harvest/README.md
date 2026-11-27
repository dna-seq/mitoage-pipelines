# The cybrid harvest, first pass (2026-11-27)

The rules are in `docs/CYBRID_HARVEST.md` and were committed before any paper was read. This folder is the
decision log that document requires, and it is **not finished**: the searches below are the first pass.

## Searched

PubMed, through its API, on 2026-11-27: `xenomitochondrial` (20 records, the whole of that term's literature),
`cybrid interspecific mitochondria` (7), `transmitochondrial mice xenomitochondrial` (4),
`xenomitochondrial cybrid respiratory chain activity` (3). Queries are recorded per row in
`candidates.tsv`. Searches still to run: the transmitochondrial and conplastic animal literature by its own
terms, forward citations of H16's five source papers, and the reviews' reference lists.

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

## Where it stands against the bar

The bar fixed in advance: at least 10 new points from at least two new papers to count as useful, reduced
points outside complex I to change H16's reading, a third host with one reduced and one not-reduced point to
open a new host.

After the first pass: **1 point, from 1 paper.** Against the bar this is nothing yet, and it is recorded as
nothing.

What the first pass does establish is where the remaining room is. The whole `xenomitochondrial` literature
is 20 records, and five of them are already H16's sources. The candidates that could still carry points are
few and named in `candidates.tsv`, and the most promising one is out of reach:

- **The hare cybrids (10.1016/j.yexcr.2025.114466, 2025)** are reciprocal brown-hare and mountain-hare
  cybrids with respirometry: a third and a fourth host species, and the only candidate that could open a new
  host by itself. The publisher refuses automated access (403 at the article page, and the DOI redirects to
  the same publisher). That gate is there on purpose and was not worked around. It is recorded as `unread`,
  which is not a negative.
- Four more are `unread` for want of a readable copy (10.1073/pnas.0303184101, 10.1023/B:JOBB.0000041778.84464.16,
  10.1074/jbc.274.23.16188, 10.1080/10495398.2013.841709), and two are `unscreened` pending the amino-acid
  check that rule 3 requires of a within-species conplastic line.

So the direction is not yet answerable either way, and the reason is access rather than absence. The
judgement of whether to ask the maintainer for copies of the five unread papers belongs to him and is not
made here.
