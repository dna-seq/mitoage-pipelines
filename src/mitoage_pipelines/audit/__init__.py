"""An audit of the build's mitogenome records (H13 in docs/HYPOTHESES.md): which strand each record deposits,
its gene order, its genetic code, and records that are near-identical to another species' or sit nearer to
another order or class than to their own. Descriptive checks of the data for a person to read; not a test of a
hypothesis about longevity.

  records.py     the RefSeq records of the build's species, from the pinned release, with their inputs' hashes
  strand.py      orientation against the vertebrate convention, beside the GC/AT balance (skews)
  order.py       gene order canonicalised by rotation; arrangements per class, order and family; flags
  code.py        declared against expected genetic code; internal stop codons
  similar.py     near-identical pairs of different species; nearest neighbours outside their order or class
  common.py      table writing (stable order) and warnings aggregated by reason
"""
