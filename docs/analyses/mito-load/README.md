# Design tally for the mitochondrial-load adaptations

`tally.tsv` is what `docs/MITO_LOAD_ADAPTATIONS.md` and the design of H21-H23 read: for each candidate adaptation, how
many species of the build carry it, how many are unknown, how often it changes on the dated tree, and how many
trait-bearing species have nuclear partners. It was computed on 2026-11-16, before H21-H23 were written (the bat hibernation row before H24 and H25, the mammal hibernation row after them, from a table that carries no lifespan), and tests
nothing: no rate, selection, composition or coupling statistic was computed for any trait.

```sh
uv run python docs/analyses/mito-load/tally.py <dna-seq/mitoage checkout> data/input/build2/fishbase/v26.06 docs/analyses/mito-load/tally.tsv \
    data/input/hibernation/wilkinson_adams_2019_bats.txt data/input/build2/taxdump/2026-11-13/names.dmp \
    data/input/hibernation/heldstab_2018_supplementary_material.xlsx
```

Inputs:
- dna-seq/mitoage at 03c0949 (the data of longevity-genie/mitoage 0f4e970; unchanged at bf2722b):
  `public/data/species.csv` (sha256 74f07932bb64..., 2,344 rows, retired records left out) and
  `src/data/generated/phylogeny.json` (sha256 bfcd086948b3..., opentree16.1, 2,314 tips);
- `public/data/mitoport/partners_species.parquet` of dna-seq/mitoage at bf2722b (D20), the species with vertebrate partners of complex IV
  (H12; Hugging Face longevity-genie/mitoage-pipelines 126d237): 177 mammals, 78 birds, 61 ray-finned fishes;
- FishBase v26.06, `species`, `synonyms` and `stocks`, md5 equal to `build2/pins.py` (`download fishbase-temperature`).
  The 980 fishes are matched by H20's rule (`longevity.thermal.fishbase_codes`): 955 by exact name, 25 by synonym,
  none ambiguous or in conflict with the build's FishBase record;
- the bat table of Wilkinson & Adams 2019 (`download hibernation`, SHA-256 pinned in `adaptations/hibernation.py`),
  its `Hibernation` column only, matched to the 44 bats by exact name (42) or by a name NCBI's taxonomy of 2026-11-13
  gives to the species' own taxid (2: *Eptesicus nilssonii*, *Rousettus lanosus*);
- the mammal sheet of Heldstab et al. 2018 (Dryad, md5 pinned, saved by hand), its `hibernation` column only, combined
  with the bat table by `hibernation.combine` (one source: its state; both alike: that state; both differing or
  neither: unknown).

Columns: `yes`, `no`, `unknown` (a missing FishBase field, and FishBase's default `WaterAssumed` for air breathing,
which means that no record exists); `families_yes`; `trait_clades`, the maximal clades of the dated tree whose known
tips are all trait-bearing; `min_changes`, the fewest changes of state on the tree (Sankoff, gains and losses
together, unknown tips free). The two bracket how many times a trait arose within the build's sample; neither counts
origins in nature. `yes_with_partners` and `scope_with_partners` count species in the partner set above.
