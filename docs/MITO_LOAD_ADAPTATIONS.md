# Adaptations that load mitochondria, and what the build can test about them

The maintainer, 2026-11-16: "brainstorm about potential biological adaptations that put load in mitochondria evolutions
in various species, so we can actually try to classify species that have such adaptations and then check if something
can be found in coevolutions".

This document ranks the candidates by what the MitoAge build can test, with five things settled for each: why the
adaptation should load mitochondria, how species are classified, how many species of the build it covers, how many
independent origins it has on the tree, and what signal to expect with which tool. The shortlist agreed with the
maintainer the same day ("Fish ecology + mammals") is written as H21, H22 and H23 in `docs/HYPOTHESES.md`; his follow-up
("check for hibernation idea as it is connected with longevity", "bats look interesting") as H24 and H25. Nothing
here has been run. Every citation was taken from a search result and its title checked; none is from memory.

Temperature, size and metabolic rate as continuous traits are H20 and are not repeated here.

## What decides testability

**Independent origins, not species.** A trait that arose once is a sample of one, whatever its species count: bats
are 44 species and one origin, so any difference between bats and other mammals is confounded with everything else
that happened on that one branch. The tally counts origins two ways on the build's dated tree:
- `min_changes`, the fewest changes of state (gains and losses together) that explain the tips;
- `trait_clades`, the maximal clades whose known tips all carry the trait.

Both count within the build's sample. Neither counts origins in nature, and either can be inflated by a species
whose record or classification is wrong.

**Nuclear partners, for anything about coevolution.** The coevolution tables (H5, H14, H19; on the site under D20)
need the nuclear subunits, and those exist for 177 mammals, 78 birds and 61 ray-finned fishes (complex IV's vertebrate
set of H12). H19's vertebrate alignments, once run, will hold a median of about 120 ray-finned fishes per gene. Only
four candidates put more than 8 trait-bearing species into the partner set: marine mammals (26), diadromous fishes
(22), bats (13) and, probably, hibernating mammals (to be counted with the mammal-wide source below; 6 hibernating bats
are in the set). Every other candidate can be tested on
the mitochondrial genes only: dN/dS, RELAX, aBSREL, substitution rate, composition.

**Unknown is not no.** FishBase leaves migration (`AnaCat`) empty for 403 of 893 ray-finned fishes and depth for 322.
Its air-breathing field defaults to `WaterAssumed`, which means no record exists. Read as unknown, air breathing has
32 known yes and no known no, so it cannot be tested from FishBase at all.

## The tally

From `docs/analyses/mito-load/` (inputs, rules and columns in its README). Every rule is an explicit taxonomic rule
or a FishBase v26.06 field, except bat hibernation, read from the pinned table of Wilkinson & Adams 2019. Traits that
need a source we do not hold yet are not tallied: hibernation beyond bats, high altitude,
migration and diving in birds.

| Candidate (rule) | Scope | Yes | No | Unknown | Families | Clades | Changes | Yes with partners |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| Mammal flight: Chiroptera | 625 mammals | 44 | 581 | 0 | 8 | 1 | 1 | 13 |
| Mammal marine: cetaceans, Otariidae, Phocidae, Odobenidae, Sirenia, *Enhydra* | 625 | 50 | 575 | 0 | 16 | 4 | 4 | 26 |
| Mammal semi-aquatic: otters, *Castor*, platypus, hippos, *Neomys*, *Myocastor*, *Condylura* | 625 | 13 | 612 | 0 | 7 | 9 | 7 | 5 |
| Mammal subterranean: Bathyergidae, Spalacidae, Talpidae, *Spalacopus*, *Ellobius* | 625 | 7 | 618 | 0 | 5 | 5 | 5 | 2 |
| Mammal myrmecophagy: anteaters, echidnas, numbat, aardvark, aardwolf, two armadillos, sloth bear | 625 | 11 | 614 | 0 | 7 | 7 | 7 | 3 |
| Bird flightless: ratites, penguins, kakapo | 427 birds | 11 | 416 | 0 | 5 | 3 | 3 | 4 |
| Fish diadromous: FishBase `AnaCat` | 893 ray-finned | 122 | 368 | 403 | 36 | 69 | 59 | 22 |
| Fish deep: FishBase `DepthRangeDeep` ≥ 1000 m | 893 | 52 | 519 | 322 | 28 | 41 | 40 | 8 |
| Fish deep: FishBase `DemersPelag` bathydemersal or bathypelagic | 893 | 38 | 855 | 0 | 24 | 32 | 32 | 2 |
| Fish polar: FishBase `EnvTemp` "polar" in any stock | 893 | 42 | 851 | 0 | 15 | 34 | 30 | 3 |
| Fish air-breathing: FishBase `AirBreathing` | 893 | 32 | 0 | 861 | 21 | 1 | 0 | 3 |
| Fish regional endothermy: *Thunnus*, *Katsuwonus*, billfishes, *Lampris* | 893 | 16 | 877 | 0 | 4 | 4 | 4 | 2 |
| Shark regional endothermy: Lamnidae, Alopiidae | 87 cartilaginous | 7 | 80 | 0 | 2 | 3 | 3 | 0 |
| Electric discharge: FishBase `Electrogenic` (19 of the 23 are skates) | 980 fishes | 23 | 957 | 0 | 6 | 14 | 9 | 0 |
| Bat hibernation: Wilkinson & Adams 2019, within Chiroptera | 44 bats | 17 | 27 | 0 | 3 | 7 | 4 | 6 |

## The candidates, ranked

### Tier 1: many origins, classification in hand (FishBase v26.06, already pinned)

**Diadromy (shortlisted: H21).**
1. *Why.* Migration between fresh and salt water couples long sustained swimming with an osmoregulatory switch.
   Across 401 teleost mitogenomes, migratory fishes had lower Ka/Ks in the mitochondrial protein genes than
   non-migratory ones, read as stronger purifying selection from higher energy demand (Sun et al. 2011). Sedentary
   fishes showed relaxed purifying selection in all 13 genes, a pattern absent from nuclear housekeeping genes
   (Strohm et al. 2015).
2. *Classification.* FishBase v26.06 `species.AnaCat`: anadromous, catadromous, amphidromous or diadromous is yes;
   potamodromous, oceanodromous, non-migratory or oceano-estuarine is no; empty is unknown.
3. *Counts.* 122 yes, 368 no, 403 unknown; 22 yes in the partner set.
4. *Origins.* 59 changes, 69 trait clades, the most of any candidate. Diadromy is labile: in Clupeiformes alone,
   13 transitions to it and 7 losses were counted (Egan et al. 2024).
5. *Signal and tool.* Lower dN/dS on diadromous branches (RELAX, K > 1); a lower ω in pairs, each diadromous lineage
   against its nearest non-diadromous relative; on the coevolution side, a joint shift of mitochondrial and nuclear
   subunit rates on diadromous branches, if H19's fish set holds enough of them.

**The deep sea (shortlisted: H22).**
1. *Why.* Two opposite expectations, which is why this one is interesting. Oxygen consumption falls with depth in
   benthopelagic and, less, in benthic fishes (Drazen & Seibel 2007), which predicts relaxed constraint. Yet 11
   mitochondrial genes carry convergent or parallel amino-acid changes across nine independent deep-sea lineages,
   read as repeated adaptation (Shen et al. 2019).
2. *Classification.* FishBase `species.DemersPelag` bathydemersal or bathypelagic, against every other value (no
   unknowns among the matched fishes). Sensitivity: `DepthRangeDeep` ≥ 1000 m against < 200 m.
3. *Counts.* 38 yes and 855 no by habitat; 52 yes, 322 unknown by depth. Only 2 (habitat) or 8 (depth) with partners.
4. *Origins.* 32 changes by habitat, 40 by depth.
5. *Signal and tool.* RELAX and paired ω, two-sided; aBSREL on the stems of deep-sea clades. The coevolution side is
   not assessable with the current partner set.

**Polar fishes (not shortlisted).**
1. *Why.* Cold water slows enzyme kinetics, and cold-bodied fishes compensate with mitochondrial density. In Antarctic
   notothenioids, ND6 was first reported lost (Papetti et al. 2007) and then found relocated into a rearranged control
   region (Zhuang & Cheng 2010), a warning about annotation in this group.
2. *Classification.* FishBase `stocks.EnvTemp` "polar" in any stock.
3. *Counts.* 42 yes, 3 with partners.
4. *Origins.* 30 changes; the Antarctic notothenioids inside it are 2.
5. *Signal and tool.* As above. Left out because it overlaps H20's temperature, a continuous trait on the same
   fishes; it can be run as H20's categorical complement later.

**Air breathing (not testable from FishBase).** 32 yes, 19 to 21 families, but no confirmed no (see above).

### Tier 2: few origins, partners in hand (mammals)

**Hibernation (shortlisted: H23, and with longevity H24).**
1. *Why.* Hibernators suppress metabolism to about 6% of basal rate during torpor bouts, daily heterotherms to about
   35% (Ruf & Geiser 2015); Staples (2016) reviews how mitochondrial metabolism takes part in that suppression.
   Their mitochondria cycle between deep suppression and the high flux of arousal.
2. *Classification.* Ruf & Geiser's species table of 214 heterotherms was the first choice, but no script can fetch it
   (the publisher, PMC and the authors' repository each answer with a bot challenge). In its place, and better for a
   two-state trait because it codes non-hibernators too: the supplementary table of Heldstab, Isler & van Schaik 2018
   (1,104 mammals; Dryad, CC0), downloaded by hand once and accepted only with the md5 of Dryad's file listing; and for
   bats the trait table of Wilkinson & Adams 2019 (DRUM, SHA-256 pinned), `Hibernation` Yes or No. Both are fetched or
   checked by `download hibernation`. A species in neither is unknown; disagreement is unknown.
3. *Counts.* Bats: all 44 matched (42 by name, 2 by NCBI synonym), 17 hibernate, 27 do not, 6 hibernators with
   partners. Mammals beyond bats: not counted until the Dryad file is placed.
4. *Origins.* Bats: 4 changes, 7 trait clades (the hibernating Vespertilionidae with one tropical *Myotis* that does
   not hibernate, *Rhinolophus*, *Mystacina*). Torpor is widespread across mammals, so origins and losses will be
   hard to tell apart.
5. *Signal and tool.* Two-sided RELAX and paired ω; on the coevolution side, the joint mitochondrial and nuclear rate
   shift on hibernator branches.

**Marine mammals (shortlisted: H23).**
1. *Why.* Diving imposes hypoxia and reoxygenation. Across 41 placental mammals, cytochrome b varies widely at
   functionally important sites in "species with more-specialized metabolic requirements" (da Fonseca et al. 2008). Across the genomes of three
   independently marine mammals, convergent substitutions were common but just as common in terrestrial sister
   taxa, so adaptive convergence was rare (Foote et al. 2015). That is the control this design copies.
2. *Classification.* Cetacean families (whether filed under Cetacea or Artiodactyla, RM15 of the site repository),
   Otariidae, Phocidae, Odobenidae, Sirenia and *Enhydra*.
3. *Counts.* 50 species, 26 with partners.
4. *Origins.* 4 changes, 4 clades: a test bounds little.
5. *Signal and tool.* RELAX across the four clades, aBSREL on their four stems; the joint rate shift on the
   coevolution side.

**Hibernation and longevity (H24).** Hibernating mammals have about 15% higher annual survival than non-hibernators of
similar size, mature later and have longer generations, yet do not live longer than non-hibernators with similar
survival, so the link runs through lower mortality (Turbill, Bieber & Ruf 2011). Hibernating bats live about 6 years
longer on average than those that do not hibernate (Wilkinson & South 2002), and body mass and hibernation predict bat
longevity again in a larger sample (Wilkinson & Adams 2019). Ageing itself slows: epigenetic age essentially stalls
through the hibernation of yellow-bellied marmots (Pinho et al. 2022), and big brown bats age more slowly by a bat
epigenetic clock during hibernation (Sullivan et al. 2022). H11 found that mtDNA composition follows the slow life
history as a whole; hibernation is a slow life history reached many times. H24 asks whether part of lifespan's
composition signal is hibernation's, and whether hibernators accumulate mtDNA substitutions more slowly per year.

**Bats (the positive control of H23, and inside bats H25).**
1. *Why.* Both the mitochondrial and the nuclear OXPHOS genes show adaptive evolution on the ancestral branch of bats
   (Shen et al. 2010).
2. *Classification.* Order Chiroptera.
3. *Counts.* 44 species, 13 with partners.
4. *Origins.* One, so not a test of flight. It is the branch on which a method that claims to see a joint
   mitochondrial and nuclear shift must see it first.
5. *Signal and tool.* Nuclear OXPHOS subunits faster than cytosolic controls on the bat stem; aBSREL on the stem.
6. *Inside bats (H25).* Bats outlive non-flying mammals of their size, which Austad & Fischer (1991) read as reduced
   vulnerability to the environment rather than slow metabolism; as later reviews summarise it, heart mitochondria of
   the little brown bat produced half to a third of the reactive oxygen species of a shrew's and a mouse's
   (Brunet-Rossinni 2004). Within bats, flight is shared by every species and lifespan is far less tied to body mass,
   the confound of every mammal-wide test. 43 of the 44 bats have a lifespan record of high or acceptable quality.
   H25 tests lifespan against the substitution rate on sister pairs and against composition, and hibernation against
   constraint (RELAX across its 7 bat clades), with the power that 43 species and 4 origins allow.

### Tier 3: several origins, few species (mitochondrial genes only, not shortlisted)

- **Subterranean mammals** (7 species, 5 changes). dN/dS was higher in subterranean tuco-tucos and the coruro for 11 of
  13 genes, none above 1 (Tomasco & Lessa 2011).
- **Myrmecophagy** (11 species, 7 changes). Ant- and termite-eaters converge on low metabolism (McNab 1984).
- **Semi-aquatic mammals** (13 species, 7 changes).

These suit the paired design on the mitochondrial genes, with 5 to 7 pairs each: below H21-H23's minimum of 8.

### Tier 4: too few origins or no source (described, not tested)

- **Regional endothermy** (16 teleosts, 4 changes; 7 sharks, 3 changes). Dickson & Graham (2004) review where endothermy
  occurs among fishes and test its proposed advantages on the phylogeny. No elevated dN/dS was found in COX2 on the tuna lineages (Dalziel et al. 2006).
- **Flightless birds** (11 species, 3 changes: the build's three ratites have no flying palaeognath among them). Weak
  fliers and flightless birds accumulated relatively more nonsynonymous substitutions (Shen et al. 2009); a later test
  of positive against relaxed selection across animal comparisons, flightless birds among them, found that elevated
  dN/dS was often not caused by the mechanism predicted (Zwonitzer et al.).
- **Electric discharge** (23 species, 19 of them skates).
- **High altitude** in mammals and birds: strong biology (a COX3 substitution in the bar-headed goose changes the
  enzyme's kinetics; Scott et al. 2011), but no pinned elevation source. The build has *Bos grunniens*, *Lama*,
  *Vicugna*, *Ochotona*, *Panthera uncia*, *Montifringilla*, *Pyrrhocorax*; FishBase marks 4 fishes "high altitude".
- **Migration and diving in birds**: need AVONET (migration, Tobias et al. 2022) or EltonTraits (foraging below the
  water surface, Wilman et al. 2014; its field names not yet checked). Adding either is the maintainer's call.
- **Insects, parasites**: fewer than five species each in the build.

## The shortlist and its design

Agreed with the maintainer on 2026-11-16 (the answer to "Which shortlist should the document and hypotheses be
written for?" was "Fish ecology + mammals"):
- **H21**, diadromous fishes;
- **H22**, deep-sea fishes;
- **H23**, hibernating and marine mammals, with bats as the positive control;
- and, after the maintainer's follow-up the same day: **H24**, hibernation and longevity across mammals; **H25**, lifespan,
  hibernation and the mitochondrial genome inside bats.

They share one design, defined in full in `docs/HYPOTHESES.md`:
- *Mitochondrial side, one contrast per origin:* each trait clade paired with its nearest relative outside it, and an
  outgroup, so that a lineage that radiated counts once. RELAX with the trait branches as test branches, against
  control draws of matched non-trait clades (RELAX on 12 linked genes needs its own calibration).
- *Coevolution side:* the branch lengths of every protein on one fixed topology (H5's `branch-lengths`, IQ-TREE). On
  trait branches, do the nuclear subunits of the respiratory complexes speed up or slow down beyond the cytosolic
  controls, with the mitochondrial genes? Bats first: if the method does not see their stem, nothing else is read.
- *Composition:* H11's 43 features against the trait, exploratory.

## Sources to add, the maintainer's call

- **Hibernation** (H23-H25): approved by the maintainer on 2026-11-16 ("Yes, write both", to writing H24 and H25 and
  pinning both hibernation tables). The bat table of Wilkinson & Adams 2019 is pinned and fetched. The mammal table of
  Heldstab et al. 2018 waits for one download by hand: "Supplementary material.xlsx" from
  https://datadryad.org/dataset/doi:10.5061/dryad.753d06g, saved as
  `data/input/hibernation/heldstab_2018_supplementary_material.xlsx`. Ruf & Geiser 2015 could not be fetched (above).
- AVONET (bird migration, hand-wing index); EltonTraits 1.0 (diving in birds and mammals); COMBINE (mammal traits;
  whether it carries a hibernation field was not confirmed). None is needed by H21-H23.

## References

- Austad SN, Fischer KE (1991) Mammalian aging, metabolism, and ecology: evidence from the bats and marsupials.
  J Gerontol 46:B47-B53. doi:10.1093/geronj/46.2.b47
- Brunet-Rossinni (2004) Reduced free-radical production and extreme longevity in the little brown bat (Myotis
  lucifugus) versus two non-flying mammals. Mech Ageing Dev 125. doi:10.1016/j.mad.2003.09.003
- Dalziel AC, Moyes CD, Fredriksson E, Lougheed SC (2006) Molecular evolution of cytochrome c oxidase in
  high-performance fish (Teleostei: Scombroidei). J Mol Evol 62:319-331. doi:10.1007/s00239-005-0110-7
- da Fonseca RR, Johnson WE, O'Brien SJ, Ramos MJ, Antunes A (2008) The adaptive evolution of the mammalian
  mitochondrial genome. BMC Genomics 9:119. doi:10.1186/1471-2164-9-119
- Dickson KA, Graham JB (2004) Evolution and consequences of endothermy in fishes. Physiol Biochem Zool 77:998-1018.
  doi:10.1086/423743
- Drazen, Seibel (2007) Depth-related trends in metabolism of benthic and benthopelagic deep-sea fishes. Limnol
  Oceanogr 52:2306-2316. doi:10.4319/lo.2007.52.5.2306
- Egan et al. (2024) Phylogenomics, lineage diversification rates, and the evolution of diadromy in Clupeiformes
  (anchovies, herrings, sardines, and relatives). Syst Biol. PMID 38756097
- Foote AD et al. (2015) Convergent evolution of the genomes of marine mammals. Nat Genet 47:272-275.
  doi:10.1038/ng.3198
- Heldstab SA, Isler K, van Schaik CP (2018) Hibernation constrains brain size evolution in mammals. J Evol Biol.
  doi:10.1111/jeb.13353. Data: Dryad doi:10.5061/dryad.753d06g
- McNab BK (1984) Physiological convergence amongst ant-eating and termite-eating mammals. J Zool 203:485-510.
  doi:10.1111/j.1469-7998.1984.tb02345.x
- Papetti et al. (2007) Antarctic fish mitochondrial genomes lack ND6 gene. J Mol Evol. doi:10.1007/s00239-007-9030-z
- Pinho et al. (2022) Hibernation slows epigenetic ageing in yellow-bellied marmots. Nat Ecol Evol 6:418-426.
  doi:10.1038/s41559-022-01679-1
- Ruf, Geiser (2015) Daily torpor and hibernation in birds and mammals. Biol Rev 90:891-926. doi:10.1111/brv.12137
- Scott GR et al. (2011) Molecular evolution of cytochrome c oxidase underlies high-altitude adaptation in the
  bar-headed goose. Mol Biol Evol 28:351-363. doi:10.1093/molbev/msq205
- Shen et al. (2009) Relaxation of selective constraints on avian mitochondrial DNA following the degeneration of
  flight ability. Genome Res 19:1760.
- Shen et al. (2010) Adaptive evolution of energy metabolism genes and the origin of flight in bats. PNAS
  107:8666-8671. PMC2889356
- Shen, Pu, Chen, Murphy, Shen (2019) Convergent evolution of mitochondrial genes in deep-sea fishes. Front
  Genet 10:925. doi:10.3389/fgene.2019.00925
- Staples JF (2016) Metabolic flexibility: hibernation, torpor, and estivation. Compr Physiol 6:737-771.
- Strohm, Gwiazdowski, Hanner (2015) Fast fish face fewer mitochondrial mutations: patterns of dN/dS across fish
  mitogenomes. Gene 572:27-34. doi:10.1016/j.gene.2015.06.074
- Sullivan IR, Adams DM, Greville LJS, Faure PA, Wilkinson GS (2022) Big brown bats experience slower epigenetic ageing
  during hibernation. Proc R Soc B 289. doi:10.1098/rspb.2022.0635
- Sun, Shen, Irwin, Zhang (2011) Evaluating the roles of energetic functional constraints on teleost
  mitochondrial-encoded protein evolution. Mol Biol Evol 28:39-44. doi:10.1093/molbev/msq256
- Tobias et al. (2022) AVONET: morphological, ecological and geographical data for all birds. Ecol Lett 25:581-597.
  doi:10.1111/ele.13898
- Tomasco IH, Lessa EP (2011) The evolution of mitochondrial genomes in subterranean caviomorph rodents: adaptation
  against a background of purifying selection. Mol Phylogenet Evol 61:64-70. doi:10.1016/j.ympev.2011.06.014
- Turbill C, Bieber C, Ruf T (2011) Hibernation is associated with increased survival and the evolution of slow life
  histories among mammals. Proc R Soc B 278:3355-3363. doi:10.1098/rspb.2011.0190
- Wilkinson GS, Adams DM (2019) Recurrent evolution of extreme longevity in bats. Biol Lett 15:20180860.
  doi:10.1098/rsbl.2018.0860. Data: DRUM hdl:1903/21501
- Wilkinson GS, South JM (2002) Life history, ecology and longevity in bats. Aging Cell 1:124-131.
  doi:10.1046/j.1474-9728.2002.00020.x
- Wilman H et al. (2014) EltonTraits 1.0: species-level foraging attributes of the world's birds and mammals. Ecology
  95:2027. doi:10.1890/13-1917.1
- Zhuang, Cheng (2010) ND6 gene "lost" and found: evolution of mitochondrial gene rearrangement in Antarctic
  notothenioids. Mol Biol Evol 27:1391-1403. doi:10.1093/molbev/msq026
- Zwonitzer et al. Disentangling positive vs. relaxed selection in animal mitochondrial genomes. bioRxiv
  doi:10.1101/2022.10.05.510972; published version PMC10955554.
