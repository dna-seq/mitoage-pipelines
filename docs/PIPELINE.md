# The pipeline in pictures

What is downloaded, which command reads and writes what, and which table each finding is read from.
Hypotheses (H0 to H6), their tests and results are in [HYPOTHESES.md](HYPOTHESES.md).

## Where the data comes from

```mermaid
flowchart LR
  subgraph sources [Original sources]
    HF[(MitoAge build<br/>Hugging Face)]
    SITE[(dated species tree<br/>site repository)]
    NCBI[(NCBI ortholog sets)]
  end
  subgraph prepared [Prepared copies on Hugging Face]
    P1[(nuclear partners)]
    P2[(embeddings)]
    P3[(aligner variants)]
    P4[(results)]
  end
  HF -- "download build" --> IN[data/input]
  SITE -- "download tree" --> IN
  NCBI -- "download orthologs" --> IN
  P1 -- "download prepared" --> IN
  P2 -- "download prepared" --> IN
  P3 -- "download prepared" --> INT[data/interim]
  P1 -- "download prepared" --> INT
  P4 -- "download prepared" --> REF[data/reference]
  IN --> STEPS{{mitoage commands}}
  INT --> STEPS
  STEPS --> INT
  STEPS --> OUT[data/output]
  OUT -. compare .- REF
```

Two ways in: fetch from the original sources and run every step, or fetch our prepared copies and start
anywhere. `data/reference` holds our result tables; a rerun writes `data/output`, and the two can be compared
file by file.

## Longevity in the mitochondrial proteins (H0 to H3)

```mermaid
flowchart TD
  B[build: species.parquet<br/>13 protein and codon alignments] --> C[longevity controls]
  B --> F[longevity fourfold]
  C --> COMP[composition.tsv]
  C --> T1[controls.tsv<br/>H1: amino acids vs lifespan]
  F --> FF[fourfold.tsv]
  COMP --> C2[longevity composition]
  FF --> C2
  C2 --> T2[controls_composition.tsv<br/>H1 with base composition]
  PJ[phylogeny.json] --> DT[longevity dated-tree]
  B --> DT
  DT --> TREE[dated.nwk<br/>quality_ok.txt]
  COMP --> S[longevity sites]
  FF --> S
  S --> T3[sites.tsv<br/>H2 without a tree: inflated]
  COMP --> PG[longevity pgls]
  FF --> PG
  TREE --> PG
  PG --> T4[pgls.Mammalia.tsv, pgls.Aves.tsv<br/>H2 as fixed]
  TREE --> RP[longevity relax-prep]
  COMP --> RP
  RP --> RX[longevity relax<br/>HyPhy]
  RX --> T5[relax/GENE.relax.json<br/>H3]
  AL[aligners/GENE.tool.afa] --> CA[longevity compare-aligners]
  CA --> T0[agreement.tsv<br/>H0]
```

Every model that tests lifespan carries body mass and base composition at fourfold-degenerate sites; every
step that needs a tree takes the dated species tree, and only lifespan records of high or acceptable quality.

## Peptides encoded inside mtDNA (H4)

```mermaid
flowchart LR
  B[build: rRNA and codon alignments<br/>species.parquet] --> PP[peptides pilot]
  PP --> T[humanin.tsv, SHLP6.tsv, MOTS-c.tsv, MTALTND4.tsv<br/>summary.txt]
```

For each published peptide: find it in the human row, project its columns onto every species, ask for a
start codon and an open frame; compare with other windows of the same gene that are open in human.

## Mitochondrial and nuclear partners: mitoport (H5, H6)

```mermaid
flowchart TD
  N[ncbi/GENE: proteins, coding sequences] --> SEL[mitoport select]
  SEL --> S1[sel: one protein per species and gene]
  S1 --> AL[mitoport align<br/>MAFFT]
  B[build: 13 protein alignments] --> AL
  TREE[dated.nwk] --> AL
  AL --> A1[aln, genes.txt, species.txt, topology.nwk]
  A1 --> BL[mitoport branch-lengths<br/>IQ-TREE, fixed topology]
  BL --> B1[bl/GENE.treefile]
  B1 --> ERC[mitoport erc]
  ERC --> T1[erc.tsv<br/>H5: rates of partners over branches]
  E[embeddings: ESM-2 mean per protein] --> CO[mitoport coembed]
  TREE --> CO
  CO --> T2[coembed.tsv<br/>H6: embedding distances of partners]
```

Both tests compare nuclear subunits of complexes III and IV with three kinds of control: cytosolic
proteins, complex II (which has no mitochondrially encoded subunit) and mitochondrial matrix enzymes.
