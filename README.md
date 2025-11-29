
# Positive selection analsyis for Nectar genomics project

associated with preprint [Convergent and lineage-specific genomic changes shape adaptations in sugar-consuming birds](https://www.biorxiv.org/content/10.1101/2024.08.30.610474v2)


This pipeline describes the following steps:

## 0. Getting orthologs from genome alignments with [TOGA](https://github.com/hillerlab/TOGA)

## 1. [Extracting alignments of one2one orthologs](https://github.com/osipovarev/absrel/blob/main/README_one2ones.md)

## 2. [Running aBSREL (HyPhy) screen for selection](https://github.com/osipovarev/absrel/blob/main/README_run_absrel.md)

for details on the method see: [HyPhy](https://stevenweaver.github.io/hyphy-site/methods/selection-methods/)


## 3. [Analysis of aBSREL results](https://github.com/osipovarev/absrel/blob/main/absrel_analysis.ipynb)

## 4. [Functional enrichment analysis with ClusterProfiler](https://github.com/osipovarev/absrel/blob/main/absrel_analysis_2024/README_enrich_analsyis.md)

### Scripts
assocaited scripts can be found in the following repos: \
- [selection_screen_tools](https://github.com/osipovarev/selection_screen_tools) 
- [Enrichment_analysis](https://github.com/osipovarev/Enrichment_analysis) 
- [Annotation_scripts](https://github.com/osipovarev/Annotation_scripts)
- [msa_tools](https://github.com/osipovarev/msa_tools/)
- [fasta_tools](https://github.com/osipovarev/fasta_tools/)

***

### Additional

#### [How to make gene MSA for Figures](https://github.com/osipovarev/absrel/blob/main/README_msa_for_figures.md)