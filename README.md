This is a custom script to process Illumina MiSeq amplicon data directly modified from the workflows published on the official 
[DADA2 website](https://benjjneb.github.io/dada2/index.html).
The script includes some modifications to fulfill specific needs and has proven to work well for several plates of Illumina sequencing run in succession throught the same pipeline before being compiled, collapsed using collapse_no_mismatch, and finally being exported as an amplicon sequence variants (ASV) table.

The code used here was adapted from:

Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP. DADA2: high-resolution sample inference from Illumina amplicon data. Nat Methods 2006;13:581–3.

Callahan, Benjamin J, Paul J McMurdie, and Susan P Holmes. 2017. “Exact Sequence Variants Should Replace Operational Taxonomic Units in Marker Gene Data Analysis.” bioRxiv. Cold Spring Harbor Labs Journals, 113597.

### Tools

R
cutadapt

R packages (dada2, tidyverse, ...)

### Inputs

directories containing PE fastqs and contols

### Outputs

ASV abundance tables for each flow_id in two formats: *seqtab format* xxx.RDS , *long format* xxx.csv for Chim1 and Chim4, i.e. 2 per sequencing run

A unified ASV table that has been fully or partially collapsed in *seqtab format* xxx.RDS and *long format* xxx.csv

# Marine-Microbes-dada2-pipeline
This is the R code used for processing the paired end Illumina reads for the all of the Marine Microbes pelagic project samples as of August 2020. Raw paired end reads were obtained from the Bioplatforms Australia data portal (Australian Microbiome data portal: https://data.bioplatforms.com/organization/about/australian-microbiome) and were run through a customised version of dada2 pipeline in order to check the read quality, remove forward and reverse primers from the reads using cutadapt, and truncate the reads to eliminate low quality terminal bases. Error rates were learned based on 1e8 bases and max consist =20. Reads were then dereplicated and merged using pseudo pooling. The chimeras were then removed using `minFoldParentOverAbundance`=1. Finally all plates run through the pipeline were saved as individual seqtab files, they were merged and together run through the `collapseNoMismatch` step in order to collapse two identical merged reads without internal mismatches, but having different terminal lengths into one ASV, combining their abundances and retaining the most abundant of the two collapsed reads. Reads were then assigned using different databases depending on their target organism, this pipeline was used for three Marine Microbes time series, each with their own pipeline, primers and truncLengths: 

| Target | F primer | R primer |
| :------------- |:------------- |:-----|
Archaea a16s rRNA | (A2F/Arch21f): 5’-TTCCGGTTGATCCYGCCGGA-3’ | (519 R*): 5’-GWATTACCGCGGCKGCTG-3’ |
Bacterial 16s rRNA | (27 F): 5’-AGAGTTTGATCMTGGCTCAG-3’| (519 R*): 5’-GWATTACCGCGGCKGCTG-3’ |
Eukaryotic 18s rRNA |  (TAReuk454FWD1): 5’-CCAGCASCYGCGGTAATTCC-3’ | (TAReuk-Rev3): 5’-ACTTTCGTTCTTGATYRATGATCTRYATC-3’|

***

## Summary of processes

01. [Preparation](../01_preparation.md)
02. [Main DaDa2 pipeline](../02_main_pipeline.md)
  * [Archaea](../02_Archaea_a16s)
  * [Bacteria](../02_Bacteria_b16s)
  * [Eukaryotes](../02_Euk_e18s)

03. Processing multiple sequencing runs on a HPC using PBS. [example PBS script](../03_pbs_script)
04. [Building a unified ASV table](../04_build_table.md) *including `collapseNoMismatch`*
05. Taxonomic Assignment
  * [GTDB (for a16s and/or b16s)](../05_Archaea_a16s)
  * [Silva v138 (b16s)](../05_Bacteria_b16s)
  * [PR2 (for e18s and chloroplasts)](../05_Euk_e18s)
06. Pre-analysis options



