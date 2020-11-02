This is a custom script to process Illumina MiSeq amplicon data directly modified from the workflows published on the official 
[DADA2 website](https://benjjneb.github.io/dada2/index.html).
The script includes some modifications to fulfill specific needs and has proven to work well for several plates of Illumina sequencing run in succession throught the same pipeline before being compiled, collapsed using collapse_no_mismatch, and finally being exported as an amplicon sequence variants (ASV) table.

The code used here was adapted from:

Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP. DADA2: high-resolution sample inference from Illumina amplicon data. Nat Methods 2006;13:581–3.

Callahan, Benjamin J, Paul J McMurdie, and Susan P Holmes. 2017. “Exact Sequence Variants Should Replace Operational Taxonomic Units in Marker Gene Data Analysis.” bioRxiv. Cold Spring Harbor Labs Journals, 113597.

# Marine-Microbes-dada2-pipeline
This is the R code used for processing the paired end Illumina reads for the all of the Marine Microbes pelagic project samples as of August 2020. Raw paired end reads were obtained from the Bioplatforms Australia data portal (Australian Microbiome data portal: https://data.bioplatforms.com/organization/about/australian-microbiome) and were run through a customised version of dada2 pipeline in order to check the read quality, remove forward and reverse primers from the reads using cutadapt, and truncate the reads to eliminate low quality terminal bases. Error rates were learned based on 1e8 bases and max consist =20. Reads were then dereplicated and merged using pseudo pooling. The chimeras were then removed using `minFoldParentOverAbundance`=1. Finally all plates run through the pipeline were saved as individual seqtab files, they were merged and together run through the `collapseNoMismatch` step in order to collapse two identical merged reads without internal mismatches, but having different terminal lengths into one ASV, combining their abundances and retaining the most abundant of the two collapsed reads. Reads were then assigned using different databases depending on their target organism, this pipeline was used for three Marine Microbes time series, each with their own pipeline, primers and truncLengths: 

| Target | F primer | R primer |
| :------------- |:------------- |:-----|
Archaea a16s rRNA | (A2F/Arch21f): 5’-TTCCGGTTGATCCYGCCGGA-3’ | (519 R*): 5’-GWATTACCGCGGCKGCTG-3’ |
Bacterial 16s rRNA | (27 F): 5’-AGAGTTTGATCMTGGCTCAG-3’| (519 R*): 5’-GWATTACCGCGGCKGCTG-3’ |
Eukaryotic 18s rRNA |  (TAReuk454FWD1): 5’-CCAGCASCYGCGGTAATTCC-3’ | (TAReuk-Rev3): 5’-ACTTTCGTTCTTGATYRATGATCTRYATC-3’|

***

## Summary of processes

01. Preparation
02. Main DaDa2 pipeline
  * [Archaea](../Archaea_a16s_rRNA)
  * [Bacteria](../Bacteria_b16s_rRNA)
  * [Eukaryotes](../Euk_e18s_rRNA)

03. Processing multiple sequencing runs on a HPC using PBS
04. Building a unified ASV table and running collapse_no_mismatch
05. Assigning taxonomy:
  * GTDB (for a16s and/or b16s)
  * Silva_v138 (for b16s)
  * PR2 (for e18s)
06. Pre-analysis options



