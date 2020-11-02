## Preparation steps 

1. Download PE fastq files from the bioplatforms data portal and store them in directory corresponding to the flow_id of the sequencing run for each locus (i.e. Bacteria, Archaea, Eukaryotes)
2. build a run list (plates.csv) which consists of the flow_ids with the follwing structure

```r
plate,run,trunL,trunR,yeild
C3VFP,1,255,250,NA
C5PBC,2,255,250,NA
C7NP7,3,255,250,NA
B9YMG,4,255,250,NA
B6FJD,5,255,250,NA
B9YMG,6,255,250,NA
BH3PP,7,255,250,NA
```

The run list can be used to generate a [PBS array job](../03_pbs_script) to process multiple sequencing runs on a HPC



