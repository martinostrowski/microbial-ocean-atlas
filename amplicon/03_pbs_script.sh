#!/bin/bash
 
#PBS -J 41-50
#PBS -N p41-50
#PBS -l ncpus=10
#PBS -l mem=100GB
#PBS -l walltime=100:00:00
#PBS -q workq


#module load devel/R-current;
module load devel/c8/R-4.0.2


cd /shared/c3/bio_db/BPA/amplicons/16s
echo "Job ID is ${PBS_JOBID}"
echo "Job Array ID is ${PBS_ARRAY_INDEX}"
echo "Timestamp is $(date +%F_%T)"
echo "Directory is $(pwd)"
echo "Running on host $(hostname)"
echo "Working directory is ${PBS_O_WORKDIR}"
echo "Job has the following nodes/cores:"
cat ${PBS_NODEFILE}


#PARAMETERS=$(awk -v line=${PBS_ARRAY_INDEX} '{if (NR == line) { print $0; };}' the.conf)

date +%F_%T


echo "$PBS_ARRAY_INDEX"


Rscript --verbose dada2.pipeline.r $PBS_ARRAY_INDEX > dd.${PBS_ARRAY_INDEX}.dada2.pipeline.r.out;
