#!/bin/bash
#PBS -N redkmer_R
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=20:mem=16gb
#PBS -e /home/nikiwind/reports
#PBS -o /home/nikiwind/reports

module purge
module load R

cd $PBS_O_WORKDIR/

#R CMD BATCH --no-save --no-restore 2R_analysis_pacBioBins.R
R CMD BATCH --no-save --no-restore 5R_analysis_kmers.R
#R CMD BATCH --no-save --no-restore 6R_analysis_kmers2genome.R
#Rscript Rtest.r




