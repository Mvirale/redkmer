#!/bin/bash
#PBS -N redkmer2R
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=16:mem=16gb

module purge
module load R

cd $EPHEMERAL/redkmer/Rscripts

R CMD BATCH --no-save --no-restore 2R_analysis_pacBioBins.R


