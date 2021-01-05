#!/bin/bash
#PBS -N redkmer5R
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=16:mem=16gb

module purge
module load R

export R_LIBS="/home/nikiwind/localRlibs"

cd $EPHEMERAL/redkmer/Rscripts/

R CMD BATCH --no-save --no-restore 5R_analysis_kmers.R



