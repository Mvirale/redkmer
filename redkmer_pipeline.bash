#!/bin/bash
#PBS -N redkmer_pipeline
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=8:mem=40gb

bash 1_redkmer_QCandMitocleanup.bash
bash 2_redkmer_pacbins.bash

cd $EPHEMERAL/redkmer/Rscripts
R CMD BATCH --no-save --no-restore 2R_analysis_pacBioBins.R
cd ..

bash 3_redkmer_kmerGenerate.bash
bash 4_redkmer_kmersbowtie.bash
bash 5_redkmer_uniqueness.bash

cd $EPHEMERAL/redkmer/Rscripts
R CMD BATCH --no-save --no-restore 5R_analysis_kmers.R
cd ..

bash 6_redkmer_kmer2genome.bash

cd $EPHEMERAL/redkmer/Rscripts
R CMD BATCH --no-save --no-restore 6R_analysis_kmers2genome.R

