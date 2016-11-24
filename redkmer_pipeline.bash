#!/bin/bash


bash 1_redkmer_QCandMitocleanup.bash
bash 2_redkmer_pacbins.bash
R CMD BATCH --no-save --no-restore 2R_analysis_pacBioBins.R
bash 3_redkmer_kmerGenerate.bash
bash 4_redkmer_kmersbowtie.bash
bash 5_redkmer_uniqueness.bash
R CMD BATCH --no-save --no-restore 5R_analysis_kmers.R
bash 6_redkmer_kmer2genome.bash
R CMD BATCH --no-save --no-restore 6R_analysis_kmers2genome.R


rm *Rout
