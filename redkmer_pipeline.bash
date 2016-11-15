#!/bin/bash

runfile="$1"
if source ${runfile}; then
printf "======= redkmer 1.0 =======\n"
printf "Obtained run data from ${runfile}\n"
printf "Working Directory: ${CWD}\n"
printf "Pacbio Read Directory: ${pacDIR}\n"
printf "Illumina Read Directory: ${illDIR}\n"
printf "Running full pipeline !\n"

else
printf 'Failed to obtain run data. Exiting!\n'
exit 0
fi

source redkmer.cfg

bash ${BASEDIR}/1_redkmer_QCandMitocleanup.bash ${runfile}
bash ${BASEDIR}/2_redkmer_pacbins.bash ${runfile}
R CMD BATCH --no-save --no-restore ${BASEDIR}/2R_analysis_pacBioBins.R
bash ${BASEDIR}/3_redkmer_kmerGenerate.bash ${runfile}
bash ${BASEDIR}/4_redkmer_kmersbowtie.bash ${runfile}
bash ${BASEDIR}/5_redkmer_uniqueness.bash ${runfile}
R CMD BATCH --no-save --no-restore ${BASEDIR}/5R_analysis_kmers.R
bash ${BASEDIR}/6_redkmer_kmer2genome.bash ${runfile}
R CMD BATCH --no-save --no-restore ${BASEDIR}/6R_analysis_kmers2genome.R


mv ${BASEDIR}/*.Rout ${CWD}/plots
