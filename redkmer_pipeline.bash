#!/bin/bash

runfile="$1"
if source ${runfile}; then
printf "======= redkmer 0.2 =======\n"
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
bash ${BASEDIR}/3_redkmer_kmerGenerate.bash ${runfile}
bash ${BASEDIR}/4_redkmer_kmersblast.bash ${runfile}
bash ${BASEDIR}/5_redkmer_uniqueness.bash ${runfile}
bash ${BASEDIR}/6_redkmer_kmer2genome.bash ${runfile}
