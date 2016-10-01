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

bash ${CWD}/1_redkmer_qc.bash ${runfile}

bash ${CWD}/2_redkmer_pacbins.bash ${runfile}

bash ${CWD}/3_redkmer_kmers.bash ${runfile}

bash ${CWD}/4_redkmer_kmersblast.bash ${runfile}

