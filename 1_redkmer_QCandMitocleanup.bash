#!/bin/bash

runfile="$1"
if source ${runfile}; then
printf "======= redkmer 0.1 =======\n"
printf "Obtained run data from ${runfile}\n"
printf "Working Directory: ${CWD}\n"
printf "Pacbio Read Directory: ${pacDIR}\n"
printf "Illumina Read Directory: ${illDIR}\n"
printf "Running script.\n"

else
printf 'Failed to obtain run data. Exiting!\n'
exit 0
fi

source redkmer.cfg

mkdir -p $CWD/QualityReports

echo "========== producing quality report for illumina libraries =========="

$FASTQC ${CWD}/${illDIR}/raw_f.fastq -o $CWD/QualityReports/
$FASTQC ${CWD}/${illDIR}/raw_m.fastq -o $CWD/QualityReports/

echo "========== removing illumina reads mapping to mitochondrial DNA =========="

mkdir -p $CWD/MitoIndex

# Build the index and map the Illumina data
$BOWTIEB $MtREF $CWD/MitoIndex/MtRef

# Map the Illumina data on the mito, the option  --un gives the unmapped read (not mitochondrial)
$BOWTIE -a -p$CORES -v 2 $CWD/MitoIndex/MtRef ${CWD}/${illDIR}/raw_f.fastq --un ${CWD}/${illDIR}/f.fastq --al ${CWD}/${illDIR}/f_mito.fastq
$BOWTIE -a -p$CORES -v 2 $CWD/MitoIndex/MtRef ${CWD}/${illDIR}/raw_m.fastq --un ${CWD}/${illDIR}/m.fastq --al ${CWD}/${illDIR}/m_mito.fastq


printf "======= Done step 1 =======\n"
