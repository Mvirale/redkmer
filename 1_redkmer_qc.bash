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

mkdir -p $CWD/qualityreports

echo "========== quality checking illumina reads and producing report =========="

$FASTQC $illF -o $CWD/qualityreports/
$FASTQC $illM -o $CWD/qualityreports/

# echo "========== filtering for read quality =========="
# 
# ##implement QC for reads
# 
# mv ${illDIR}/ultraraw_f.fastq ${illDIR}/raw_f.fastq
# mv ${illDIR}/ultraraw_m.fastq ${illDIR}/raw_m.fastq
# 
# 
# echo "========== filter for mitochondrial DNA =========="
# 
# mkdir -p $CWD/index
# 
# # Build the index and map the Illumina data
# $BOWTIEB $mtREF $CWD/index/MtRef
# 
# # Map the Illumina data on the mito, the option  --un gives the unmapped read (not mitochondrial)
# $BOWTIE -a -p20 -v 0 $CWD/index/MtRef ${illDIR}/raw_f.fastq --un ${illDIR}/f.fastq --al ${illDIR}/f_mito.fastq
# $BOWTIE -a -p20 -v 0 $CWD/index/MtRef ${illDIR}/raw_m.fastq --un ${illDIR}/m.fastq --al ${illDIR}/m_mito.fastq


