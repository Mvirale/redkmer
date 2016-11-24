#!/bin/bash
#PBS -N redkmer_step1
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=20:mem=16gb
#PBS -e /home/nikiwind/reports
#PBS -o /home/nikiwind/reports

if [[ "$RUNINCLUSTER" -eq "1" ]]; then
source $PBS_O_WORKDIR/redkmer.cfg
module load fastqc
module load bowtie/1.1.1
else
source redkmer.cfg
fi

echo "========== setting up =========="

mkdir -p $CWD/QualityReports
mkdir -p $CWD/pacBio_illmapping
mkdir -p $CWD/pacBio_illmapping/logs
mkdir -p $CWD/pacBio_illmapping/mapping_rawdata
mkdir -p $CWD/pacBio_illmapping/index
mkdir -p $CWD/pacBio_bins
mkdir -p $CWD/pacBio_bins/fasta
mkdir -p $CWD/temp
mkdir -p $CWD/kmers
mkdir -p $CWD/kmers/rawdata
mkdir -p $CWD/kmers/fasta/
mkdir -p $CWD/plots
mkdir -p $CWD/kmers/Refgenome_blast
mkdir -p $CWD/kmers/bowtie
mkdir -p $CWD/kmers/bowtie/index
mkdir -p $CWD/kmers/bowtie/mapping
mkdir -p $CWD/kmers/bowtie/mapping/logs
mkdir -p $CWD/kmers/Refgenome_blast
mkdir -p $CWD/kmers/bowtie/offtargets
mkdir -p $CWD/kmers/bowtie/offtargets/logs

echo "========== producing quality report for illumina libraries =========="

$FASTQC ${CWD}/${illDIR}/raw_f.fastq -o $CWD/QualityReports/
$FASTQC ${CWD}/${illDIR}/raw_m.fastq -o $CWD/QualityReports/

echo "========== removing illumina reads mapping to mitochondrial DNA =========="

mkdir -p $CWD/MitoIndex

# Build the index and map the Illumina data
$BOWTIEB $MtREF $CWD/MitoIndex/MtRef
#$BOWITE2B $MtREF $CWD/MitoIndex/MtRef_bowtie2

# Map the Illumina data on the mito, the option  --un gives the unmapped read (not mitochondrial)
$BOWTIE -p $CORES $CWD/MitoIndex/MtRef ${CWD}/${illDIR}/raw_f.fastq --un ${CWD}/${illDIR}/f.fastq 2> ${CWD}/${illDIR}/f_bowtie.log
$BOWTIE -p $CORES $CWD/MitoIndex/MtRef ${CWD}/${illDIR}/raw_m.fastq --un ${CWD}/${illDIR}/m.fastq 2> ${CWD}/${illDIR}/m_bowtie.log
#($BOWTIE2 -p $CORES -x $CWD/MitoIndex/MtRef_bowtie2 ${CWD}/${illDIR}/raw_f.fastq --un ${CWD}/${illDIR}/f_bowtie2.fastq) 2> ${CWD}/${illDIR}/f_bowtie2.log

printf "======= Done step 1 =======\n"

