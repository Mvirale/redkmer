#!/bin/bash
#PBS -N redkmer1
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=24:mem=16gb



if [ -z ${PBS_ENVIRONMENT+x} ]
then
echo "---> running on the Perugia numbercruncher..."
source redkmer.cfg
else
echo "---> running on HPC cluster..."
source $PBS_O_WORKDIR/redkmer.cfg
module load bowtie/1.1.1
fi


echo "========== setting up =========="

mkdir -p $CWD/qsubscripts
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
mkdir -p $CWD/MitoIndex
mkdir -p $CWD/reports


if [ -z ${PBS_ENVIRONMENT+x} ]
then
	# Build the index and map the Illumina data
	$BOWTIEB $MtREF ${CWD}/MitoIndex/MtRef
	# Map the Illumina data on the mito, the option  --un gives the unmapped read (not mitochondrial)
	$BOWTIE -p $CORES $CWD/MitoIndex/MtRef ${illDIR}/raw_f.fastq --un ${illDIR}/f.fastq 2> ${illDIR}/f_bowtie.log
	$BOWTIE -p $CORES $CWD/MitoIndex/MtRef ${illDIR}/raw_m.fastq --un ${illDIR}/m.fastq 2> ${illDIR}/m_bowtie.log
else

echo "========== generating mitochondrial index =========="

$BOWTIEB $MtREF ${CWD}/MitoIndex/MtRef

cat > ${CWD}/qsubscripts/femalemito.bashX <<EOF
#!/bin/bash
#PBS -N redkmer_f_mito
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=24:mem=16gb
#PBS -e $CWD/reports
#PBS -o $CWD/reports

module load bowtie/1.1.1
module load fastqc
module load anaconda3/personal

cp ${illDIR}/raw_f.fastq XXXXX/raw_f.fastq
echo "========== producing quality report for female illumina library =========="
$FASTQC XXXXX/raw_f.fastq -o ${CWD}/QualityReports -d ${CWD}/QualityReports
echo "========== removing female illumina reads mapping to mitochondrial DNA =========="
$BOWTIE -p $CORES $CWD/MitoIndex/MtRef XXXXX/raw_f.fastq --un XXXXX/f.fastq 2> ${illDIR}/f_bowtie.log
cp XXXXX/f.fastq ${illDIR}
	
EOF
sed 's/XXXXX/$TMPDIR/g' ${CWD}/qsubscripts/femalemito.bashX > ${CWD}/qsubscripts/femalemito.bash


cat > ${CWD}/qsubscripts/malemito.bashX <<EOF
#!/bin/bash
#PBS -N redkmer_m_mito
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=24:mem=32gb:tmpspace=500gb
#PBS -e $CWD/reports
#PBS -o $CWD/reports

module load bowtie/1.1.1
module load fastqc
module load anaconda3/personal

cp ${illDIR}/raw_m.fastq XXXXX/raw_m.fastq
echo "========== producing quality report for male illumina library =========="
$FASTQC XXXXX/raw_m.fastq -o ${CWD}/QualityReports -d ${CWD}/QualityReports
echo "========== removing male illumina reads mapping to mitochondrial DNA =========="
$BOWTIE -p $CORES $CWD/MitoIndex/MtRef XXXXX/raw_m.fastq --un XXXXX/m.fastq 2> ${illDIR}/m_bowtie.log
cp XXXXX/m.fastq ${illDIR}
	
EOF
sed 's/XXXXX/$TMPDIR/g' ${CWD}/qsubscripts/malemito.bashX > ${CWD}/qsubscripts/malemito.bash


	MMALEJOB=$(qsub ${CWD}/qsubscripts/malemito.bash)
	echo $MMALEJOB
	MFEMALEJOB=$(qsub ${CWD}/qsubscripts/femalemito.bash)
	echo $MFEMALEJOB	
	
	while qstat $MFEMALEJOB &> /dev/null; do
	    sleep 10;
	done;

	while qstat $MMALEJOB &> /dev/null; do
	    sleep 10;
	done;


fi

#cp ${illDIR}/raw_f.fastq $TMPDIR/raw_f.fastq
#$BOWTIE -p $CORES $CWD/MitoIndex/MtRef $TMPDIR/raw_f.fastq --un $TMPDIR/f.fastq 2> ${illDIR}/f_bowtie.log
#cp $TMPDIR/f.fastq ${illDIR}

#cp ${illDIR}/raw_m.fastq $TMPDIR/raw_m.fastq
#$BOWTIE -p $CORES $CWD/MitoIndex/MtRef $TMPDIR/raw_m.fastq --un $TMPDIR/m.fastq 2> ${illDIR}/m_bowtie.log
#cp $TMPDIR/m.fastq ${illDIR}



printf "======= Done step 1 =======\n"

