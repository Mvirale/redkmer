#!/bin/bash
#PBS -N redkmer_step1
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=20:mem=16gb
#PBS -e /home/nikiwind/reports
#PBS -o /home/nikiwind/reports


if [ -z ${PBS_ENVIRONMENT+x} ]
then
echo "---> running on the Perugia numbercruncher..."
source redkmer.cfg
else
echo "---> running on HPC cluster..."
source $PBS_O_WORKDIR/redkmer.cfg
module load samtools
module load bowtie/1.1.1
fi

mkdir -p $CWD/read_analysis
mkdir -p $CWD/read_analysis/database
mkdir -p $CWD/read_analysis/database/index
mkdir -p $CWD/read_analysis/database/mapping/sam
mkdir -p $CWD/read_analysis/database/mapping/counts
mkdir -p $CWD/read_analysis/database/mapping/bed

printf "======= calculating library sizes =======\n"

illLIBMsize=$(wc -l $illM | awk '{print ($1/4)}')
illLIBFsize=$(wc -l $illF | awk '{print ($1/4)}')
illnorm=$((($illLIBMsize+$illLIBFsize)/2))

####~~~~~~~~ IMPORTANT:~~~~~~~~~~~####
####~~~~~~~~ REQUIRES FILE read2check.list containing names of pacBIO reads to be checked####
 
 printf "======= 1. extracting fasta from list =======\n"

cat $CWD/read_analysis/reads2check.list | xargs $SAMTOOLS faidx $pacM > $CWD/read_analysis/database/reads2check.fasta
 
printf "======= 2. building bowtie database =======\n"
$BOWTIEB $CWD/read_analysis/database/reads2check.fasta $CWD/read_analysis/database/index/reads2check

printf "======= 3. mapping data from females =======\n"
 
$BOWTIE -a -t -p$CORES -v 0 $CWD/read_analysis/database/index/reads2check $illF -S $CWD/read_analysis/database/mapping/sam/female.sam
 
printf "======= 4. converting sam into bed and deleting files for females =======\n"
 
 $SAMTOOLS view -bS $CWD/read_analysis/database/mapping/sam/female.sam | $SAMTOOLS sort - $CWD/read_analysis/database/mapping/sam/female
 $SAMTOOLS index $CWD/read_analysis/database/mapping/sam/female.bam
 $SAMTOOLS idxstats $CWD/read_analysis/database/mapping/sam/female.bam > $CWD/read_analysis/database/mapping/counts/female.counts
 rm $CWD/read_analysis/database/mapping/sam/female.sam
 $BEDTOOLS genomecov -ibam $CWD/read_analysis/database/mapping/sam/female.bam -d > $CWD/read_analysis/database/mapping/bed/female.bed
 
 
printf "======= 5. mapping data from males =======\n"
 
 $BOWTIE -a -t -p$CORES -v 0 $CWD/read_analysis/database/index/reads2check $illM -S $CWD/read_analysis/database/mapping/sam/male.sam
 
printf "======= 6. converting sam into bed and deleting files for males =======\n"
 
 $SAMTOOLS view -bS $CWD/read_analysis/database/mapping/sam/male.sam | $SAMTOOLS sort - $CWD/read_analysis/database/mapping/sam/male
 $SAMTOOLS index $CWD/read_analysis/database/mapping/sam/male.bam
 $SAMTOOLS idxstats $CWD/read_analysis/database/mapping/sam/male.bam > $CWD/read_analysis/database/mapping/counts/male.counts
 rm $CWD/read_analysis/database/mapping/sam/male.sam
$BEDTOOLS genomecov -ibam $CWD/read_analysis/database/mapping/sam/male.bam -d > $CWD/read_analysis/database/mapping/bed/male.bed
 

printf "======= 7. combining bed files =======\n"

 awk '{print $1"_"$2,$0}' $CWD/read_analysis/database/mapping/bed/female.bed > tmpfile; mv tmpfile $CWD/read_analysis/database/mapping/bed/female.bed
 awk '{print $1"_"$2,$0}' $CWD/read_analysis/database/mapping/bed/male.bed > tmpfile; mv tmpfile $CWD/read_analysis/database/mapping/bed/male.bed

 join -a1 -a2 -1 1 -2 1 -o '1.2,1.3,1.4,2.4' $CWD/read_analysis/database/mapping/bed/female.bed $CWD/read_analysis/database/mapping/bed/male.bed > $CWD/read_analysis/database/mapping/bed/merge.bed

printf "======= normalizing to library size =======\n"
 awk -v ma="$illLIBMsize" -v fema="$illLIBFsize" -v le="$illnorm" '{print $1,$2,($3*fema/le),($4*ma/le)}' $CWD/read_analysis/database/mapping/bed/merge.bed > tmpfile; mv tmpfile $CWD/read_analysis/database/mapping/bed/merge.bed
 
printf "======= calculating diff of coverage =======\n"
 
 awk '{print $0, ($4-$3)}' $CWD/read_analysis/database/mapping/bed/merge.bed > tmpfile; mv tmpfile $CWD/read_analysis/database/mapping/bed/merge.bed

printf "======= bring female to negative scale =======\n"

 awk '{print $1,$2,($3*-1),$4,$5}' $CWD/read_analysis/database/mapping/bed/merge.bed > tmpfile; mv tmpfile $CWD/read_analysis/database/mapping/bed/merge.bed

# Replace space with tabs
awk -v OFS="\t" '$1=$1' $CWD/read_analysis/database/mapping/bed/merge.bed > tmpfile; mv tmpfile $CWD/read_analysis/database/mapping/bed/merge.bed
 
printf "======= generating reads2check_spatial.txt file  =======\n"
 
# # Add column header
awk 'BEGIN {print "pacbio_read\tpos\tfemale\tmale\tdiff"} {print}' $CWD/read_analysis/database/mapping/bed/merge.bed > $CWD/read_analysis/reads2check_spatial.txt

printf "======= done  step 7 =======\n"

