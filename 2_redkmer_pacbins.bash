#!/bin/bash

runfile="$1"
if source ${runfile}; then
printf "======= redkmer 0.2 =======\n"
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

mkdir -p $CWD/counts
mkdir -p $CWD/index
mkdir -p $CWD/bins

## Use this to skip sections if needed
#if false; then
## everything below "if false; then" is skipped
#fi
## everything below "fi" is run


# Build the index and map the Illumina data
$BOWTIEB $pacM $CWD/index/m_pac



# Map the Illumina data on the pacbio reads
$BOWTIE -a -p20 -v 0 $CWD/index/m_pac --suppress 1,2,4,5,6,7,8,9 $illF $CWD/counts/female.txt
$BOWTIE -a -p20 -v 0 $CWD/index/m_pac --suppress 1,2,4,5,6,7,8,9 $illM $CWD/counts/male.txt

printf "======= sort and counting files =======\n"

sort $CWD/counts/female.txt | uniq -c > $CWD/counts/female
sort $CWD/counts/male.txt | uniq -c >$CWD/counts/male

printf "======= merging female and male libraries =======\n"

join -a1 -a2 -1 2 -2 2 -o '0,1.1,2.1' -e "0" $CWD/counts/female $CWD/counts/male > $CWD/counts/merge
awk '$4=(($2+0.1)/($3+0.1))' $CWD/counts/merge > $CWD/counts/merge_cq

# Replace space with tabs
awk -v OFS="\t" '$1=$1' $CWD/counts/merge_cq > $CWD/counts/merge_cq.2

printf "======= generating pacBio_MappedReads.txt file  =======\n"

# Add column header
awk 'BEGIN {print "pacbio_read\tfemale\tmale\tCQ"} {print}' $CWD/counts/merge_cq.2 > $CWD/counts/pacBio_MappedReads.txt

printf "======= creating chromosomal bins of pacbio reads =======\n"

# Create chromosomal bins
awk '{if($4>1.5) print $1}' $CWD/counts/merge_cq.2 > $CWD/bins/X_reads
awk '{if($4<1.5 && $4>0.2) print $1}' $CWD/counts/merge_cq.2 > $CWD/bins/A_reads
awk '{if($4<0.2) print $1}' $CWD/counts/merge_cq.2 > $CWD/bins/Y_reads

# Get sequences of pacBio bins

cat $CWD/bins/X_reads | xargs $SAMTOOLS faidx $pacM > $CWD/bins/X_fasta
cat $CWD/bins/A_reads | xargs $SAMTOOLS faidx $pacM > $CWD/bins/A_fasta
cat $CWD/bins/Y_reads | xargs $SAMTOOLS faidx $pacM > $CWD/bins/Y_fasta

printf "======= cleaning up  =======\n"

rm $CWD/counts/merge
rm $CWD/counts/merge_cq
rm $CWD/counts/merge_cq.2
rm $CWD/counts/female.txt
rm $CWD/counts/male.txt
rm $CWD/counts/female
rm $CWD/counts/male
rm $CWD/bins/X_reads
rm $CWD/bins/A_reads
rm $CWD/bins/Y_reads

