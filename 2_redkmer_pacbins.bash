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
mkdir -p $CWD/counts/logs
mkdir -p $CWD/counts/mapping_rawdata
mkdir -p $CWD/index
mkdir -p $CWD/bins

## Use this to skip sections if needed
#if false; then
## everything below "if false; then" is skipped
#fi
## everything below "fi" is run

printf "======= calculating library sizes =======\n"

illLIBMsize=$(wc -l $illM | awk '{print ($1/4)}')
illLIBFsize=$(wc -l $illF | awk '{print ($1/4)}')
illnorm=$((($illLIBMsize+$illLIBFsize)/2))

printf "======= Building index of pacBIO reads =======\n"

$BOWTIEB $pacM $CWD/index/m_pac

printf "======= Mapping illumina reads to pacBIO reads =======\n"

$BOWTIE -a -t -p$CORES -v 0 $CWD/index/m_pac --suppress 1,2,4,5,6,7,8,9 $illF 1> $CWD/counts/mapping_rawdata/female.txt 2> $CWD/counts/logs/female_log.txt
$BOWTIE -a -t -p$CORES -v 0 $CWD/index/m_pac --suppress 1,2,4,5,6,7,8,9 $illM 1> $CWD/counts/mapping_rawdata/male.txt 2> $CWD/counts/logs/male_log.txt

printf "======= sort and counting files =======\n"

sort -k1b,1 $CWD/counts/mapping_rawdata/female.txt | uniq -c > $CWD/counts/mapping_rawdata/female_uniq
sort -k1b,1 $CWD/counts/mapping_rawdata/male.txt | uniq -c > $CWD/counts/mapping_rawdata/male_uniq

printf "======= merging female and male counts =======\n"

join -a1 -a2 -1 2 -2 2 -o '0,1.1,2.1' -e "0" $CWD/counts/mapping_rawdata/female_uniq $CWD/counts/mapping_rawdata/male_uniq > $CWD/counts/mapping_rawdata/merge

printf "======= normalizing to library size =======\n"
awk -v ma="$illLIBMsize" -v fema="$illLIBFsize" -v le="$illnorm" '{print $1, ($2*fema/le), ($3*ma/le)}' $CWD/counts/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/counts/mapping_rawdata/merge

printf "======= calculating CQ of pacBIO reads =======\n"

awk '{_div1= $3 ? ($2 / $3) : 0 ; print $0, _div1 }' $CWD/counts/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/counts/mapping_rawdata/merge

printf "======= calculating sum of counts on pacBIO reads =======\n"

awk '{print $0, ($2+$3)}' $CWD/counts/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/counts/mapping_rawdata/merge

# Replace space with tabs
awk -v OFS="\t" '$1=$1' $CWD/counts/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/counts/mapping_rawdata/merge

printf "======= generating pacBio_MappedReads.txt file  =======\n"

# Add column header
awk 'BEGIN {print "pacbio_read\tfemale\tmale\tCQ\tSum"} {print}' $CWD/counts/mapping_rawdata/merge > $CWD/counts/pacBio_MappedReads.txt

printf "======= creating chromosomal bins of pacbio reads =======\n"

awk '{if($4>=1.5) print $1}' $CWD/counts/mapping_rawdata/merge > $CWD/bins/X_reads
awk '{if($4<1.5 && $4>0.2) print $1}' $CWD/counts/mapping_rawdata/merge > $CWD/bins/A_reads
awk '{if($4<0.2) print $1}' $CWD/counts/mapping_rawdata/merge > $CWD/bins/Y_reads

# Get sequences of pacBio bins
rm -f $CWD/testreadspac/m_pac.fasta.fai
cat $CWD/bins/X_reads | xargs $SAMTOOLS faidx $pacM > $CWD/bins/X_fasta
cat $CWD/bins/A_reads | xargs $SAMTOOLS faidx $pacM > $CWD/bins/A_fasta
cat $CWD/bins/Y_reads | xargs $SAMTOOLS faidx $pacM > $CWD/bins/Y_fasta

printf "======= done  step 2 =======\n"
