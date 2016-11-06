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

mkdir -p $CWD/pacBio_illmapping
mkdir -p $CWD/pacBio_illmapping/logs
mkdir -p $CWD/pacBio_illmapping/mapping_rawdata
mkdir -p $CWD/pacBio_illmapping/index
mkdir -p $CWD/pacBio_bins
mkdir -p $CWD/pacBio_bins/fasta
mkdir -p $CWD/temp

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

$BOWTIEB $pacM $CWD/pacBio_illmapping/index/m_pac

printf "======= Mapping illumina reads to pacBIO reads =======\n"

$BOWTIE -a -t -p$CORES -v 0 $CWD/pacBio_illmapping/index/m_pac --suppress 1,2,4,5,6,7,8,9 $illF 1> $CWD/pacBio_illmapping/mapping_rawdata/female.txt 2> $CWD/pacBio_illmapping/logs/female_log.txt
$BOWTIE -a -t -p$CORES -v 0 $CWD/pacBio_illmapping/index/m_pac --suppress 1,2,4,5,6,7,8,9 $illM 1> $CWD/pacBio_illmapping/mapping_rawdata/male.txt 2> $CWD/pacBio_illmapping/logs/male_log.txt

printf "======= sort and counting files =======\n"

time sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G $CWD/pacBio_illmapping/mapping_rawdata/female.txt | uniq -c > $CWD/pacBio_illmapping/mapping_rawdata/female_uniq &
time sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G $CWD/pacBio_illmapping/mapping_rawdata/male.txt | uniq -c > $CWD/pacBio_illmapping/mapping_rawdata/male_uniq

#old sort
#sort -k1b,1 $CWD/pacBio_illmapping/mapping_rawdata/female.txt | uniq -c > $CWD/pacBio_illmapping/mapping_rawdata/female_uniq
#sort -k1b,1 $CWD/pacBio_illmapping/mapping_rawdata/male.txt | uniq -c > $CWD/pacBio_illmapping/mapping_rawdata/male_uniq

printf "======= merging female and male pacBio_illmapping =======\n"

join -a1 -a2 -1 2 -2 2 -o '0,1.1,2.1' -e "0" $CWD/pacBio_illmapping/mapping_rawdata/female_uniq $CWD/pacBio_illmapping/mapping_rawdata/male_uniq > $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= normalizing to library size =======\n"
awk -v ma="$illLIBMsize" -v fema="$illLIBFsize" -v le="$illnorm" '{print $1, ($2*fema/le), ($3*ma/le)}' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= calculating CQ of pacBIO reads =======\n"

awk '{_div1= $3 ? ($2 / $3) : 0 ; print $0, _div1 }' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= calculating sum of pacBio_illmapping on pacBIO reads =======\n"

awk '{print $0, ($2+$3)}' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

# Replace space with tabs
awk -v OFS="\t" '$1=$1' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= generating pacBio_MappedReads.txt file  =======\n"

# Add column header
awk 'BEGIN {print "pacbio_read\tfemale\tmale\tCQ\tSum"} {print}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_illmapping/pacBio_MappedReads.txt

printf "======= creating chromosomal bins of pacbio reads =======\n"

awk '{if($4>=1.5) print $1}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_bins/X_reads
awk '{if($4<1.5 && $4>0.2) print $1}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_bins/A_reads
awk '{if($4<0.2) print $1}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_bins/Y_reads

# Get sequences of pacBio bins
rm -f $CWD/testreadspac/m_pac.fasta.fai
cat $CWD/pacBio_bins/X_reads | xargs $SAMTOOLS faidx $pacM > $CWD/pacBio_bins/fasta/Xbin.fasta
cat $CWD/pacBio_bins/A_reads | xargs $SAMTOOLS faidx $pacM > $CWD/pacBio_bins/fasta/Abin.fasta
cat $CWD/pacBio_bins/Y_reads | xargs $SAMTOOLS faidx $pacM > $CWD/pacBio_bins/fasta/Ybin.fasta

printf "======= done  step 2 =======\n"
