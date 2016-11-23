#!/bin/bash


runfile="$1"
if source ${runfile}; then
printf "======= redkmer v1.0 =======\n"
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


printf "======= calculating library sizes =======\n"

illLIBMsize=$(wc -l $illM | awk '{print ($1/4)}')
illLIBFsize=$(wc -l $illF | awk '{print ($1/4)}')
illnorm=$((($illLIBMsize+$illLIBFsize)/2))

printf "======= Building index of pacBIO reads & Mapping illumina reads to pacBIO reads =======\n"

#$BOWTIEB $pacM $CWD/pacBio_illmapping/index/m_pac

#cp $pacM ${CWD}/${pacDIR}/m_pac_multiline.fasta
#awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < $pacM > ${CWD}/${pacDIR}/m_pac_online.fasta
#sed '1d' ${CWD}/${pacDIR}/m_pac_online.fasta > $pacM

rm -f $CWD/pacBio_illmapping/logs/female_log.txt
rm -f $CWD/pacBio_illmapping/logs/male_log.txt
rm -f $CWD/pacBio_illmapping/mapping_rawdata/female.txt
rm -f $CWD/pacBio_illmapping/mapping_rawdata/male.txt

while read pacbioname && read pacbioseq
do

nextpacbio=$pacbioname$'\n'$pacbioseq
#echo "$nextpacbio"

(
$BOWTIEB -c "$nextpacbio" $CWD/pacBio_illmapping/index/$pacbioname
$BOWTIE -a -t -p 2 -v 0 $CWD/pacBio_illmapping/index/$pacbioname --suppress 1,2,4,5,6,7,8,9 $illF 1>> $CWD/pacBio_illmapping/mapping_rawdata/female.txt 2>> $CWD/pacBio_illmapping/logs/female_log.txt
$BOWTIE -a -t -p 2 -v 0 $CWD/pacBio_illmapping/index/$pacbioname --suppress 1,2,4,5,6,7,8,9 $illM 1>> $CWD/pacBio_illmapping/mapping_rawdata/male.txt 2>> $CWD/pacBio_illmapping/logs/male_log.txt
rm $CWD/pacBio_illmapping/index/$pacbioname
)&

#$BOWTIEB -c $linekmerseq $CWD/pacBio_illmapping/index/m_pac
#$BOWTIE -a -t -p $CORES -v 0 $CWD/pacBio_illmapping/index/m_pac --suppress 1,2,4,5,6,7,8,9 $illF 1>> $CWD/pacBio_illmapping/mapping_rawdata/female.txt 2>> $CWD/pacBio_illmapping/logs/female_log.txt
#$BOWTIE -a -t -p $CORES -v 0 $CWD/pacBio_illmapping/index/m_pac --suppress 1,2,4,5,6,7,8,9 $illM 1>> $CWD/pacBio_illmapping/mapping_rawdata/male.txt 2>> $CWD/pacBio_illmapping/logs/male_log.txt

done < $pacM

exit 1

printf "======= sort and counting files =======\n"

time sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G $CWD/pacBio_illmapping/mapping_rawdata/female.txt | uniq -c > $CWD/pacBio_illmapping/mapping_rawdata/female_uniq
time sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G $CWD/pacBio_illmapping/mapping_rawdata/male.txt | uniq -c > $CWD/pacBio_illmapping/mapping_rawdata/male_uniq

printf "======= merging female and male pacBio_illmapping =======\n"

join -a1 -a2 -1 2 -2 2 -o '0,1.1,2.1' -e "0" $CWD/pacBio_illmapping/mapping_rawdata/female_uniq $CWD/pacBio_illmapping/mapping_rawdata/male_uniq > $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= normalizing to library size =======\n"
awk -v ma="$illLIBMsize" -v fema="$illLIBFsize" -v le="$illnorm" '{print $1, ($2*fema/le), ($3*ma/le)}' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= calculating CQ of pacBIO reads =======\n"

awk '{print $0, (($2+1)/($3+1))}' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= calculating sum of pacBio_illmapping on pacBIO reads =======\n"

awk '{print $0, ($2+$3)}' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= calculating LSum (Sum/length of PBreads * median PBread length  =======\n"

rm -f $pacM.fai
$SAMTOOLS faidx $pacM
awk '{print $1, $2}' $pacM.fai | sort -k1b,1 > $pacM.lengths
join -a1 -a2 -1 1 -1 1 -o'0,2.2,1.2,1.3,1.4,1.5' -e "0" $CWD/pacBio_illmapping/mapping_rawdata/merge $pacM.lengths > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

medianlength=$(awk '{print $2}' $pacM.lengths | sort -n | awk '
  BEGIN {
    c = 0;
    sum = 0;
  }
  $1 ~ /^[0-9]*(\.[0-9]*)?$/ {
    a[c++] = $1;
    sum += $1;
  }
  END {
    ave = sum / c;
    if( (c % 2) == 1 ) {
      median = a[ int(c/2) ];
    } else {
      median = ( a[c/2] + a[c/2-1] ) / 2;
    }
    OFS="\t";
    print median;
  }
')

awk -v ml="$medianlength" '{print $0, ($6 / $2 * ml)}' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= filter length and LSum (>=2000bp and LSum>=50)  =======\n"

awk -v pl="$pac_length" '{if($2>=pl)print $0}' $CWD/pacBio_illmapping/mapping_rawdata/merge | awk -v ls="$LSum" '{if ($7>=ls) print $0}' > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge 

# Replace space with tabs
awk -v OFS="\t" '$1=$1' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= generating pacBio_MappedReads.txt file  =======\n"

# Add column header
awk 'BEGIN {print "pacbio_read\tbp\tfemale\tmale\tCQ\tSum\tLSum"} {print}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_illmapping/pacBio_MappedReads.txt

printf "======= creating chromosomal bins of pacbio reads =======\n"

awk '{if($5>=1.5 && $5<10) print $1}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_bins/X_reads
awk '{if($5<1.5 && $5>0.2) print $1}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_bins/A_reads
awk '{if($5<0.2) print $1}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_bins/Y_reads
awk '{if($5>=10) print $1}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_bins/GA_reads

# Get sequences of pacBio bins

cat $CWD/pacBio_bins/X_reads | xargs $SAMTOOLS faidx $pacM > $CWD/pacBio_bins/fasta/Xbin.fasta
cat $CWD/pacBio_bins/A_reads | xargs $SAMTOOLS faidx $pacM > $CWD/pacBio_bins/fasta/Abin.fasta
cat $CWD/pacBio_bins/Y_reads | xargs $SAMTOOLS faidx $pacM > $CWD/pacBio_bins/fasta/Ybin.fasta
cat $CWD/pacBio_bins/GA_reads | xargs $SAMTOOLS faidx $pacM > $CWD/pacBio_bins/fasta/GAbin.fasta


printf "======= done  step 2 =======\n"

if [[ "$EMAIL" -eq "1" ]]; then
echo "done step 2 ${runfile}" > done2.txt
sudo ssmtp $emailaddress < done2.txt
rm done2.txt
else
echo "done step 2"
fi