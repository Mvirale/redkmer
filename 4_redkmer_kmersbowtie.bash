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

source ${BASEDIR}/redkmer.cfg


kmers=$CWD/kmers/fasta/allkmers.fasta

printf "======= Creating bowtie index for PacBio bins =======\n"

$BOWTIEB $CWD/pacBio_bins/fasta/Xbin.fasta $CWD/kmers/bowtie/index/Xbin &
$BOWTIEB $CWD/pacBio_bins/fasta/Abin.fasta $CWD/kmers/bowtie/index/Abin &
$BOWTIEB $CWD/pacBio_bins/fasta/Ybin.fasta $CWD/kmers/bowtie/index/Ybin &
if [ -s "$CWD/pacBio_bins/fasta/GAbin.fasta" ];then
$BOWTIEB $CWD/pacBio_bins/fasta/GAbin.fasta $CWD/kmers/bowtie/index/GAbin & 
fi

wait $(jobs -p) 

printf "======= Running bowtie against X chromosome bin =======\n"
$BOWTIE -a -t -p$CORES -v 0 $CWD/kmers/bowtie/index/Xbin --suppress 2,3,4,5,6,7,8,9 -f $kmers  1> $CWD/kmers/bowtie/mapping/Xbin.txt 2> $CWD/kmers/bowtie/mapping/logs/Xbin_log.txt

printf "======= Running bowtie against autosome bin =======\n"
$BOWTIE -a -t -p$CORES -v 0 $CWD/kmers/bowtie/index/Abin --suppress 2,3,4,5,6,7,8,9 -f $kmers 1> $CWD/kmers/bowtie/mapping/Abin.txt 2> $CWD/kmers/bowtie/mapping/logs/Abin_log.txt

printf "======= Running bowtie against Y chromosome bin =======\n"
$BOWTIE -a -t -p$CORES -v 0 $CWD/kmers/bowtie/index/Ybin --suppress 2,3,4,5,6,7,8,9 -f $kmers 1> $CWD/kmers/bowtie/mapping/Ybin.txt 2> $CWD/kmers/bowtie/mapping/logs/Ybin_log.txt

printf "======= Running bowtie against GA chromosome bin =======\n"
if [ -s "$CWD/pacBio_bins/fasta/GAbin.fasta" ];then
$BOWTIE -a -t -p$CORES -v 0 $CWD/kmers/bowtie/index/GAbin --suppress 2,3,4,5,6,7,8,9 -f $kmers 1> $CWD/kmers/bowtie/mapping/GAbin.txt 2> $CWD/kmers/bowtie/mapping/logs/GAbin_log.txt
else
touch $CWD/kmers/bowtie/mapping/GAbin.txt
fi

printf "======= extracting blast results =======\n"

sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G $CWD/kmers/bowtie/mapping/Xbin.txt | uniq -c  | awk '{print $2, $1}' >  $CWD/kmers/bowtie/mapping/kmer_hits_Xbin &
sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G $CWD/kmers/bowtie/mapping/Abin.txt | uniq -c  | awk '{print $2, $1}' >  $CWD/kmers/bowtie/mapping/kmer_hits_Abin &
sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G $CWD/kmers/bowtie/mapping/Ybin.txt | uniq -c  | awk '{print $2, $1}' >  $CWD/kmers/bowtie/mapping/kmer_hits_Ybin &
sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G $CWD/kmers/bowtie/mapping/GAbin.txt | uniq -c  | awk '{print $2, $1}' >  $CWD/kmers/bowtie/mapping/kmer_hits_GAbin &

wait $(jobs -p)


join -a1 -a2 -1 1 -2 1 -o '0,1.2,2.2' -e "0" $CWD/kmers/bowtie/mapping/kmer_hits_Xbin $CWD/kmers/bowtie/mapping/kmer_hits_Abin > $CWD/kmers/bowtie/mapping/kmer_hits_XAbin
join -a1 -a2 -1 1 -2 1 -o '0,1.2,1.3,2.2' -e "0" $CWD/kmers/bowtie/mapping/kmer_hits_XAbin $CWD/kmers/bowtie/mapping/kmer_hits_Ybin > $CWD/kmers/bowtie/mapping/kmer_hits_XAYbin
join -a1 -a2 -1 1 -2 1 -o '0,1.2,1.3,1.4,2.2' -e "0" $CWD/kmers/bowtie/mapping/kmer_hits_XAYbin $CWD/kmers/bowtie/mapping/kmer_hits_GAbin > $CWD/kmers/bowtie/mapping/kmer_hits_bins

rm $CWD/kmers/bowtie/mapping/kmer_hits_XAbin
rm $CWD/kmers/bowtie/mapping/kmer_hits_XAYbin

awk '{print $0, ($2+$3+$4+$5)}' $CWD/kmers/bowtie/mapping/kmer_hits_bins > tmpfile; mv tmpfile $CWD/kmers/bowtie/mapping/kmer_hits_bins
awk '{print $0, ($2/$6)}' $CWD/kmers/bowtie/mapping/kmer_hits_bins > tmpfile; mv tmpfile $CWD/kmers/bowtie/mapping/kmer_hits_bins


printf "======= merging bowtie bin results to kmer_counts data =======\n"

sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G $CWD/kmers/bowtie/mapping/kmer_hits_bins > tmpfile1; mv tmpfile1 $CWD/kmers/bowtie/mapping/kmer_hits_bins & 
sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G $CWD/kmers/rawdata/kmers_to_merge > tmpfile2; mv tmpfile2 $CWD/kmers/rawdata/kmers_to_merge & 

wait $(jobs -p)


join -a1 -a2 -1 1 -2 1 -o '0,2.2,2.3,2.4,2.5,2.6,1.2,1.3,1.4,1.5,1.6,1.7' -e "0"  $CWD/kmers/bowtie/mapping/kmer_hits_bins $CWD/kmers/rawdata/kmers_to_merge > $CWD/kmers/rawdata/kmers_hits_results
awk '{print $0, "0"}'  $CWD/kmers/rawdata/kmers_hits_results > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_hits_results 
awk -v xsi="$XSI" '{if ($12>=xsi) {$13="pass"}; print}' $CWD/kmers/rawdata/kmers_hits_results > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_hits_results
awk -v xsi="$XSI" '{if ($12<xsi) {$13="fail"}; print}' $CWD/kmers/rawdata/kmers_hits_results > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_hits_results
awk '{if ($11==0) {$13="nohits"}; print}' $CWD/kmers/rawdata/kmers_hits_results > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_hits_results



printf "======= generating Xkmers.fasta file for off-target analysis =======\n"

awk '{if ($13=="pass") print $1, $2}' $CWD/kmers/rawdata/kmers_hits_results |awk '{print ">"$1"\n"$2}' > $CWD/kmers/fasta/Xkmers.fasta

printf "======= generating kmers_all_results file =======\n"

awk -v OFS="\t" '$1=$1' $CWD/kmers/rawdata/kmers_hits_results > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_hits_results

#Add column header
awk 'BEGIN {print "kmer_id\tseq\tfemale\tmale\tCQ\tsum\thits_X\thits_A\thits_Y\thits_GA\thits_sum\tperchitsX\thits_threshold"} {print}' $CWD/kmers/rawdata/kmers_hits_results > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_hits_results

printf "======= done step 4 =======\n"


if [[ "$EMAIL" -eq "1" ]]; then
echo "done step 4 ${runfile}" > done4.txt
sudo ssmtp $emailaddress < done4.txt
rm done4.txt
else
echo "done step 4"
fi


