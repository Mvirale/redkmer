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


printf "======= Running offtargets analysis against autosome bin 2mismatches =======\n"
$BOWTIE -a -t -p$CORES -v 2 $CWD/kmers/bowtie/index/Abin --suppress 2,3,4,5,6,7,8,9 -f $CWD/kmers/fasta/Xkmers.fasta 1> $CWD/kmers/bowtie/offtargets/Abin.txt 2> $CWD/kmers/bowtie/offtargets/logs/Abin_log.txt

printf "======= Running offtargets analysis against Y chromosome bin 2mismatches =======\n"
$BOWTIE -a -t -p$CORES -v 2 $CWD/kmers/bowtie/index/Ybin --suppress 2,3,4,5,6,7,8,9 -f $CWD/kmers/fasta/Xkmers.fasta 1> $CWD/kmers/bowtie/offtargets/Ybin.txt 2> $CWD/kmers/bowtie/offtargets/logs/Ybin_log.txt

printf "======= Running offtargets analysis against GA bin 2mismatches =======\n"
if [ -s "$CWD/pacBio_bins/fasta/GAbin.fasta" ];then
$BOWTIE -a -t -p$CORES -v 2 $CWD/kmers/bowtie/index/GAbin --suppress 2,3,4,5,6,7,8,9 -f $CWD/kmers/fasta/Xkmers.fasta 1> $CWD/kmers/bowtie/offtargets/GAbin.txt 2> $CWD/kmers/bowtie/offtargets/logs/GAbin_log.txt
else
touch $CWD/kmers/bowtie/offtargets/GAbin.txt
fi

cat $CWD/kmers/bowtie/offtargets/Abin.txt $CWD/kmers/bowtie/offtargets/Ybin.txt $CWD/kmers/bowtie/offtargets/GAbin.txt | sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G | uniq -c  | awk '{print $2, $1}' >  $CWD/kmers/bowtie/offtargets/kmer_hits_AYGAbin

sed '1d' $CWD/kmers/rawdata/kmers_hits_results > tmpfile
join -a1 -a2 -1 1 -2 1 -o '0,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,1.2' -e "0" $CWD/kmers/bowtie/offtargets/kmer_hits_AYGAbin tmpfile > $CWD/kmers/kmer_results.txt; rm tmpfile

printf "======= generating kmers_all_results_withofftargets file =======\n"

awk -v OFS="\t" '$1=$1' $CWD/kmers/kmer_results.txt > tmpfile; mv tmpfile $CWD/kmers/kmer_results.txt

#Add column header
awk 'BEGIN {print "kmer_id\tseq\tfemale\tmale\tCQ\tsum\thits_X\thits_A\thits_Y\thits_GA\thits_sum\tperchitsX\thits_threshold\tofftargets"} {print}' $CWD/kmers/kmer_results.txt > tmpfile; mv tmpfile $CWD/kmers/kmer_results.txt

printf "======= done step 5 =======\n"

if [[ "$EMAIL" -eq "1" ]]; then
echo "done step 5 ${runfile}" > done5.txt
sudo ssmtp $emailaddress < done5.txt
rm done5.txt
else
echo "done step 5"
fi



