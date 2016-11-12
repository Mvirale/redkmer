!/bin/bash

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

mkdir -p $CWD/kmers/offtargets
mkdir -p $CWD/plots
mkdir -p $CWD/kmers/bowtie
mkdir -p $CWD/kmers/bowtie/index
mkdir -p $CWD/kmers/bowtie/mapping
mkdir -p $CWD/kmers/bowtie/offtargets
mkdir -p $CWD/kmers/bowtie/offtargets/logs


printf "======= performing off-target analysis at 90perc (3bp mismatch) =======\n"
#$BOWTIEB $CWD/pacBio_bins/fasta/Xbin.fasta $CWD/kmers/bowtie/index/Xbin &
#$BOWTIEB $CWD/pacBio_bins/fasta/Abin.fasta $CWD/kmers/bowtie/index/Abin &
#$BOWTIEB $CWD/pacBio_bins/fasta/Ybin.fasta $CWD/kmers/bowtie/index/Ybin &
#$BOWTIEB $CWD/pacBio_bins/fasta/GAbin.fasta $CWD/kmers/bowtie/index/GAbin &

wait $(jobs -p)


printf "======= Mapping illumina reads to pacBIO reads =======\n"

$BOWTIE -a -t -p$CORES -v 2 $CWD/kmers/bowtie/index/Abin --suppress 2,3,4,5,6,7,8,9 -f $CWD/kmers/fasta/Xkmers.fasta 1> $CWD/kmers/bowtie/offtargets/Abin.txt 2> $CWD/kmers/bowtie/offtargets/logs/Abin_log.txt
$BOWTIE -a -t -p$CORES -v 2 $CWD/kmers/bowtie/index/Ybin --suppress 2,3,4,5,6,7,8,9 -f $CWD/kmers/fasta/Xkmers.fasta 1> $CWD/kmers/bowtie/offtargets/Ybin.txt 2> $CWD/kmers/bowtie/offtargets/logs/Ybin_log.txt
$BOWTIE -a -t -p$CORES -v 2 $CWD/kmers/bowtie/index/GAbin --suppress 2,3,4,5,6,7,8,9 -f $CWD/kmers/fasta/Xkmers.fasta 1> $CWD/kmers/bowtie/offtargets/GAbin.txt 2> $CWD/kmers/bowtie/offtargets/logs/GAbin_log.txt


cat $CWD/kmers/bowtie/offtargets/Abin.txt $CWD/kmers/bowtie/offtargets/Ybin.txt $CWD/kmers/bowtie/offtargets/GAbin.txt | sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G | uniq -c  | awk '{print $2, $1}' >  $CWD/kmers/bowtie/offtargets/kmer_hits_AYGAbin

sed '1d' $CWD/kmers/rawdata/kmers_blastresults > tmpfile
join -a1 -a2 -1 1 -2 1 -o '0,2.2,2.3,2.4,2.5,2.6,2.7,1.2' -e "0" $CWD/kmers/bowtie/offtargets/kmer_hits_AYGAbin tmpfile > $CWD/kmers/kmer_results.txt; rm tmpfile

printf "======= generating kmers_all_results_withofftargets file =======\n"

awk -v OFS="\t" '$1=$1' $CWD/kmers/kmer_results.txt > tmpfile; mv tmpfile $CWD/kmers/kmer_results.txt

#Add column header
awk 'BEGIN {print "kmer_id\tseq\tfemale\tmale\tCQ\tsum\tbin\tofftargets"} {print}' $CWD/kmers/kmer_results.txt > tmpfile; mv tmpfile $CWD/kmers/kmer_results.txt

printf "======= done step 5 =======\n"
