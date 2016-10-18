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

mkdir -p $CWD/kmers/offtargets
mkdir -p $CWD/plots

printf "======= performing off-target analysis at 90perc (3bp mismatch) =======\n"


$BLAST -db $CWD/kmers/blast/index/blastdb_Abin -query $CWD/kmers/fasta/Xkmers.fasta -out $CWD/kmers/offtargets/Xkmers_offtargets_Abin -perc_identity 90 -outfmt 6 -num_threads $CORES
$BLAST -db $CWD/kmers/blast/index/blastdb_Ybin -query $CWD/kmers/fasta/Xkmers.fasta -out $CWD/kmers/offtargets/Xkmers_offtargets_Ybin -perc_identity 90 -outfmt 6 -num_threads $CORES

awk '{print $1}' $CWD/kmers/offtargets/Xkmers_offtargets_Abin  >  $CWD/kmers/offtargets/kmer_hits_Abin
awk '{print $1}' $CWD/kmers/offtargets/Xkmers_offtargets_Ybin  >  $CWD/kmers/offtargets/kmer_hits_Ybin

cat $CWD/kmers/offtargets/kmer_hits_Abin $CWD/kmers/offtargets/kmer_hits_Ybin | sort -k1b,1 | uniq -c  | awk '{print $2, $1}' >  $CWD/kmers/offtargets/kmer_hits_AYbin

sed '1d' $CWD/kmers/rawdata/kmers_blastresults > tmpfile
join -a1 -a2 -1 1 -2 1 -o '0,2.2,2.3,2.4,2.5,2.6,2.7,1.2' -e "0" $CWD/kmers/offtargets/kmer_hits_AYbin tmpfile > $CWD/kmers/kmer_results.txt; rm tmpfile

printf "======= generating kmers_all_results_withofftargets file =======\n"

awk -v OFS="\t" '$1=$1' $CWD/kmers/kmer_results.txt > tmpfile; mv tmpfile $CWD/kmers/kmer_results.txt

#Add column header
awk 'BEGIN {print "kmer_id\tseq\tfemale\tmale\tCQ\tsum\tbin\tofftargets"} {print}' $CWD/kmers/kmer_results.txt > tmpfile; mv tmpfile $CWD/kmers/kmer_results.txt

printf "======= done step 5 =======\n"



	