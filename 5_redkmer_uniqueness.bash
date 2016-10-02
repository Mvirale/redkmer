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

mkdir -p $CWD/kmers/scores

kmers=$CWD/kmers/candidateXkmers.fasta

#Add a line with kmer_increasing_num to name the kmers
awk '{print ">"$1"\n"$2}' $CWD/kmers/candidateXkmers.seq > $CWD/kmers/candidateXkmers.fasta


# $BLAST_DB -in $kmers -dbtype nucl -out $CWD/blast/index/blastdb_candidateXkmers
# 
# $BLAST -db $CWD/blast/index/blastdb_Abin -query $kmers -out $CWD/blast/hits_Abin -perc_identity 100 -outfmt 6 -num_threads 2

$BLAST -db $CWD/blast/index/blastdb_Abin -query $kmers -out $CWD/kmers/scores/hits_Abin -perc_identity 90 -outfmt 6 -num_threads 2
$BLAST -db $CWD/blast/index/blastdb_Ybin -query $kmers -out $CWD/kmers/scores/hits_Ybin -perc_identity 90 -outfmt 6 -num_threads 2

awk '{print $1}' $CWD/kmers/scores/hits_Abin  >  $CWD/kmers/scores/kmer_hits_Abin
awk '{print $1}' $CWD/kmers/scores/hits_Ybin  >  $CWD/kmers/scores/kmer_hits_Ybin

cat $CWD/kmers/scores/kmer_hits_Abin $CWD/kmers/scores/kmer_hits_Ybin | sort -k1b,1 | uniq -c  | awk '{print $2, $1}' >  $CWD/kmers/scores/kmer_hits_AYbin

sed '1d' $CWD/kmers/kmers_all_results > tmpfile
join -a1 -a2 -1 1 -2 1 -o '0,2.2,2.3,2.4,2.5,2.6,2.7,1.2' -e "0" $CWD/kmers/scores/kmer_hits_AYbin tmpfile > $CWD/kmers/scores/kmers_all_results_withofftargets; rm tmpfile

printf "======= generating kmers_all_results_withofftargets file =======\n"

awk -v OFS="\t" '$1=$1' $CWD/kmers/scores/kmers_all_results_withofftargets > tmpfile; mv tmpfile $CWD/kmers/scores/kmers_all_results_withofftargets

#Add column header
awk 'BEGIN {print "kmer_id\tseq\tfemale\tmale\tCQ\tsum\tbin\tofftargets"} {print}' $CWD/kmers/scores/kmers_all_results_withofftargets > tmpfile; mv tmpfile $CWD/kmers/scores/kmers_all_results_withofftargets

printf "======= done step 5 =======\n"

	