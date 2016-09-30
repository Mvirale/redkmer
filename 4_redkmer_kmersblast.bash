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

mkdir -p $CWD/blast/
mkdir -p $CWD/blast/index/

kmers=$CWD/kmers/kmer.fasta


#make blast index

$BLAST_DB -in $pacM -dbtype nucl -out $CWD/blast/index/blastdb_allpacBio
$BLAST_DB -in $CWD/bins/X_fasta -dbtype nucl -out $CWD/blast/index/blastdb_Xbin
$BLAST_DB -in $CWD/bins/A_fasta -dbtype nucl -out $CWD/blast/index/blastdb_Abin
$BLAST_DB -in $CWD/bins/Y_fasta -dbtype nucl -out $CWD/blast/index/blastdb_Ybin

#run blasts
$BLAST -db $CWD/blast/index/blastdb_allpacBio -query $kmers -out $CWD/blast/hits_allpacBio -perc_identity 100 -outfmt 6 -num_threads 2
$BLAST -db $CWD/blast/index/blastdb_Xbin -query $kmers -out $CWD/blast/hits_Xbin -perc_identity 100 -outfmt 6 -num_threads 2
$BLAST -db $CWD/blast/index/blastdb_Abin -query $kmers -out $CWD/blast/hits_Abin -perc_identity 100 -outfmt 6 -num_threads 2
$BLAST -db $CWD/blast/index/blastdb_Ybin -query $kmers -out $CWD/blast/hits_Ybin -perc_identity 100 -outfmt 6 -num_threads 2

# add this later
# -max_target_seqs 50000000

printf "======= extract blast results =======\n"

#extract blast results
awk '{print $1}' $CWD/blast/hits_allpacBio | sort | uniq >  $CWD/blast/kmer_hits_allpacBio
awk '{print $1}' $CWD/blast/hits_Xbin | sort | uniq >  $CWD/blast/kmer_hits_Xbin
awk '{print $1}' $CWD/blast/hits_Abin | sort | uniq >  $CWD/blast/kmer_hits_Abin
awk '{print $1}' $CWD/blast/hits_Ybin | sort | uniq >  $CWD/blast/kmer_hits_Ybin

cat $CWD/blast/kmer_hits_Abin $CWD/blast/kmer_hits_Ybin | sort | uniq > $CWD/blast/kmers_hits_AYbin

awk '{print $0, "1"}' $CWD/blast/kmer_hits_Xbin > $CWD/blast/kmer_hits_Xbin_lab 
awk '{print $0, "0.1"}' $CWD/blast/kmers_hits_AYbin > $CWD/blast/kmers_hits_AYbin_lab

#joining
join -a1 -a2 -1 1 -2 1 -o '0,1.2,2.2' -e "0" $CWD/blast/kmer_hits_Xbin_lab $CWD/blast/kmers_hits_AYbin_lab | awk '{print $1, ($2+$3)}' > $CWD/blast/kmers_merge

awk '{if ($2==0.1){$2="AY"};print}' $CWD/blast/kmers_merge > $CWD/blast/kmers_merge.b
awk '{if ($2==1){$2="X"};print}' $CWD/blast/kmers_merge.b > $CWD/blast/kmers_merge.c
awk '{if ($2==1.1){$2="XAY"};print}' $CWD/blast/kmers_merge.c > $CWD/blast/kmers_blast_hits


join -a1 -a2 -1 1 -2 1 -o '0,2.2,2.3,2.4,2.5,2.6,1.2' -e "nohits" $CWD/blast/kmers_blast_hits $CWD/kmers/kmers_to_merge > $CWD/kmers/kmers_all_results_a

awk -v OFS="\t" '$1=$1' $CWD/kmers/kmers_all_results_a > $CWD/kmers/kmers_all_results_b

#Add column header
awk 'BEGIN {print "kmer_id\tseq\tfemale\tmale\tCQ\tsum\tbin"} {print}' $CWD/kmers/kmers_all_results_b > $CWD/kmers/kmers_all_results

printf "======= cleaning up =======\n"

rm $CWD/blast/kmer_hits_allpacBio
rm $CWD/blast/kmer_hits_Xbin
rm $CWD/blast/kmer_hits_Abin
rm $CWD/blast/kmer_hits_Ybin
rm $CWD/blast/kmers_hits_AYbin
rm $CWD/blast/kmer_hits_Xbin_lab 
rm $CWD/blast/kmers_hits_AYbin_lab
rm $CWD/blast/kmers_merge
rm $CWD/blast/kmers_merge.b
rm $CWD/blast/kmers_merge.c
rm $CWD/kmers/kmers_to_merge
rm $CWD/kmers/kmers_all_results_a
rm $CWD/kmers/kmers_all_results_b