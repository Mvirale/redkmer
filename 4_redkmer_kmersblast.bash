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
mkdir -p $CWD/blast/bin_rawdata

kmers=$CWD/kmers/kmer.fasta


printf "======= Creating bowtie index for PacBio bins =======\n"

$BLAST_DB -in $pacM -dbtype nucl -out $CWD/blast/index/blastdb_allpacBio
$BLAST_DB -in $CWD/bins/X_fasta -dbtype nucl -out $CWD/blast/index/blastdb_Xbin
$BLAST_DB -in $CWD/bins/A_fasta -dbtype nucl -out $CWD/blast/index/blastdb_Abin
$BLAST_DB -in $CWD/bins/Y_fasta -dbtype nucl -out $CWD/blast/index/blastdb_Ybin

printf "======= Running blast against X chromosome bin =======\n"
$BLAST -db $CWD/blast/index/blastdb_Xbin -query $kmers -out $CWD/blast/bin_rawdata/hits_Xbin -perc_identity 100 -outfmt 6 -num_threads 2
printf "======= Running blast against autosomal bin =======\n"
$BLAST -db $CWD/blast/index/blastdb_Abin -query $kmers -out $CWD/blast/bin_rawdata/hits_Abin -perc_identity 100 -outfmt 6 -num_threads 2
printf "======= Running blast against Y chromosome bin =======\n"
$BLAST -db $CWD/blast/index/blastdb_Ybin -query $kmers -out $CWD/blast/bin_rawdata/hits_Ybin -perc_identity 100 -outfmt 6 -num_threads 2

# add this later
# -max_target_seqs 50000000

printf "======= extracting blast results =======\n"

awk '{print $1}' $CWD/blast/bin_rawdata/hits_Xbin | sort -k1b,1 | uniq >  $CWD/blast/bin_rawdata/kmer_hits_Xbin
awk '{print $1}' $CWD/blast/bin_rawdata/hits_Abin | sort -k1b,1 | uniq >  $CWD/blast/bin_rawdata/kmer_hits_Abin
awk '{print $1}' $CWD/blast/bin_rawdata/hits_Ybin | sort -k1b,1 | uniq >  $CWD/blast/bin_rawdata/kmer_hits_Ybin

cat $CWD/blast/bin_rawdata/kmer_hits_Abin $CWD/blast/bin_rawdata/kmer_hits_Ybin | sort -k1b,1 | uniq > $CWD/blast/bin_rawdata/kmers_hits_AYbin

awk '{print $0, "1"}' $CWD/blast/bin_rawdata/kmer_hits_Xbin > $CWD/blast/bin_rawdata/kmer_hits_Xbin_lab 
awk '{print $0, "0.1"}' $CWD/blast/bin_rawdata/kmers_hits_AYbin > $CWD/blast/bin_rawdata/kmers_hits_AYbin_lab

join -a1 -a2 -1 1 -2 1 -o '0,1.2,2.2' -e "0" $CWD/blast/bin_rawdata/kmer_hits_Xbin_lab $CWD/blast/bin_rawdata/kmers_hits_AYbin_lab | awk '{print $1, ($2+$3)}' > $CWD/blast/bin_rawdata/kmers_blast_hits

awk '{if ($2==0.1){$2="AY"};print}' $CWD/blast/bin_rawdata/kmers_blast_hits > tmpfile; mv tmpfile $CWD/blast/bin_rawdata/kmers_blast_hits
awk '{if ($2==1){$2="X"};print}' $CWD/blast/bin_rawdata/kmers_blast_hits > tmpfile; mv tmpfile $CWD/blast/bin_rawdata/kmers_blast_hits
awk '{if ($2==1.1){$2="XAY"};print}' $CWD/blast/bin_rawdata/kmers_blast_hits > tmpfile; mv tmpfile $CWD/blast/bin_rawdata/kmers_blast_hits

printf "======= merging blast result to kmer_counts data =======\n"

sort -k 1b,1 $CWD/blast/bin_rawdata/kmers_blast_hits > tmpfile; mv tmpfile $CWD/blast/bin_rawdata/kmers_blast_hits
sort -k 1b,1 $CWD/kmers/kmer_rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/kmers_to_merge

join -a1 -a2 -1 1 -2 1 -o '0,2.2,2.3,2.4,2.5,2.6,1.2' -e "nohits" $CWD/blast/bin_rawdata/kmers_blast_hits $CWD/kmers/kmer_rawdata/kmers_to_merge > $CWD/kmers/kmers_all_results

printf "======= generating Xkmers.fasta file for off-target analysis =======\n"

awk '{if ($7=="X") print $1, $2}' $CWD/kmers/kmers_all_results |awk '{print ">"$1"\n"$2}' > $CWD/kmers/Xkmers.fasta

printf "======= generating kmers_all_results file =======\n"

awk -v OFS="\t" '$1=$1' $CWD/kmers/kmers_all_results > tmpfile; mv tmpfile $CWD/kmers/kmers_all_results

#Add column header
awk 'BEGIN {print "kmer_id\tseq\tfemale\tmale\tCQ\tsum\tbin"} {print}' $CWD/kmers/kmers_all_results > tmpfile; mv tmpfile $CWD/kmers/kmers_all_results

printf "======= done step 4 =======\n"



