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

mkdir -p $CWD/kmers/blast/
mkdir -p $CWD/kmers/blast/index/
mkdir -p $CWD/kmers/blast/rawdata
mkdir -p $CWD/kmers/fasta/

kmers=$CWD/kmers/fasta/allkmers.fasta



printf "======= Creating bowtie index for PacBio bins =======\n"

#$BLAST_DB -in $pacM -dbtype nucl -out $CWD/kmers/blast/index/blastdb_allpacBio
$BLAST_DB -in $CWD/pacBio_bins/fasta/Xbin.fasta -dbtype nucl -out $CWD/kmers/blast/index/blastdb_Xbin
$BLAST_DB -in $CWD/pacBio_bins/fasta/Abin.fasta -dbtype nucl -out $CWD/kmers/blast/index/blastdb_Abin
$BLAST_DB -in $CWD/pacBio_bins/fasta/Ybin.fasta -dbtype nucl -out $CWD/kmers/blast/index/blastdb_Ybin

printf "======= Running blast against X chromosome bin =======\n"
$BLAST -db $CWD/kmers/blast/index/blastdb_Xbin -query $kmers -out $CWD/kmers/blast/rawdata/hits_Xbin -perc_identity 100 -outfmt 6 -num_threads $CORES -max_target_seqs 10
printf "======= Running blast against autosomal bin =======\n"
$BLAST -db $CWD/kmers/blast/index/blastdb_Abin -query $kmers -out $CWD/kmers/blast/rawdata/hits_Abin -perc_identity 100 -outfmt 6 -num_threads $CORES -max_target_seqs 10
printf "======= Running blast against Y chromosome bin =======\n"
$BLAST -db $CWD/kmers/blast/index/blastdb_Ybin -query $kmers -out $CWD/kmers/blast/rawdata/hits_Ybin -perc_identity 100 -outfmt 6 -num_threads $CORES -max_target_seqs 10

# add this later
# -max_target_seqs 50000000

printf "======= extracting blast results =======\n"

awk '{print $1}' $CWD/kmers/blast/rawdata/hits_Xbin | sort -k1b,1 | uniq >  $CWD/kmers/blast/rawdata/kmer_hits_Xbin
awk '{print $1}' $CWD/kmers/blast/rawdata/hits_Abin | sort -k1b,1 | uniq >  $CWD/kmers/blast/rawdata/kmer_hits_Abin
awk '{print $1}' $CWD/kmers/blast/rawdata/hits_Ybin | sort -k1b,1 | uniq >  $CWD/kmers/blast/rawdata/kmer_hits_Ybin

cat $CWD/kmers/blast/rawdata/kmer_hits_Abin $CWD/kmers/blast/rawdata/kmer_hits_Ybin | sort -k1b,1 | uniq > $CWD/kmers/blast/rawdata/kmers_hits_AYbin

awk '{print $0, "1"}' $CWD/kmers/blast/rawdata/kmer_hits_Xbin > $CWD/kmers/blast/rawdata/kmer_hits_Xbin_lab 
awk '{print $0, "0.1"}' $CWD/kmers/blast/rawdata/kmers_hits_AYbin > $CWD/kmers/blast/rawdata/kmers_hits_AYbin_lab

join -a1 -a2 -1 1 -2 1 -o '0,1.2,2.2' -e "0" $CWD/kmers/blast/rawdata/kmer_hits_Xbin_lab $CWD/kmers/blast/rawdata/kmers_hits_AYbin_lab | awk '{print $1, ($2+$3)}' > $CWD/kmers/blast/rawdata/kmers_blast_hits

awk '{if ($2==0.1){$2="AY"};print}' $CWD/kmers/blast/rawdata/kmers_blast_hits > tmpfile; mv tmpfile $CWD/kmers/blast/rawdata/kmers_blast_hits
awk '{if ($2==1){$2="X"};print}' $CWD/kmers/blast/rawdata/kmers_blast_hits > tmpfile; mv tmpfile $CWD/kmers/blast/rawdata/kmers_blast_hits
awk '{if ($2==1.1){$2="XAY"};print}' $CWD/kmers/blast/rawdata/kmers_blast_hits > tmpfile; mv tmpfile $CWD/kmers/blast/rawdata/kmers_blast_hits

printf "======= merging blast result to kmer_counts data =======\n"

sort -k 1b,1 $CWD/kmers/blast/rawdata/kmers_blast_hits > tmpfile; mv tmpfile $CWD/kmers/blast/rawdata/kmers_blast_hits
sort -k 1b,1 $CWD/kmers/rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_to_merge

join -a1 -a2 -1 1 -2 1 -o '0,2.2,2.3,2.4,2.5,2.6,1.2' -e "nohits" $CWD/kmers/blast/rawdata/kmers_blast_hits $CWD/kmers/rawdata/kmers_to_merge > $CWD/kmers/rawdata/kmers_blastresults

printf "======= generating Xkmers.fasta file for off-target analysis =======\n"

awk '{if ($7=="X") print $1, $2}' $CWD/kmers/rawdata/kmers_blastresults |awk '{print ">"$1"\n"$2}' > $CWD/kmers/fasta/Xkmers.fasta

printf "======= generating kmers_all_results file =======\n"

awk -v OFS="\t" '$1=$1' $CWD/kmers/rawdata/kmers_blastresults > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_blastresults

#Add column header
awk 'BEGIN {print "kmer_id\tseq\tfemale\tmale\tCQ\tsum\tbin"} {print}' $CWD/kmers/rawdata/kmers_blastresults > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_blastresults

printf "======= done step 4 =======\n"



