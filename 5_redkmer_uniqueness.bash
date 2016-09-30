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

$BLAST -db $CWD/blast/index/blastdb_Abin -query $kmers -out $CWD/kmers/scores/hits_Abin -perc_identity 90 -outfmt 7 -num_threads 2
$BLAST -db $CWD/blast/index/blastdb_Ybin -query $kmers -out $CWD/kmers/scores/hits_Ybin -perc_identity 90 -outfmt 6 -num_threads 2

awk '{print $1}' $CWD/kmers/scores/hits_Abin | sort | uniq >  $CWD/kmers/scores/kmer_hits_Abin
awk '{print $1}' $CWD/kmers/scores/hits_Ybin | sort | uniq >  $CWD/kmers/scores/kmer_hits_Ybin

cat $CWD/kmers/scores/kmer_hits_Abin $CWD/kmers/scores/kmer_hits_Ybin | sort | uniq >  $CWD/kmers/scores/kmer_hits_AYbin
	