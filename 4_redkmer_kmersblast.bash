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

awk '{print $1}' $CWD/blast/hits_allpacBio | sort | uniq >  $CWD/blast/kmer_hits_allpacBio
awk '{print $1}' $CWD/blast/hits_Xbin | sort | uniq >  $CWD/blast/kmer_hits_Xbin
awk '{print $1}' $CWD/blast/hits_Abin | sort | uniq >  $CWD/blast/kmer_hits_Abin
awk '{print $1}' $CWD/blast/hits_Ybin | sort | uniq >  $CWD/blast/kmer_hits_Ybin

cat $CWD/blast/kmer_hits_Abin $CWD/blast/kmer_hits_Ybin | sort | uniq > $CWD/blast/kmers_hits_AYbin

#cat $CWD/blast/kmer_hits_Xbin $CWD/blast/kmers_hits_AYbin | sort | uniq -c | a

# printf "======= cleaning up =======\n"

