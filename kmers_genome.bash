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

mkdir -p $CWD/kmers/Refgenome_blast

genome=$CWD/refgenome/MaleGenome_FIX.fasta

printf "======= making fasta file from targetXkmers =======\n"

awk '{print ">"$1"\n"$2}' $CWD/kmers/candidateXkmers.seq > $CWD/kmers/candidateXkmers.fasta

printf "======= making blastDB from genome file =======\n"

$BLAST_DB -in $genome -dbtype nucl -out $CWD/blast/index/refGenome

printf "======= running blast of targetXkmers to genome =======\n"

#$BLAST -db $CWD/blast/index/refGenome -query $CWD/kmers/candidateXkmers.fasta -out $CWD/kmers/Refgenome_blast/candidateXkmers_vs_genome -perc_identity 90 -outfmt 6 -num_threads $CORES
#$BLAST -db $CWD/blast/index/refGenome -query $CWD/kmers/Xkmers.fasta -out $CWD/kmers/Refgenome_blast/all_Xkmers_vs_genome -perc_identity 90 -outfmt 6 -num_threads $CORES
$BLAST -db $CWD/blast/index/refGenome -query $CWD/bins/X_fasta -out $CWD/kmers/Refgenome_blast/X_pacBio_bins_vs_genome -perc_identity 90 -outfmt 6 -num_threads $CORES

#awk -v OFS="\t" '$1=$1' $CWD/kmers/Refgenome_blast/candidateXkmers_vs_genome > tmpfile; mv tmpfile $CWD/kmers/Refgenome_blast/candidateXkmers_vs_genome
#awk -v OFS="\t" '$1=$1' $CWD/kmers/Refgenome_blast/all_Xkmers_vs_genome > tmpfile; mv tmpfile $CWD/kmers/Refgenome_blast/all_Xkmers_vs_genome
awk -v OFS="\t" '$1=$1' $CWD/kmers/Refgenome_blast/X_pacBio_bins_vs_genome > tmpfile; mv tmpfile $CWD/kmers/Refgenome_blast/X_pacBio_bins_vs_genome

#Add column header
#awk 'BEGIN {print "queryid\tchromosome\tidentity\talignmentlength\tmismatches\tgapopens\tq.start\tq.end\ts.start\ts.end\tevalue\tbitscore"} {print}' $CWD/kmers/Refgenome_blast/candidateXkmers_vs_genome > tmpfile; mv tmpfile $CWD/kmers/Refgenome_blast/candidateXkmers_vs_genome
#awk 'BEGIN {print "queryid\tchromosome\tidentity\talignmentlength\tmismatches\tgapopens\tq.start\tq.end\ts.start\ts.end\tevalue\tbitscore"} {print}' $CWD/kmers/Refgenome_blast/all_Xkmers_vs_genome > tmpfile; mv tmpfile $CWD/kmers/Refgenome_blast/all_Xkmers_vs_genome
awk 'BEGIN {print "queryid\tchromosome\tidentity\talignmentlength\tmismatches\tgapopens\tq.start\tq.end\ts.start\ts.end\tevalue\tbitscore"} {print}' $CWD/kmers/Refgenome_blast/X_pacBio_bins_vs_genome > tmpfile; mv tmpfile $CWD/kmers/Refgenome_blast/X_pacBio_bins_vs_genome


$SAMTOOLS faidx $genome
printf "======= done =======\n"
	

	