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


printf "======= calculating library sizes =======\n"

illLIBMsize=$(wc -l $illM | awk '{print ($1/4)}')
illLIBFsize=$(wc -l $illF | awk '{print ($1/4)}')
illnorm=$((($illLIBMsize+$illLIBFsize)/2))

mkdir -p $CWD/kmers
mkdir -p $CWD/kmers/kmer_rawdata

printf "======= using jellyfish to create kmers of lenght 30 from male and female illumina libraries =======\n"

$JFISH count -C -L 2 -m 30 $illM -o $CWD/kmers/kmer_rawdata/m -c 3 -s 1000000000 -t 2
$JFISH count -C -L 2 -m 30 $illF -o $CWD/kmers/kmer_rawdata/f -c 3 -s 1000000000 -t 2
$JFISH dump $CWD/kmers/kmer_rawdata/m -c -L 2 -o $CWD/kmers/kmer_rawdata/m.counts
$JFISH dump $CWD/kmers/kmer_rawdata/f -c -L 2 -o $CWD/kmers/kmer_rawdata/f.counts

printf "======= sorting and counting kmer libraries =======\n"

sort -k1b,1 $CWD/kmers/kmer_rawdata/m.counts > $CWD/kmers/kmer_rawdata/m.sorted
sort -k1b,1 $CWD/kmers/kmer_rawdata/f.counts > $CWD/kmers/kmer_rawdata/f.sorted

printf "======= merging kmer libraries =======\n"

join -a1 -a2 -1 1 -2 1 -o '0,1.2,2.2' -e "0" $CWD/kmers/kmer_rawdata/f.sorted $CWD/kmers/kmer_rawdata/m.sorted > $CWD/kmers/kmer_rawdata/kmers_to_merge

printf "======= removing kmers absent in male library (from seq errors or low read depth) =======\n"
awk '{if ($3>0) print}' $CWD/kmers/kmer_rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/kmer_rawdata/kmers_to_merge

printf "======= creating unique ids for all kmers =======\n"
awk '{printf("%.1f %s\n", 1+(NR-1), $0)}' $CWD/kmers/kmer_rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/kmer_rawdata/kmers_to_merge
awk '{print "kmer_"$0}' $CWD/kmers/kmer_rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/kmer_rawdata/kmers_to_merge

sort -k1b,1 $CWD/kmers/kmer_rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/kmer_rawdata/kmers_to_merge

printf "======= normalizing to library sizes =======\n"
awk -v ma="$illLIBMsize" -v fema="$illLIBFsize" -v le="$illnorm" '{print $1, $2, ($3*fema/le), ($4*ma/le)}' $CWD/kmers/kmer_rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/kmer_rawdata/kmers_to_merge

printf "======= calculating kmer CQ and sum of kmer occurences in both libraries =======\n"
awk '{_div1= $4 ? ($3 / $4) : 0 ; print $0, _div1 }' $CWD/kmers/kmer_rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/kmer_rawdata/kmers_to_merge

#Adding sum
awk '{print $0, ($3+$4)}' $CWD/kmers/kmer_rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/kmer_rawdata/kmers_to_merge

#Replace space with tabs
awk -v OFS="\t" '$1=$1' $CWD/kmers/kmer_rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/kmer_rawdata/kmers_to_merge

printf "======= generating final kmer_counts file =======\n"

#Add column header
awk 'BEGIN {print "kmer_id\tseq\tfemale\tmale\tCQ\tsum"} {print}' $CWD/kmers/kmer_rawdata/kmers_to_merge > $CWD/kmers/kmer_counts

printf "======= generating fasta file for next blast =======\n"

# make fasta file from kmers for blast
awk '{print ">"$1"\n"$2}' $CWD/kmers/kmer_rawdata/kmers_to_merge > $CWD/kmers/kmer.fasta


printf "======= Done step 3=======\n"

