#!/bin/bash
#PBS -N redkmer_step1
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=20:mem=16gb
#PBS -e /home/nikiwind/reports
#PBS -o /home/nikiwind/reports

if [[ "$RUNINCLUSTER" -eq "1" ]]; then
source $PBS_O_WORKDIR/redkmer.cfg
module load samtools
module load jellyfish
else
source redkmer.cfg
fi

printf "======= calculating library sizes =======\n"

illLIBMsize=$(wc -l $illM | awk '{print ($1/4)}')
illLIBFsize=$(wc -l $illF | awk '{print ($1/4)}')
illnorm=$((($illLIBMsize+$illLIBFsize)/2))


printf "======= using jellyfish to create kmers of lenght 30 from male and female illumina libraries =======\n"

$JFISH count -C -L 2 -m 25 $illM -o $CWD/kmers/rawdata/m -c 3 -s 1000000000 -t $CORES
$JFISH count -C -L 2 -m 25 $illF -o $CWD/kmers/rawdata/f -c 3 -s 1000000000 -t $CORES
$JFISH dump $CWD/kmers/rawdata/m -c -L 2 -o $CWD/kmers/rawdata/m.counts
$JFISH dump $CWD/kmers/rawdata/f -c -L 2 -o $CWD/kmers/rawdata/f.counts

printf "======= sorting and counting kmer libraries =======\n"

time sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G $CWD/kmers/rawdata/m.counts > $CWD/kmers/rawdata/m.sorted &
time sort -k1b,1 --parallel=8 -T $CWD/temp --buffer-size=5G $CWD/kmers/rawdata/f.counts > $CWD/kmers/rawdata/f.sorted &

wait $(jobs -p)

printf "======= merging kmer libraries =======\n"

join -a1 -a2 -1 1 -2 1 -o '0,1.2,2.2' -e "0" $CWD/kmers/rawdata/f.sorted $CWD/kmers/rawdata/m.sorted > $CWD/kmers/rawdata/kmers_to_merge

printf "======= removing kmers absent in male library (from seq errors or low read depth) =======\n"
awk '{if ($3>0) print}' $CWD/kmers/rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_to_merge

printf "======= creating unique ids for all kmers =======\n"
awk '{printf("%.1f %s\n", 1+(NR-1), $0)}' $CWD/kmers/rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_to_merge
awk '{print "kmer_"$0}' $CWD/kmers/rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_to_merge

sort -k1b,1 $CWD/kmers/rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_to_merge

printf "======= normalizing to library sizes =======\n"
awk -v ma="$illLIBMsize" -v fema="$illLIBFsize" -v le="$illnorm" '{print $1, $2, ($3*fema/le), ($4*ma/le)}' $CWD/kmers/rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_to_merge

printf "======= calculating kmer CQ and sum of kmer occurences in both libraries =======\n"
awk '{_div1= $4 ? ($3 / $4) : 0 ; print $0, _div1 }' $CWD/kmers/rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_to_merge

#Adding sum
awk '{print $0, ($3+$4)}' $CWD/kmers/rawdata/kmers_to_merge > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmers_to_merge

#Replace space with tabs
awk -v OFS="\t" '$1=$1' $CWD/kmers/rawdata/kmers_to_merge > $CWD/kmers/rawdata/kmer_counts

printf "======= generating final kmer_counts file =======\n"

#Add column header
awk 'BEGIN {print "kmer_id\tseq\tfemale\tmale\tCQ\tsum"} {print}' $CWD/kmers/rawdata/kmer_counts > tmpfile; mv tmpfile $CWD/kmers/rawdata/kmer_counts

printf "======= generating fasta file for next blast =======\n"

# make fasta file from kmers for blast
awk '{print ">"$1"\n"$2}' $CWD/kmers/rawdata/kmers_to_merge > $CWD/kmers/fasta/allkmers.fasta

printf "======= Done step 3=======\n"

