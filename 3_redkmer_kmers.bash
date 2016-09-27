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

mkdir -p $CWD/kmers

## Use this to skip sections if needed
#if false; then
## everything below "if false; then" is skipped
#fi
## everything below "fi" is run

printf "======= using jellyfish to create kmers of lenght 30 from each Illumina library =======\n"

# make kmer libraries

$JFISH count -C -L 2 -m 30 $illM -o $CWD/kmers/m -c 3 -s 1000000000 -t 2
$JFISH count -C -L 2 -m 30 $illF -o $CWD/kmers/f -c 3 -s 1000000000 -t 2
$JFISH dump $CWD/kmers/m -c -L 2 -o $CWD/kmers/m.counts
$JFISH dump $CWD/kmers/f -c -L 2 -o $CWD/kmers/f.counts


# sort kmers libraries

sort $CWD/kmers/m.counts > $CWD/kmers/m.sorted
sort $CWD/kmers/f.counts > $CWD/kmers/f.sorted

# Join libraries reporting all results and adding a 0 

join -a1 -a2 -1 1 -2 1 -o '0,1.2,2.2' -e "0" $CWD/kmers/f.sorted $CWD/kmers/m.sorted > $CWD/kmers/kmer_counts

# Remove kmers with no male match
awk '{if ($3>0) print}' $CWD/kmers/kmer_counts > $CWD/kmers/kmer_shared_counts_a


# Add a line with kmer_increasing_num to name the kmers
awk '{printf("%.1f %s\n", 1+(NR-1), $0)}' $CWD/kmers/kmer_shared_counts_a > $CWD/kmers/kmer_shared_counts_b
awk '{print "kmer_"$0}' $CWD/kmers/kmer_shared_counts_b > $CWD/kmers/kmers


printf "======= cleaning up =======\n"

#rm $CWD/kmers/m
#rm $CWD/kmers/f
#rm $CWD/kmers/m.counts
#rm $CWD/kmers/f.counts
#rm $CWD/kmers/kmer_counts
#rm $CWD/kmers/kmer_shared_counts_a
#rm $CWD/kmers/kmer_shared_counts_b

