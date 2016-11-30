#!/bin/bash
#PBS -N redkmer2
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=16:mem=16gb
#PBS -e /home/nikiwind/reports
#PBS -o /home/nikiwind/reports

if [ -z ${PBS_ENVIRONMENT+x} ]
then
echo "---> running on the Perugia numbercruncher..."
source redkmer.cfg
else
echo "---> running on HPC cluster..."
source $PBS_O_WORKDIR/redkmer.cfg
module load samtools
module load bowtie/1.1.1
fi

printf "======= calculating library sizes =======\n"

illLIBMsize=$(wc -l $illM | awk '{print ($1/4)}')
illLIBFsize=$(wc -l $illF | awk '{print ($1/4)}')
illnorm=$((($illLIBMsize+$illLIBFsize)/2))

printf "======= Building index of pacBIO reads =======\n"


if [ -z ${PBS_ENVIRONMENT+x} ]
then
$BOWTIEB $pacM $CWD/pacBio_illmapping/index/m_pac
else
cp $pacM $TMPDIR
$BOWTIEB $TMPDIR/m_pac.fasta $TMPDIR/m_pac
cp $TMPDIR/*ebwt $CWD/pacBio_illmapping/index/  2>/dev/null || :
cp $TMPDIR/*ebwtl $CWD/pacBio_illmapping/index/  2>/dev/null || :
fi


printf "======= Mapping illumina reads to pacBIO reads =======\n"

if [ -z ${PBS_ENVIRONMENT+x} ]
then 
	$BOWTIE -a -t -p $CORES -v 0 $CWD/pacBio_illmapping/index/m_pac --suppress 1,2,4,5,6,7,8,9 $illF 1> $CWD/pacBio_illmapping/mapping_rawdata/female.txt 2> $CWD/pacBio_illmapping/logs/female_log.txt
	$BOWTIE -a -t -p $CORES -v 0 $CWD/pacBio_illmapping/index/m_pac --suppress 1,2,4,5,6,7,8,9 $illM 1> $CWD/pacBio_illmapping/mapping_rawdata/male.txt 2> $CWD/pacBio_illmapping/logs/male_log.txt
else
	cp $illF $TMPDIR
	cp $illM $TMPDIR
	$BOWTIE -a -t -p $CORES -v 0 $TMPDIR/m_pac --suppress 1,2,4,5,6,7,8,9 $TMPDIR/f.fastq 1> $TMPDIR/female.txt 2> $CWD/pacBio_illmapping/logs/female_log.txt
	$BOWTIE -a -t -p $CORES -v 0 $TMPDIR/m_pac --suppress 1,2,4,5,6,7,8,9 $TMPDIR/m.fastq 1> $TMPDIR/male.txt 2> $CWD/pacBio_illmapping/logs/male_log.txt
fi

printf "======= sort and counting files =======\n"

if [ -z ${PBS_ENVIRONMENT+x} ]
then 
	split --number=l/$CORES $CWD/pacBio_illmapping/mapping_rawdata/female.txt $CWD/pacBio_illmapping/mapping_rawdata/_sorttmp;
	ls -1 $CWD/pacBio_illmapping/mapping_rawdata/_sorttmp* | (while read SORTFILE; do sort -k1b,1 $SORTFILE -o $SORTFILE & done;
	wait
	)
	sort -m $CWD/pacBio_illmapping/mapping_rawdata/_sorttmp* | uniq -c > $CWD/pacBio_illmapping/mapping_rawdata/female_uniq
	rm $CWD/pacBio_illmapping/mapping_rawdata/_sorttmp*

	split --number=l/$CORES $CWD/pacBio_illmapping/mapping_rawdata/male.txt $CWD/pacBio_illmapping/mapping_rawdata/_sorttmp;
	ls -1 $CWD/pacBio_illmapping/mapping_rawdata/_sorttmp* | (while read SORTFILE; do sort -k1b,1 $SORTFILE -o $SORTFILE & done;
	wait
	)
	sort -m $CWD/pacBio_illmapping/mapping_rawdata/_sorttmp* | uniq -c > $CWD/pacBio_illmapping/mapping_rawdata/male_uniq
	rm $CWD/pacBio_illmapping/mapping_rawdata/_sorttmp*
	
	#time sort -k1b,1 --parallel=$LESSCORES -T $CWD/temp --buffer-size=5G $CWD/pacBio_illmapping/mapping_rawdata/female.txt | uniq -c > $CWD/pacBio_illmapping/mapping_rawdata/female_uniq
	#time sort -k1b,1 --parallel=$LESSCORES -T $CWD/temp --buffer-size=5G $CWD/pacBio_illmapping/mapping_rawdata/male.txt | uniq -c > $CWD/pacBio_illmapping/mapping_rawdata/male_uniq
	
else

	split --number=l/$CORES $TMPDIR/female.txt $TMPDIR/_sorttmp;
	ls -1 $TMPDIR/_sorttmp* | (while read SORTFILE; do sort -k1b,1 -T $TMPDIR/temp $SORTFILE -o $SORTFILE & done;
	wait
	)
	sort -m $TMPDIR/_sorttmp* | uniq -c > $TMPDIR/female_uniq
	rm $TMPDIR/_sorttmp*

	split --number=l/$CORES $TMPDIR/male.txt $TMPDIR/_sorttmp;
	ls -1 $TMPDIR/_sorttmp* | (while read SORTFILE; do sort -k1b,1 -T $TMPDIR/temp $SORTFILE -o $SORTFILE & done;
	wait
	)
	sort -m $TMPDIR/_sorttmp* | uniq -c > $TMPDIR/male_uniq
	rm $CWD/pacBio_illmapping/mapping_rawdata/_sorttmp*

	#sort -k1b,1 -T $TMPDIR/temp --buffer-size=10G $TMPDIR/female.txt | uniq -c > $TMPDIR/female_uniq
	#sort -k1b,1 -T $TMPDIR/temp --buffer-size=10G $TMPDIR/male.txt | uniq -c > $TMPDIR/male_uniq

	cp $TMPDIR/female_uniq $CWD/pacBio_illmapping/mapping_rawdata/female_uniq
	cp $TMPDIR/male_uniq $CWD/pacBio_illmapping/mapping_rawdata/male_uniq
fi

printf "======= merging female and male pacBio_illmapping =======\n"

join -a1 -a2 -1 2 -2 2 -o '0,1.1,2.1' -e "0" $CWD/pacBio_illmapping/mapping_rawdata/female_uniq $CWD/pacBio_illmapping/mapping_rawdata/male_uniq > $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= normalizing to library size =======\n"
awk -v ma="$illLIBMsize" -v fema="$illLIBFsize" -v le="$illnorm" '{print $1, ($2*fema/le), ($3*ma/le)}' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= calculating CQ of pacBIO reads =======\n"

awk '{print $0, (($2+1)/($3+1))}' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= calculating sum of pacBio_illmapping on pacBIO reads =======\n"

awk '{print $0, ($2+$3)}' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= calculating LSum (Sum/length of PBreads * median PBread length  =======\n"

rm -f $pacM.fai
$SAMTOOLS faidx $pacM
awk '{print $1, $2}' $pacM.fai | sort -k1b,1 > $pacM.lengths
join -a1 -a2 -1 1 -1 1 -o'0,2.2,1.2,1.3,1.4,1.5' -e "0" $CWD/pacBio_illmapping/mapping_rawdata/merge $pacM.lengths > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

medianlength=$(awk '{print $2}' $pacM.lengths | sort -n | awk '
  BEGIN {
    c = 0;
    sum = 0;
  }
  $1 ~ /^[0-9]*(\.[0-9]*)?$/ {
    a[c++] = $1;
    sum += $1;
  }
  END {
    ave = sum / c;
    if( (c % 2) == 1 ) {
      median = a[ int(c/2) ];
    } else {
      median = ( a[c/2] + a[c/2-1] ) / 2;
    }
    OFS="\t";
    print median;
  }
')

awk -v ml="$medianlength" '{print $0, ($6 / $2 * ml)}' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= filter length and LSum (>=2000bp and LSum>=50)  =======\n"

awk -v pl="$pac_length" '{if($2>=pl)print $0}' $CWD/pacBio_illmapping/mapping_rawdata/merge | awk -v ls="$LSum" '{if ($7>=ls) print $0}' > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge 

# Replace space with tabs
awk -v OFS="\t" '$1=$1' $CWD/pacBio_illmapping/mapping_rawdata/merge > tmpfile; mv tmpfile $CWD/pacBio_illmapping/mapping_rawdata/merge

printf "======= generating pacBio_MappedReads.txt file  =======\n"

# Add column header
awk 'BEGIN {print "pacbio_read\tbp\tfemale\tmale\tCQ\tSum\tLSum"} {print}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_illmapping/pacBio_MappedReads.txt

printf "======= creating chromosomal bins of pacbio reads =======\n"

awk '{if($5>=1.5 && $5<10) print $1}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_bins/X_reads
awk '{if($5<1.5 && $5>0.2) print $1}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_bins/A_reads
awk '{if($5<0.2) print $1}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_bins/Y_reads
awk '{if($5>=10) print $1}' $CWD/pacBio_illmapping/mapping_rawdata/merge > $CWD/pacBio_bins/GA_reads

# Get sequences of pacBio bins

cat $CWD/pacBio_bins/X_reads | xargs $SAMTOOLS faidx $pacM > $CWD/pacBio_bins/fasta/Xbin.fasta
cat $CWD/pacBio_bins/A_reads | xargs $SAMTOOLS faidx $pacM > $CWD/pacBio_bins/fasta/Abin.fasta
cat $CWD/pacBio_bins/Y_reads | xargs $SAMTOOLS faidx $pacM > $CWD/pacBio_bins/fasta/Ybin.fasta
cat $CWD/pacBio_bins/GA_reads | xargs $SAMTOOLS faidx $pacM > $CWD/pacBio_bins/fasta/GAbin.fasta


printf "======= done  step 2 =======\n"

