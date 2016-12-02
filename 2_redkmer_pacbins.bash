#!/bin/bash
#PBS -N redkmer2
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=24:mem=16gb:tmpspace=400gb
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

	printf "======= sort and counting files =======\n"
	
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
	
else

cat > ${CWD}/qsubscripts/malepacbins.bash <<EOF
#!/bin/bash
#PBS -N redkmer_m_pacb
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=24:mem=32gb:tmpspace=500gb
#PBS -e /home/nikiwind/reports
#PBS -o /home/nikiwind/reports
module load bowtie/1.1.1

	echo "==================================== Working on male pacbins ======================================="
		cp $illM XXXXX
		$BOWTIE -a -t -p $CORES -v 0 $CWD/pacBio_illmapping/index/m_pac --suppress 1,2,4,5,6,7,8,9 XXXXX/m.fastq 1> XXXXX/male.txt 2> $CWD/pacBio_illmapping/logs/male_log.txt

	echo "==================================== Done male pacbins, sorting ===================================="

		split --number=l/$CORES XXXXX/male.txt XXXXX/_sorttmp;
		ls -1 XXXXX/_sorttmp* | (while read SORTFILE; do sort -k1b,1 -T XXXXX/temp YYYYY -o YYYYY & done;
		wait
		)
		sort -m XXXXX/_sorttmp* | uniq -c > XXXXX/male_uniq
		cp XXXXX/male_uniq $CWD/pacBio_illmapping/mapping_rawdata/male_uniq

	echo "==================================== Done sorting ! ===================================="
EOF

sed 's/XXXXX/$TMPDIR/g' ${CWD}/qsubscripts/malepacbins.bash > ${CWD}/qsubscripts/malepacbins.bashX
sed 's/YYYYY/$SORTFILE/g' ${CWD}/qsubscripts/malepacbins.bashX > ${CWD}/qsubscripts/malepacbins.bash


cat > ${CWD}/qsubscripts/femalepacbins.bash <<EOF
#!/bin/bash
#PBS -N redkmer_f_pacb
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=24:mem=32gb:tmpspace=500gb
#PBS -e /home/nikiwind/reports
#PBS -o /home/nikiwind/reports
module load bowtie/1.1.1

	echo "==================================== Working on female pacbins ======================================="
		cp $illF XXXXX
		$BOWTIE -a -t -p $CORES -v 0 $CWD/pacBio_illmapping/index/m_pac --suppress 1,2,4,5,6,7,8,9 XXXXX/f.fastq 1> XXXXX/female.txt 2> $CWD/pacBio_illmapping/logs/female_log.txt

	echo "==================================== Done female pacbins, sorting ===================================="

		split --number=l/$CORES XXXXX/female.txt XXXXX/_sorttmp;
		ls -1 XXXXX/_sorttmp* | (while read SORTFILE; do sort -k1b,1 -T XXXXX/temp YYYYY -o YYYYY & done;
		wait
		)
		sort -m XXXXX/_sorttmp* | uniq -c > XXXXX/female_uniq
		cp XXXXX/female_uniq $CWD/pacBio_illmapping/mapping_rawdata/female_uniq

	echo "==================================== Done sorting ! ===================================="

EOF

sed 's/XXXXX/$TMPDIR/g' ${CWD}/qsubscripts/femalepacbins.bash > ${CWD}/qsubscripts/femalepacbins.bashX
sed 's/YYYYY/$SORTFILE/g' ${CWD}/qsubscripts/femalepacbins.bashX > ${CWD}/qsubscripts/femalepacbins.bash


	MALEJOB=$(qsub ${CWD}/qsubscripts/malepacbins.bash)
	echo $MALEJOB
	FEMALEJOB=$(qsub ${CWD}/qsubscripts/femalepacbins.bash)
	echo $FEMALEJOB	
	
	while qstat $FEMALEJOB &> /dev/null; do
	    sleep 10;
	done;

	while qstat $MALEJOB &> /dev/null; do
	    sleep 10;
	done;

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

