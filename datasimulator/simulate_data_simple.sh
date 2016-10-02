#!/bin/bash

PY2=/usr/bin/python2
SIMUDIR=/home/nikolai/Software/redkmer/datasimulator
DATASETID=simple

mkdir -p ../simulateddatasets/
mkdir -p ../simulateddatasets/${DATASETID}/
mkdir -p ../simulateddatasets/${DATASETID}/refgenome/
mkdir -p ../simulateddatasets/${DATASETID}/testreadspac/
mkdir -p ../simulateddatasets/${DATASETID}/testreadsill/

chromTARGETDIR=../simulateddatasets/${DATASETID}/refgenome/
pacTARGETDIR=../simulateddatasets/${DATASETID}/testreadspac/
illTARGETDIR=../simulateddatasets/${DATASETID}/testreadsill/

#size in MB
GSIZE=0.05
MSIZE=0.01

#Generate random reference chromosomes
./readgenerators/rangen $GSIZE > A.fasta
./readgenerators/rangen $GSIZE > X.fasta
./readgenerators/rangen $GSIZE > Y.fasta
./readgenerators/rangen $MSIZE > M.fasta


#X repeats
tail X.fasta -n 2 > rep.fasta
for ((c=0; c<=6; c++))
do
cat rep.fasta rep.fasta > rep2.fasta
#cat head.text rep2.fasta > rep3.fasta
#$PY2 mutate.py rep3.fasta 0.001 > rep.fasta
mv rep2.fasta rep.fasta
#sed '1d' rep.fasta > tmpfile; mv tmpfile rep.fasta
done
#fold -w 80 rep.fasta > tmpfile; mv tmpfile rep.fasta
cat X.fasta rep.fasta > rep2.fasta
mv rep2.fasta X.fasta
mv rep.fasta Xrep.fasta
#rm rep3.fasta

chromname=">A"
sed "1s/.*/${chromname}/" A.fasta > tmpfile; mv tmpfile A.fasta
chromname=">X"
sed "1s/.*/${chromname}/" X.fasta > tmpfile; mv tmpfile X.fasta
chromname=">Y"
sed "1s/.*/${chromname}/" Y.fasta > tmpfile; mv tmpfile Y.fasta
chromname=">M"
sed "1s/.*/${chromname}/" M.fasta > tmpfile; mv tmpfile M.fasta

cat A.fasta X.fasta Y.fasta M.fasta > ${chromTARGETDIR}Genome.fasta
cp A.fasta $chromTARGETDIR
cp X.fasta $chromTARGETDIR
cp Y.fasta $chromTARGETDIR
cp M.fasta $chromTARGETDIR

#------------------------------- Illumina ------------------------------------------

# error rates
ILLsnp=0.001
ILLins=0.001
ILLdel=0.001
#0,001
# length
ILLsize=100
# number
ILLnum=8k

#generate reads
./readgenerators/randomreads.sh ref=A.fasta out=A1.illm snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=A.fasta out=A2.illm snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=X.fasta out=X.illm snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=Y.fasta out=Y.illm snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=M.fasta out=M.illm snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1

./readgenerators/randomreads.sh ref=A.fasta out=A1.illf snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=A.fasta out=A2.illf snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=X.fasta out=X1.illf snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=X.fasta out=X2.illf snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=M.fasta out=M1.illf snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=M.fasta out=M2.illf snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1

cat *.illm > raw_m.fastq
cat *.illf > raw_f.fastq

rm *.illm
rm *.illf

cp raw_m.fastq ${illTARGETDIR}
cp raw_f.fastq ${illTARGETDIR}

#-------------------------------- Pacbio --------------------------------------------

# error rates
PACerror=0.001
#0.001
# coverage
PACcov=5
# length
PACsize=3000

./readgenerators/fasta2DAM A A.fasta
./readgenerators/fasta2DAM X X.fasta
./readgenerators/fasta2DAM Y Y.fasta
./readgenerators/fasta2DAM M M.fasta

./readgenerators/simulator A.dam -c${PACcov}. -m${PACsize} -e${PACerror} > A1.pacreads
./readgenerators/simulator A.dam -c${PACcov}. -m${PACsize} -e${PACerror} > A2.pacreads
./readgenerators/simulator X.dam -c${PACcov}. -m${PACsize} -e${PACerror} > X.pacreads
./readgenerators/simulator Y.dam -c${PACcov}. -m${PACsize} -e${PACerror} > Y.pacreads
./readgenerators/simulator M.dam -c${PACcov}. -m${PACsize} -e${PACerror} > M.pacreads

cat *.pacreads* > m_pac.fasta
rm *.pacreads
rm *.dam
cp m_pac.fasta ${pacTARGETDIR}

rm ${SIMUDIR}/*.fasta
rm ${SIMUDIR}/*.fastq

