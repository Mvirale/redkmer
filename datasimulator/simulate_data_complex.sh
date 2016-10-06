#!/bin/bash

PY2=/usr/bin/python2
SIMUDIR=/home/nikolai/Software/redkmer/datasimulator
DATASETID=complex

mkdir -p ../simulateddatasets/
mkdir -p ../simulateddatasets/${DATASETID}/
mkdir -p ../simulateddatasets/${DATASETID}/refgenome/
mkdir -p ../simulateddatasets/${DATASETID}/refgenome/sequences
mkdir -p ../simulateddatasets/${DATASETID}/testreadspac/
mkdir -p ../simulateddatasets/${DATASETID}/testreadsill/

chromTARGETDIR=../simulateddatasets/${DATASETID}/refgenome/
pacTARGETDIR=../simulateddatasets/${DATASETID}/testreadspac/
illTARGETDIR=../simulateddatasets/${DATASETID}/testreadsill/

#size in MB
GSIZE=0.1
MSIZE=0.01

#Generate random reference chromosomes
./readgenerators/rangen $GSIZE > A.fasta
./readgenerators/rangen $GSIZE > X.fasta
./readgenerators/rangen $GSIZE > Y.fasta
./readgenerators/rangen $MSIZE > M.fasta

#-------------------------------- chromsome repeats and seqs ----------------------------------------------

function repeatgen {
./readgenerators/rangen 0.1 > ${1}.s
tail ${1}.s -n ${2} > tmpfile; mv tmpfile ${1}.s
for ((c=0; c<=${3}; c++))
do
cat ${1}.s ${1}.s > tmpfile; mv tmpfile ${1}.s
cat head.txt ${1}.s > tmpfile; mv tmpfile ${1}.s
$PY2 mutate.py ${1}.s ${4} > tmpfile; mv tmpfile ${1}.s
sed '1d' ${1}.s > tmpfile; mv tmpfile ${1}.s
done
fold -w 80 ${1}.s > tmpfile; mv tmpfile ${1}.s
}

function seqgen {
./readgenerators/rangen ${2} > ${1}.s
sed '1d' ${1}.s > tmpfile; mv tmpfile ${1}.s
sed '$ d' ${1}.s > tmpfile; mv tmpfile ${1}.s
}

repeatgen Xrep1 2 6 0.0025
repeatgen Xrep2 2 6 0.0025
repeatgen Xrep3 2 6 0.0025
repeatgen Xrep4 2 6 0.0025
repeatgen Arep 2 6 0.0025
repeatgen Yrep 2 6 0.0025
repeatgen XYrep 2 6 0.0025
repeatgen AXrep 2 6 0.0025
repeatgen AYrep 2 6 0.0025
repeatgen ALLrep 2 6 0.0025

seqgen XYshared 0.025
seqgen XAshared 0.025
seqgen AYshared 0.025
seqgen ALLshared 0.025

chromname=">A" ; sed "1s/.*/${chromname}/" A.fasta > tmpfile; mv tmpfile A.fasta
chromname=">X"; sed "1s/.*/${chromname}/" X.fasta > tmpfile; mv tmpfile X.fasta
chromname=">Y" ; sed "1s/.*/${chromname}/" Y.fasta > tmpfile; mv tmpfile Y.fasta
chromname=">M"; sed "1s/.*/${chromname}/" M.fasta > tmpfile; mv tmpfile M.fasta

cat X.fasta Xrep1.s Xrep2.s Xrep3.s Xrep4.s XYrep.s AXrep.s ALLrep.s XYshared.s XAshared.s ALLshared.s > tmpfile; mv tmpfile X.fasta
cat A.fasta Arep.s AXrep.s AYrep.s ALLrep.s XAshared.s AYshared.s ALLshared.s > tmpfile; mv tmpfile A.fasta
cat Y.fasta Yrep.s XYrep.s AYrep.s ALLrep.s XYshared.s AYshared.s ALLshared.s > tmpfile; mv tmpfile Y.fasta

cat A.fasta A.fasta X.fasta Y.fasta M.fasta > MaleGenome.fasta
cat A.fasta A.fasta X.fasta X.fasta M.fasta M.fasta > FemaleGenome.fasta

#------------------------------- Illumina ------------------------------------------

# error rates
ILLsnp=0.002
ILLins=0.002
ILLdel=0.002
#0,001
# length
ILLsize=100
# number
ILLnum=60k

#generate reads
./readgenerators/randomreads.sh ref=MaleGenome.fasta out=${illTARGETDIR}raw_m.fastq snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=FemaleGenome.fasta out=${illTARGETDIR}raw_f.fastq snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1

#-------------------------------- Pacbio --------------------------------------------

# error rates
PACerror=0.002
# coverage
PACcov=8
# length
PACsize=3000

./readgenerators/fasta2DAM MaleGenome MaleGenome.fasta
./readgenerators/simulator MaleGenome.dam -c${PACcov}. -m${PACsize} -e${PACerror} > ${pacTARGETDIR}m_pac.fasta
rm ${SIMUDIR}/*.dam

#------------------------------- Copy reference sequences ------------------------------------------

cp *.fasta $chromTARGETDIR
cp *.s ${chromTARGETDIR}sequences

rm ${SIMUDIR}/*.fasta
rm ${SIMUDIR}/*.s



