#!/bin/bash

PY2=/usr/bin/python2
SIMUDIR=/home/nikolai/Software/redkmer/datasimulator

#full complexity
DATASETID=complex

mkdir -p ../simulateddatasets/
mkdir -p ../simulateddatasets/${DATASETID}/
mkdir -p ../simulateddatasets/${DATASETID}/refgenome/
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

head A.fasta -n 1 > head.text

#-------------------------------- chromsome repeats ----------------------------------------------

#X repeats
tail X.fasta -n 2 > rep.fasta
for ((c=0; c<=8; c++))
do
cat rep.fasta rep.fasta > rep2.fasta
cat head.text rep2.fasta > rep3.fasta
$PY2 mutate.py rep3.fasta 0.005 > rep.fasta
sed '1d' rep.fasta > tmpfile; mv tmpfile rep.fasta
done
fold -w 80 rep.fasta > tmpfile; mv tmpfile rep.fasta
cat X.fasta rep.fasta > rep2.fasta
mv rep2.fasta X.fasta
mv rep.fasta Xrep.fasta
rm rep3.fasta

#A repeats
tail A.fasta -n 2 > rep.fasta
for ((c=0; c<=6; c++))
do
cat rep.fasta rep.fasta > rep2.fasta
cat head.text rep2.fasta > rep3.fasta
$PY2 mutate.py rep3.fasta 0.005 > rep.fasta
sed '1d' rep.fasta > tmpfile; mv tmpfile rep.fasta
done
fold -w 80 rep.fasta > tmpfile; mv tmpfile rep.fasta
cat A.fasta rep.fasta > rep2.fasta
mv rep2.fasta A.fasta
mv rep.fasta Arep.fasta
rm rep3.fasta

#Y repeats
tail Y.fasta -n 2 > rep.fasta
for ((c=0; c<=6; c++))
do
cat rep.fasta rep.fasta > rep2.fasta
cat head.text rep2.fasta > rep3.fasta
$PY2 mutate.py rep3.fasta 0.005 > rep.fasta
sed '1d' rep.fasta > tmpfile; mv tmpfile rep.fasta
done
fold -w 80 rep.fasta > tmpfile; mv tmpfile rep.fasta
cat Y.fasta rep.fasta > rep2.fasta
mv rep2.fasta Y.fasta
mv rep.fasta Yrep.fasta
rm rep3.fasta

#AY repeats
./readgenerators/rangen $MSIZE > AYrep.fasta
tail AYrep.fasta -n 2 > rep.fasta
for ((c=0; c<=6; c++))
do
cat rep.fasta rep.fasta > rep2.fasta
cat head.text rep2.fasta > rep3.fasta
$PY2 mutate.py rep3.fasta 0.005 > rep.fasta
sed '1d' rep.fasta > tmpfile; mv tmpfile rep.fasta
done
fold -w 80 rep.fasta > tmpfile; mv tmpfile rep.fasta
cat Y.fasta rep.fasta > tmpfile; mv tmpfile Y.fasta
cat A.fasta rep.fasta > tmpfile; mv tmpfile A.fasta
mv rep.fasta AYrep.fasta
rm rep3.fasta

#AX repeats
./readgenerators/rangen $MSIZE > AXrep.fasta
tail AXrep.fasta -n 2 > rep.fasta
for ((c=0; c<=6; c++))
do
cat rep.fasta rep.fasta > rep2.fasta
cat head.text rep2.fasta > rep3.fasta
$PY2 mutate.py rep3.fasta 0.005 > rep.fasta
sed '1d' rep.fasta > tmpfile; mv tmpfile rep.fasta
done
fold -w 80 rep.fasta > tmpfile; mv tmpfile rep.fasta
cat X.fasta rep.fasta > tmpfile; mv tmpfile X.fasta
cat A.fasta rep.fasta > tmpfile; mv tmpfile A.fasta
mv rep.fasta AXrep.fasta
rm rep3.fasta

#XY repeats
./readgenerators/rangen $MSIZE > XYrep.fasta
tail XYrep.fasta -n 2 > rep.fasta
for ((c=0; c<=6; c++))
do
cat rep.fasta rep.fasta > rep2.fasta
cat head.text rep2.fasta > rep3.fasta
$PY2 mutate.py rep3.fasta 0.005 > rep.fasta
sed '1d' rep.fasta > tmpfile; mv tmpfile rep.fasta
done
fold -w 80 rep.fasta > tmpfile; mv tmpfile rep.fasta
cat X.fasta rep.fasta > tmpfile; mv tmpfile X.fasta
cat Y.fasta rep.fasta > tmpfile; mv tmpfile Y.fasta
mv rep.fasta XYrep.fasta
rm rep3.fasta

#all repeats
./readgenerators/rangen $MSIZE > ALLrep.fasta
tail ALLrep.fasta -n 2 > rep.fasta
for ((c=0; c<=6; c++))
do
cat rep.fasta rep.fasta > rep2.fasta
cat head.text rep2.fasta > rep3.fasta
$PY2 mutate.py rep3.fasta 0.005 > rep.fasta
sed '1d' rep.fasta > tmpfile; mv tmpfile rep.fasta
done
fold -w 80 rep.fasta > tmpfile; mv tmpfile rep.fasta
cat X.fasta rep.fasta > tmpfile; mv tmpfile X.fasta
cat Y.fasta rep.fasta > tmpfile; mv tmpfile Y.fasta
cat A.fasta rep.fasta > tmpfile; mv tmpfile A.fasta
rm rep3.fasta
mv rep.fasta ALLrep.fasta

#-------------------------------- chromsome shared regions ----------------------------------------------

#AX
./readgenerators/rangen 0.05 > AXshared.fasta
$PY2 mutate.py AXshared.fasta 0.005 > AXshared_A.fasta
$PY2 mutate.py AXshared.fasta 0.005 > AXshared_X.fasta

tail AXshared_A.fasta > temp.fasta; mv temp.fasta AXshared_A.fasta
fold -w 80 AXshared_A.fasta > tmpfile; mv tmpfile AXshared_A.fasta
cat A.fasta AXshared_A.fasta > tmpfile; mv tmpfile A.fasta

tail AXshared_X.fasta > temp.fasta; mv temp.fasta AXshared_X.fasta
fold -w 80 AXshared_X.fasta > tmpfile; mv tmpfile AXshared_X.fasta
cat X.fasta AXshared_X.fasta > tmpfile; mv tmpfile X.fasta

rm AXshared.fasta


#AY
./readgenerators/rangen 0.05 > AYshared.fasta
$PY2 mutate.py AYshared.fasta 0.005 > AYshared_A.fasta
$PY2 mutate.py AYshared.fasta 0.005 > AYshared_Y.fasta

tail AYshared_A.fasta > temp.fasta; mv temp.fasta AYshared_A.fasta
fold -w 80 AYshared_A.fasta > tmpfile; mv tmpfile AYshared_A.fasta
cat A.fasta AYshared_A.fasta > tmpfile; mv tmpfile A.fasta

tail AYshared_Y.fasta > temp.fasta; mv temp.fasta AYshared_Y.fasta
fold -w 80 AYshared_Y.fasta > tmpfile; mv tmpfile AYshared_Y.fasta
cat Y.fasta AYshared_Y.fasta > tmpfile; mv tmpfile Y.fasta

rm AYshared.fasta

#XY
./readgenerators/rangen 0.05 > XYshared.fasta
$PY2 mutate.py XYshared.fasta 0.005 > XYshared_X.fasta
$PY2 mutate.py XYshared.fasta 0.005 > XYshared_Y.fasta

tail XYshared_X.fasta > temp.fasta; mv temp.fasta XYshared_X.fasta
fold -w 80 XYshared_X.fasta > tmpfile; mv tmpfile XYshared_X.fasta
cat X.fasta XYshared_X.fasta > tmpfile; mv tmpfile X.fasta

tail XYshared_Y.fasta > temp.fasta; mv temp.fasta XYshared_Y.fasta
fold -w 80 XYshared_Y.fasta > tmpfile; mv tmpfile XYshared_Y.fasta
cat Y.fasta XYshared_Y.fasta > tmpfile; mv tmpfile Y.fasta

rm XYshared.fasta

#ALL
./readgenerators/rangen 0.05 > ALLshared.fasta
$PY2 mutate.py ALLshared.fasta 0.005 > ALLshared_X.fasta
$PY2 mutate.py ALLshared.fasta 0.005 > ALLshared_Y.fasta
$PY2 mutate.py ALLshared.fasta 0.005 > ALLshared_A.fasta

tail ALLshared_X.fasta > temp.fasta; mv temp.fasta ALLshared_X.fasta
fold -w 80 ALLshared_X.fasta > tmpfile; mv tmpfile ALLshared_X.fasta
cat X.fasta ALLshared_X.fasta > tmpfile; mv tmpfile X.fasta

tail ALLshared_Y.fasta > temp.fasta; mv temp.fasta ALLshared_Y.fasta
fold -w 80 ALLshared_Y.fasta > tmpfile; mv tmpfile ALLshared_Y.fasta
cat Y.fasta ALLshared_Y.fasta > tmpfile; mv tmpfile Y.fasta

tail ALLshared_A.fasta > temp.fasta; mv temp.fasta ALLshared_A.fasta
fold -w 80 ALLshared_A.fasta > tmpfile; mv tmpfile ALLshared_A.fasta
cat A.fasta ALLshared_A.fasta > tmpfile; mv tmpfile A.fasta

rm ALLshared.fasta



#------------------------------- Copy reference sequences ------------------------------------------

cp *.fasta $chromTARGETDIR

#------------------------------- Illumina ------------------------------------------

# error rates
ILLsnp=0.0025
ILLins=0.0025
ILLdel=0.0025
#0,001
# length
ILLsize=100
# number
ILLnum=15k

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
PACerror=0.0025
#0.001
# coverage
PACcov=8
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

