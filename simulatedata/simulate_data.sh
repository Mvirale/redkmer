#!/bin/bash

PY2=/usr/bin/python2
SIMUDIR=/home/nikolai/Software/redkmer/simulatedata

mkdir -p ../testprojectNEW/
mkdir -p ../testprojectNEW/
mkdir -p ../testprojectNEW/testreadspac/
mkdir -p ../testprojectNEW/testreadsill/

chromTARGETDIR=../testprojectNEW/
pacTARGETDIR=../testprojectNEW/testreadspac/
illTARGETDIR=../testprojectNEW/testreadsill/

#size in MB
GSIZE=0.08
MSIZE=0.01

#Generate random reference chromosomes
./readgenerators/rangen $GSIZE > A.fasta
./readgenerators/rangen $GSIZE > X.fasta
./readgenerators/rangen $GSIZE > Y.fasta
#./readgenerators/rangen $MSIZE > M.fasta


#XR=$(((RANDOM % 10) + 1))
#XS=$(((RANDOM % 3) + 1))
XR=8
XS=2
REPEATDEGEN=0.001

#introduce X repeats
tail X.fasta -n $XS > rep.fasta
for ((c=0; c<=XR; c++))
do
cat rep.fasta rep.fasta > rep2.fasta

head X.fasta -n 1 > head.text
cat head.text rep2.fasta > rep3.fasta
$PY2 mutate.py rep3.fasta $REPEATDEGEN > rep.fasta

#mv rep2.fasta rep.fasta

sed '1d' rep.fasta > tmpfile; mv tmpfile rep.fasta
done
fold -w 80 rep.fasta > tmpfile; mv tmpfile rep.fasta
cat X.fasta rep.fasta > rep2.fasta

mv rep2.fasta X.fasta
mv rep.fasta Xrep.fasta
rm head.text
rm rep3.fasta

cp A.fasta $chromTARGETDIR
cp X.fasta $chromTARGETDIR
cp Y.fasta $chromTARGETDIR
cp Xrep.fasta $chromTARGETDIR
#cp M.fasta $chromTARGETDIR



#------------------------------- Illumina ------------------------------------------

# error rates
ILLsnp=0
ILLins=0
ILLdel=0
#0,001
# length
ILLsize=100
# number
ILLnum=10k

#generate reads
./readgenerators/randomreads.sh ref=A.fasta out=A1.illm snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=A.fasta out=A2.illm snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=X.fasta out=X.illm snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=Y.fasta out=Y.illm snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
#./readgenerators/randomreads.sh ref=M.fasta out=M.illm snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1

./readgenerators/randomreads.sh ref=A.fasta out=A1.illf snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=A.fasta out=A2.illf snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=X.fasta out=X1.illf snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
./readgenerators/randomreads.sh ref=X.fasta out=X2.illf snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
#./readgenerators/randomreads.sh ref=M.fasta out=M1.illf snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1
#./readgenerators/randomreads.sh ref=M.fasta out=M2.illf snprate=$ILLsnp insrate=$ILLins delrate=$ILLdel reads=$ILLnum length=$ILLsize gaussian seed=-1

cat *.illm > m.fastq
cat *.illf > f.fastq

rm *.illm
rm *.illf

cp m.fastq ${illTARGETDIR}
cp f.fastq ${illTARGETDIR}

#-------------------------------- Pacbio --------------------------------------------

# error rates
PACerror=0
#0.001
# coverage
PACcov=5
# length
PACsize=3000

./readgenerators/fasta2DAM A A.fasta
./readgenerators/fasta2DAM X X.fasta
./readgenerators/fasta2DAM Y Y.fasta
#/readgenerators/fasta2DAM M M.fasta

./readgenerators/simulator A.dam -c${PACcov}. -m${PACsize} -e${PACerror} > A1.pacreads
./readgenerators/simulator A.dam -c${PACcov}. -m${PACsize} -e${PACerror} > A2.pacreads
./readgenerators/simulator X.dam -c${PACcov}. -m${PACsize} -e${PACerror} > X.pacreads
./readgenerators/simulator Y.dam -c${PACcov}. -m${PACsize} -e${PACerror} > Y.pacreads
#./readgenerators/simulator M.dam -c${PACcov}. -m${PACsize} -e${PACerror} > M.pacreads

cat *.pacreads* > m_pac.fasta
rm *.pacreads
rm *.dam
cp m_pac.fasta ${pacTARGETDIR}

