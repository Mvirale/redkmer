#!/usr/bin/env Rscript
library (ggplot2)
setwd(dirname(sys.frame(1)$ofile)) 
source("path.R")

kmer<-read.table(paste(Rworkdir,"/kmers/scores/kmers_all_results_withofftargets",sep=""), header=T, sep="\t")

kmer$candidate[kmer$CQ>=1.5]<-"X"
kmer$candidate[kmer$CQ<1.5]<-"A"
kmer$candidate[kmer$CQ<0.2]<-"Y"
kmer$candidate<-as.factor(kmer$candidate)

summary(kmer$bin)
summary(kmer$CQ)

g1 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ,color=candidate),alpha=0.4)+
  scale_color_manual(values=c("springgreen4","red2","dodgerblue2"))
ggsave((paste(Rworkdir,"/kmers/kmer_analysis_1.png",sep="")),width=13, height=10)

g2<- ggplot(kmer)+geom_histogram(aes(x=CQ),binwidth = 0.05)
ggsave((paste(Rworkdir,"/kmers/kmer_analysis_2.png",sep="")),width=13, height=10)

g3<- ggplot(kmer)+geom_histogram(aes(x=log10(sum)),binwidth = 0.05)
ggsave((paste(Rworkdir,"/kmers/kmer_analysis_3.png",sep="")),width=13, height=10)

g4<- ggplot(kmer)+geom_histogram(aes(x=sum),binwidth = 1)
ggsave((paste(Rworkdir,"/kmers/kmer_analysis_4.png",sep="")),width=13, height=10)

g5 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ,color=bin),alpha=0.8)+
  scale_color_manual(values=c("grey","black","red2","grey"))+
  ylim(0,3)
ggsave((paste(Rworkdir,"/kmers/kmer_analysis_5.png",sep="")),width=13, height=10)

g6 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ, size=log10(offtargets+1),color=bin),alpha=0.8)+
  scale_color_manual(values=c("grey","black","red2","grey"))+
  ylim(0,3)
ggsave((paste(Rworkdir,"/kmers/kmer_analysis_6.png",sep="")),width=13, height=10)


selectSum<-as.numeric(quantile(kmer$sum,c(.999)))
kmer$selection<-"bad candidates"
kmer$selection[kmer$sum>selectSum & kmer$CQ>1.5 & kmer$bin=="X" & kmer$offtargets<1 ]<-"good kmers"
kmer$selection<-as.factor(kmer$selection)

g7 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ, color=selection),alpha=0.8)+
  scale_color_manual(values=c("grey","red2"))+
  ylim(0,3)
ggsave((paste(Rworkdir,"/kmers/kmer_analysis_7.png",sep="")),width=13, height=10)

candidateXkmers<-subset(kmer,kmer$selection=="good kmers")
candidateXkmersSeq<-candidateXkmers[,c(1,2)]
write.table(candidateXkmers,file=(paste(Rworkdir,"/kmers/candidateXkmers.txt",sep="")),sep="\t",row.names = F,quote=F,col.names=F)
write.table(candidateXkmersSeq,file=(paste(Rworkdir,"/kmers/candidateXkmers.seq",sep="")),sep="\t",row.names = F,quote=F,col.names=F)


