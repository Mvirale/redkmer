#!/usr/bin/env Rscript
library (ggplot2)
source("path.R")
setwd(dirname(Rworkdir))

kmer<-read.table(paste(Rworkdir,"/kmers/kmer_results.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE)

kmer$candidate[kmer$CQ>=1.5]<-"X"
kmer$candidate[kmer$CQ<1.5]<-"A"
kmer$candidate[kmer$CQ<0.2]<-"Y"
kmer$candidate<-as.factor(kmer$candidate)
kmer$label<-kmer$bin
kmer$label[kmer$offtargets>0]<-"offtargets"
kmer$bin<-as.factor(kmer$bin)
kmer$label<-as.factor(kmer$label)
summary(kmer$bin)
summary(kmer$CQ)
summary(kmer$offtargets)
summary(kmer$label)


g1 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ,color=candidate),alpha=0.4)+
  scale_color_manual(values=c("springgreen4","red2","dodgerblue2"))
ggsave((paste(Rworkdir,"/plots/kmer_analysis_1.png",sep="")),width=13, height=10)

g2<- ggplot(kmer)+geom_histogram(aes(x=CQ),binwidth = 0.05)
ggsave((paste(Rworkdir,"/plots/kmer_analysis_2.png",sep="")),width=13, height=10)

g3<- ggplot(kmer)+geom_histogram(aes(x=log10(sum)),binwidth = 0.05)
ggsave((paste(Rworkdir,"/plots/kmer_analysis_3.png",sep="")),width=13, height=10)

g4<- ggplot(kmer)+geom_histogram(aes(x=sum),binwidth = 1)
ggsave((paste(Rworkdir,"/plots/kmer_analysis_4.png",sep="")),width=13, height=10)

g5 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ,color=bin),alpha=0.8)+
  scale_color_manual(values=c("grey","black","red2","grey"))+
  ylim(0,3)
ggsave((paste(Rworkdir,"/plots/kmer_analysis_5.png",sep="")),width=13, height=10)

g6 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ,color=label),alpha=0.8)+
  scale_color_manual(values=c("grey","black","blue","red","orange"))+
  ylim(0,3)
ggsave((paste(Rworkdir,"/plots/kmer_analysis_6.png",sep="")),width=13, height=10)

kmersX<-subset(kmer,kmer$bin=="X")
selectSum<-as.numeric(quantile(kmersX$sum,c(.995)))
kmer$selection<-"bad candidates"
kmer$selection[kmer$sum>selectSum & kmer$CQ>1.5 & kmer$bin=="X" & kmer$offtargets<1 ]<-"good kmers"
kmer$selection<-as.factor(kmer$selection)

g7 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ, color=selection),alpha=0.8)+
  scale_color_manual(values=c("grey","red2"))+
  ylim(0,3)
ggsave((paste(Rworkdir,"/plots/kmer_analysis_7.png",sep="")),width=13, height=10)

candidateXkmers<-subset(kmer,kmer$selection=="good kmers")
candidateXkmersSeq<-candidateXkmers[,c(1,2)]
write.table(candidateXkmers,file=(paste(Rworkdir,"/kmers/candidateXkmers.txt",sep="")),sep="\t",row.names = F,quote=F,col.names=F)
write.table(candidateXkmersSeq,file=(paste(Rworkdir,"/kmers/candidateXkmers.seq",sep="")),sep="\t",row.names = F,quote=F,col.names=F)


