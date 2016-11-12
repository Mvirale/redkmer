#!/usr/bin/env Rscript
library (ggplot2)
source("path.R")
setwd(dirname(Rworkdir))

kmer<-read.table(paste(Rworkdir,"/kmers/kmer_results.txt",sep=""), header=T, sep="\t",stringsAsFactors=FALSE)

kmer$candidate[kmer$CQ>=1.5]<-"X"
kmer$candidate[kmer$CQ<1.5]<-"A"
kmer$candidate[kmer$CQ<0.2]<-"Y"
kmer$candidate[kmer$CQ>=5]<-"GA"
kmer$label<-kmer$hits_threshold
kmer$label[kmer$offtargets>0]<-"offtargets"
kmer$candidate<-as.factor(kmer$candidate)

kmer$hits_threshold<-as.factor(kmer$hits_threshold)
summary(kmer$hits_threshold)

summary(kmer$bin)
summary(kmer$CQ)
summary(kmer$offtargets)



g1 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ,color=candidate),alpha=0.4)+
  scale_color_manual(values=c("springgreen4","black","red2","dodgerblue2"))+
  theme_bw()
ggsave((paste(Rworkdir,"/plots/kmer_analysis_1.png",sep="")),width=13, height=10)

g2<- ggplot(kmer)+geom_histogram(aes(x=CQ),binwidth = 0.05)+
  theme_bw()
ggsave((paste(Rworkdir,"/plots/kmer_analysis_2.png",sep="")),width=13, height=10)

g3<- ggplot(kmer)+geom_histogram(aes(x=log10(sum)),binwidth = 0.05)+
  theme_bw()
ggsave((paste(Rworkdir,"/plots/kmer_analysis_3.png",sep="")),width=13, height=10)

g4<- ggplot(kmer)+geom_histogram(aes(x=sum),binwidth = 1)+
  theme_bw()
ggsave((paste(Rworkdir,"/plots/kmer_analysis_4.png",sep="")),width=13, height=10)

g5 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ,color=hits_threshold),alpha=0.8)+
  scale_color_manual(values=c("grey","black","red"))+
  ylim(0,3)+
  theme_bw()
ggsave((paste(Rworkdir,"/plots/kmer_analysis_5.png",sep="")),width=13, height=10)

g6 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ,color=label),alpha=0.8)+
  scale_color_manual(values=c("grey","black","blue","red"))+
  ylim(0,3)
ggsave((paste(Rworkdir,"/plots/kmer_analysis_6.png",sep="")),width=13, height=10)

kmersX<-subset(kmer,kmer$hits_threshold=="pass")
selectSum<-as.numeric(quantile(kmersX$sum,c(.995)))
kmer$selection<-"bad candidates"
kmer$selection[kmer$sum>selectSum & kmer$CQ>1.5 & kmer$hits_threshold=="pass" & kmer$offtargets<1 ]<-"good kmers"
kmer$selection<-as.factor(kmer$selection)

g7 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ, color=selection),alpha=0.8)+
  scale_color_manual(values=c("grey","red2"))+
  ylim(0,3)
ggsave((paste(Rworkdir,"/plots/kmer_analysis_7.png",sep="")),width=13, height=10)

#log10s

g8 <- ggplot(kmer) + geom_point(aes(x=log10(sum), y=CQ,color=candidate),alpha=0.4)+
  scale_color_manual(values=c("springgreen4","black","red2","dodgerblue2"))+
  theme_bw()+
  ylim(0,3)
ggsave((paste(Rworkdir,"/plots/kmer_analysis_8.png",sep="")),width=13, height=10)

g9 <-  ggplot(kmer) + geom_point(aes(x=log10(sum), y=CQ,color=hits_threshold),alpha=0.8)+
  scale_color_manual(values=c("grey","black","red"))+
  ylim(0,3)+
  theme_bw()
ggsave((paste(Rworkdir,"/plots/kmer_analysis_9.png",sep="")),width=13, height=10)

g10 <- ggplot(kmer) + geom_point(aes(x=log10(sum), y=CQ,color=label),alpha=0.8)+
  scale_color_manual(values=c("grey","black","blue","red"))+
  ylim(0,3)+
  theme_bw()
ggsave((paste(Rworkdir,"/plots/kmer_analysis_10.png",sep="")),width=13, height=10)

g11 <- ggplot(kmer) + geom_point(aes(x=log10(sum), y=CQ, color=selection),alpha=0.8)+
  scale_color_manual(values=c("grey","red2"))+
  ylim(0,3)+
  theme_bw()
ggsave((paste(Rworkdir,"/plots/kmer_analysis_11.png",sep="")),width=13, height=10)

g11 <- ggplot(kmer)+
  geom_point(aes(x=sum,y=offtargets,color=hits_threshold))+
  theme_bw()
ggsave((paste(Rworkdir,"/plots/kmer_analysis_11.png",sep="")),width=13, height=10)


candidateXkmers<-subset(kmer,kmer$selection=="good kmers")
candidateXkmersSeq<-candidateXkmers[,c(1,2)]
write.table(candidateXkmers,file=(paste(Rworkdir,"/kmers/candidateXkmers.txt",sep="")),sep="\t",row.names = F,quote=F,col.names=F)
write.table(candidateXkmersSeq,file=(paste(Rworkdir,"/kmers/candidateXkmers.seq",sep="")),sep="\t",row.names = F,quote=F,col.names=F)
