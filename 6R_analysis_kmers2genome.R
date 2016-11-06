#!/usr/bin/env Rscript
library (ggplot2)
source("path.R")
setwd(dirname(Rworkdir)) 


kmer_to_genome<-read.table(paste(Rworkdir,"/kmers/Refgenome_blast/candidateXkmers_vs_genome",sep=""), header=T, sep="\t")
Xbin_pacM<-read.table(paste(Rworkdir,"/kmers/Refgenome_blast/X_pacBio_bins_vs_genome",sep=""), header=T, sep="\t")

head(kmer_to_genome)
head(pacM_to_genome)

coordinates<-(read.table(paste(Rworkdir,"/refgenome/MaleGenome_FIX.fasta.fai",sep=""), header=F, sep="\t"))
coordinates<-coordinates[,c(1,2)]
colnames(coordinates)<-c("chromosome","end")
coordinates$start<-0

coordinates2<-(read.table("~/Software/coordinates2.txt", header=T, sep="\t"))
coordinates2$size<-coordinates2$end-coordinates2$start

ggplot()+
  geom_rect(data=coordinates,aes(xmin=start, xmax=end, ymin=0, ymax=1),alpha=0.5,color="black")+
  geom_jitter(data=kmer_to_genome,aes(x=s.start,y=0.5),size=0.05, alpha=0.2,color="red")+
  facet_grid(chromosome~.)+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


ggplot()+
  geom_rect(data=coordinates2,aes(xmin=start, xmax=end, ymin=0, ymax=1,fill=id),alpha=0.5,color="black")+
  geom_jitter(data=subset(kmer_to_genome,kmer_to_genome$chromosome=="X"),aes(x=s.start,y=0.5),size=0.05, alpha=0.2,color="red")+
  facet_grid(chromosome~.)+
  scale_fill_manual(values=c("grey","grey","dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4"))+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
  

ggplot()+
  geom_rect(data=coordinates,aes(xmin=start, xmax=end, ymin=0, ymax=1),alpha=0.5,color="black")+
  geom_jitter(data=Xbin_pacM,aes(x=s.start,y=0.5,size=log10(alignmentlength)), alpha=0.2,color="red")+
  facet_grid(chromosome~.)+
  scale_fill_manual(values=c("grey","grey","dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4"))+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggplot()+
  geom_rect(data=coordinates2,aes(xmin=start, xmax=end, ymin=0, ymax=1,fill=id),alpha=0.5,color="black")+
  geom_jitter(data=Xbin_pacM,aes(x=s.start,y=0.5,size=alignmentlength),alpha=0.2,color="red")+
  facet_grid(chromosome~.)+
  scale_fill_manual(values=c("grey","grey","dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue2","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4","dodgerblue4"))+
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())






