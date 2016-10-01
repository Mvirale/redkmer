#!/usr/bin/env Rscript
library (ggplot2)
setwd(dirname(sys.frame(1)$ofile)) 
source("path.R")
kmer<-read.table(paste(Rworkdir,"/kmers/kmer_counts",sep=""), header=T, sep="\t")

kmer$candidate[kmer$CQ>=1.5]<-"X"
kmer$candidate[kmer$CQ<1.5]<-"A"
kmer$candidate[kmer$CQ<0.2]<-"Y"
kmer$candidate<-as.factor(kmer$candidate)

summary(kmer$candidate)
summary(kmer$CQ)

g1 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ,color=candidate),alpha=0.4)
plot(g1)
ggsave(paste(Rworkdir,"/kmers/kmer_points.png",sep=""))

g2<- ggplot(kmer)+geom_histogram(aes(x=CQ),binwidth = 0.05)
plot(g2)
ggsave(paste(Rworkdir,"/kmers/kmer_historgram.png",sep=""))
