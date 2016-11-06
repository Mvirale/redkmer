#!/usr/bin/env Rscript
library (ggplot2)
source("path.R")
setwd(dirname(Rworkdir)) 

pacbio<-read.table(paste(Rworkdir,"/pacBio_illmapping/pacBio_MappedReads.txt",sep=""), header=T, sep="\t")

pacbio$candidate[pacbio$CQ>=1.5]<-"X"
pacbio$candidate[pacbio$CQ<1.5]<-"A"
pacbio$candidate[pacbio$CQ<0.2]<-"Y"
pacbio$sum<-pacbio$female+pacbio$male
pacbio$mean<-(pacbio$female+pacbio$male)/2
pacbio$candidate<-as.factor(pacbio$candidate)

summary(pacbio$candidate)
summary(pacbio$CQ)

g1 <- ggplot(pacbio) + geom_point(aes(x=log10(sum), y=CQ,color=candidate),alpha=0.4)
plot(g1)
ggsave(paste(Rworkdir,"/plots/pacBIO_points.png",sep=""))

g2<- ggplot(pacbio)+geom_histogram(aes(x=CQ),binwidth = 0.05)
plot(g2)
ggsave(paste(Rworkdir,"/plots/pacBIO_histogram.png",sep=""))

