#!/usr/bin/env Rscript
library (ggplot2)
setwd(dirname(sys.frame(1)$ofile)) 
source("path.R")

kmer_blast<-read.table(paste(Rworkdir,"/kmers/kmers_all_results",sep=""), header=T, sep="\t")
head (kmer_blast)
kmer_blast$candidate[kmer_blast$CQ>=1.5]<-"X"
kmer_blast$candidate[kmer_blast$CQ<1.5]<-"A"
kmer_blast$candidate[kmer_blast$CQ<0.2]<-"Y"
kmer_blast$candidate<-as.factor(kmer_blast$candidate)

print("Perform selection")

selectSum<-as.numeric(quantile(kmer_blast$sum,c(.995)))
kmer_blast$selection<-"not selected"
kmer_blast$selection[kmer_blast$bin=="X"]<-"X-kmers"
kmer_blast$selection[kmer_blast$sum>selectSum & kmer_blast$CQ>1.5]<-"candidate kmers"
kmer_blast$selection<-as.factor(kmer_blast$selection)
summary(kmer_blast$selection)
candidateXkmers<-subset(kmer_blast,kmer_blast$selection=="candidate kmers")

print("exporting data")

write.table(candidateXkmers,file=(paste(Rworkdir,"/kmers/candidateXkmers.table",sep="")),sep="\t",row.names=F,quote=F)

candidates<-candidateXkmers[,c(1,2)]
write.table(candidates,file=(paste(Rworkdir,"/kmers/candidateXkmers.seq",sep="")),sep="\t",row.names = F,quote=F,col.names=F)

print("making plot")

g1 <- ggplot()+ 
  geom_point(data=kmer_blast, aes(x=log10(sum), y=CQ, color=selection))+
  ylim(0,3)+
  scale_color_manual(values=c("springgreen4", "red2","dodgerblue2"))
plot(g1)
ggsave((paste(Rworkdir,"/kmers/kmer_points_postblast.png",sep="")),width=13, height=10)


