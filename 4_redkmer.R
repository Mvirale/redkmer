library (ggplot2)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
kmer_blast<-read.table("./testproject_XRdegen_err/kmers/kmers_all_results", header=T, sep="\t")
head (kmer_blast)
kmer_blast$candidate[kmer_blast$CQ>=1.5]<-"X"
kmer_blast$candidate[kmer_blast$CQ<1.5]<-"A"
kmer_blast$candidate[kmer_blast$CQ<0.2]<-"Y"
kmer_blast$candidate<-as.factor(kmer_blast$candidate)

kmer_X<-subset(kmer_blast,kmer_blast$bin=="X")
dim(kmer_X)
dim(kmer_blast)
summary(kmer_blast$candidate)
summary(kmer_blast$CQ)
summary(kmer_blast$bin)

g1 <- ggplot()+ 
  geom_point(data=kmer_blast, aes(x=log10(sum), y=CQ),alpha=0.4, color="grey")+
  geom_point(data=kmer_X, aes(x=log10(sum), y=CQ),alpha=0.4, color="red")+
  ylim(0,3)
plot(g1)
ggsave("./testproject_XRdegen_err/kmers/kmer_points_postblast.png")
