library (ggplot2)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
kmer<-read.table("./testproject/kmers/kmer_counts", header=T, sep="\t")

kmer$candidate[kmer$CQ>=1.5]<-"X"
kmer$candidate[kmer$CQ<1.5]<-"A"
kmer$candidate[kmer$CQ<0.2]<-"Y"
kmer$candidate<-as.factor(kmer$candidate)

summary(kmer$candidate)
summary(kmer$CQ)

g1 <- ggplot(kmer) + geom_point(aes(x=sum, y=CQ,color=candidate),alpha=0.4)
plot(g1)
ggsave("./testproject/kmers/kmer_points.png")

g2<- ggplot(kmer)+geom_histogram(aes(x=CQ),binwidth = 0.05)
plot(g2)
ggsave("./testproject/kmers/kmer_histogram.png")

