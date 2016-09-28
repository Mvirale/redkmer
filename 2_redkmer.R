library (ggplot2)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
pacbio<-read.table("./testproject_XRdegen_err/counts/pacBio_MappedReads.txt", header=T, sep="\t")

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
ggsave("./testproject_XRdegen_err/counts/pacBIO_points.png")

g2<- ggplot(pacbio)+geom_histogram(aes(x=CQ),binwidth = 0.05)
plot(g2)
ggsave("./testproject_XRdegen_err/counts/pacBIO_histogram.png")

