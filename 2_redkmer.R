## Export the joint_kmer_ratio_name.txt produced in the pipeline

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

pacbio<-read.table("./testproject/counts/pacBio_MappedReads.txt", header=T, sep="\t")
head(pacbio)

pacbio$candidate[pacbio$CQ>1.5]<-"X"
pacbio$candidate[pacbio$CQ<1.5]<-"A"
pacbio$candidate[pacbio$CQ<0.2]<-"Y"

pacbio$adjfactor[pacbio$CQ>1.5]<-1.3333333333
pacbio$adjfactor[pacbio$CQ<1.5]<-1
pacbio$adjfactor[pacbio$CQ<0.2]<-4

pacbio$sum<-pacbio$female+pacbio$male
pacbio$mean<-(pacbio$female+pacbio$male)/2
pacbio$candidate<-as.factor(pacbio$candidate)

pacbio$sumadj <- pacbio$sum * pacbio$adjfactor


summary(pacbio$candidate)
library (ggplot2)

ggplot(pacbio)+
  geom_point(aes(x=log10(sumadj), y=CQ,color=candidate),alpha=0.4)

ggplot(pacbio)+
  geom_point(aes(x=log10(sum), y=CQ,color=candidate),alpha=0.4)

ggplot(pacbio)+
  geom_point(aes(x=log10(mean), y=CQ,color=candidate),alpha=0.4)


summary(pacbio$CQ)
