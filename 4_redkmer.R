library (ggplot2)
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
kmer_blast<-read.table("./testproject_XRdegen_err/kmers/kmers_all_results", header=T, sep="\t")
head (kmer_blast)
kmer_blast$candidate[kmer_blast$CQ>=1.5]<-"X"
kmer_blast$candidate[kmer_blast$CQ<1.5]<-"A"
kmer_blast$candidate[kmer_blast$CQ<0.2]<-"Y"
kmer_blast$candidate<-as.factor(kmer_blast$candidate)

print("Perform selection")

kmer_blast$selection<-"not selected"
kmer_blast$selection[kmer_blast$bin=="X"]<-"X-kmers"
kmer_blast$selection[kmer_blast$sum>selectSum & kmer_blast$CQ]<-"candidate kmers"
kmer_blast$selection<-as.factor(kmer_blast$selection)
summary(kmer_blast$selection)
candidateXkmers<-subset(kmer_blast,kmer_blast$selection=="candidate kmers")

print("exporting data")

write.table(candidateXkmers,file= "./testproject_XRdegen_err/kmers/candidateXkmers.table",sep="\t",row.names = F,quote=F)
candidates<-candidateXkmers[,c(1,2)]
write.table(candidates,file= "./testproject_XRdegen_err/kmers/candidateXkmers.seq",sep="\t",row.names = F,quote=F,col.names=F)

print("making plot")

g1 <- ggplot()+ 
  geom_point(data=kmer_blast, aes(x=log10(sum), y=CQ, color=selection))+
  ylim(0,3)+
  scale_color_manual(values=c("springgreen4", "red2","dodgerblue2"))
plot(g1)
ggsave("./testproject_XRdegen_err/kmers/kmer_points_postblast.png",width=13, height=10)


