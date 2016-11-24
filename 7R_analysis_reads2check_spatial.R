#!/usr/bin/env Rscript
library (ggplot2)
#install.packages("reshape")
library(reshape2)
source("redkmer.cfg.R")
setwd(dirname(Rworkdir)) 

spatial<-read.table(paste(Rworkdir,"/read_analysis/reads2check_spatial.txt",sep=""), header=T, sep="\t")
spatial_m<-melt(spatial,id.vars=c("pacbio_read","pos"))

plotspatial<-ggplot(spatial_m)+
  geom_ribbon(aes(x=pos,ymin=0,ymax=value,fill=variable),alpha=0.5,color="black",size=0.1)+
  ylab("")+
  xlab("")+
  facet_grid(pacbio_read~.,scales="free")+
  theme_bw()+
  theme(strip.text.y=element_text(angle=0))+
  scale_fill_manual(values=c("orangered2","dodgerblue1","black"))
plot(plotspatial)
ggsave(paste(Rworkdir,"/read_analysis/reads2check_spatial.png",sep=""),width=20, height=80,limitsize = FALSE)






