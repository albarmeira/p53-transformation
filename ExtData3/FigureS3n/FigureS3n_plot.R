##################################
#Plotting differentially enriched pathways in IF0131 cytogenetic subclones
#Author: Alba Rodriguez-Meira, DPhil
#Date created: 14th October 2021
#Date last modified: 14th October 2021
##################################

reactome<-read.table(file="output/gsea_report_for_cloneB_REACTOME.txt",header = T,sep="\t")

reactome$Pathway<-gsub("REACTOME_","",reactome$Pathway)
reactome$logFDR<-(-log10(reactome$FDR_q_val))

#Re-order for plotting
reactome$Order<-length(reactome$NES):1
reactome$Pathway <- factor(reactome$Pathway, levels = reactome$Pathway[order(reactome$Order)])

#####################################
library(ggplot2)
ggplot(reactome,aes(x=logFDR,y=Pathway))+
  geom_segment(aes(x=2, xend=logFDR, y=Pathway, yend=Pathway), colour="black")+
  geom_point(aes(colour=Group,size=NES))+
  theme_classic()+
  scale_size(range = c(4,10))+
  scale_x_continuous(limits = c(2,6), expand = c(0, 0))+
  #xlim(values=c(2,6))+
  scale_color_manual(values=c("blue","darkorange1","darkgoldenrod1","darkorchid"))

ggsave(file="FigureS3n.png",width = 6,height = 4)
