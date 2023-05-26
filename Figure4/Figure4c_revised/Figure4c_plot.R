##################################
#source("/Users/ameira/OneDrive - Nexus365/TARGETseq_functions.R")
##################################
# GSEA plotting
##################################
library(ggplot2)
#####################################

# Read GSEA results
preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA<-read.table(file="Figure4c_GSEA_Hallmarks.txt",sep="\t",header = T)
preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$NAME<-gsub("HALLMARK_","",preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$NAME)
rownames(preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA)<-preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$NAME
preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$Group<-"CP_TP53_MPN"
preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$Group[preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$NES>0]<-"preTP53sAML"
preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$Size<-as.numeric(abs(preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$NES))
preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$Order<-length(preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$NES):1

#Re-order for plotting
preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$NAME <- factor(preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$NAME, levels = preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$NAME[order(preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA$Order)])

#####################################
ggplot(preTP53sAML_vs_CP_TP53_MPN_p53HET_GSEA,aes(x=-(log10(FDR_qval+0.001)),y=NAME,colour=Group))+
  theme_classic()+
  geom_segment(aes(x=0.55, xend=-(log10(FDR_qval+0.001)), y=NAME, yend=NAME), colour="black",size=0.5)+
  geom_point(aes(size=Size))+
  scale_color_manual(values=c("forestgreen","darkorange"))+
  scale_x_continuous(limits = c(0.5,3.5), expand = c(0, 0))

ggsave(filename="Figure4c_GSEA.png",width = 6,height = 4)
  