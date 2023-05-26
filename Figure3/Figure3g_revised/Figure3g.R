#Figure3g

################################################
library(SingleCellExperiment)
library(ggplot2)
library(SingCellaR)
load(file = "../../../Sources/preleukemic.integration.revised.rdata")
################################################
object<-preleukemic.integration
######################################
#enrichment analysis GSVA
# 
normalized.umi<-get_normalized_umi(object)

genes.cellcycle<-c("ASPM","CENPE","CENPF","DLGAP5","MKI67","NUSAP1","STMN1","TOP2A","TUBB")
genes<-as.character(genes.cellcycle)

exprs<-as.matrix(normalized.umi[rownames(normalized.umi) %in% genes,])

#####
MM<-Matrix::colMeans(exprs)
MM.f<-data.frame(score=MM)
colnames(MM.f)<-c("cell_cycle_score")
MM2<-subset(MM.f,rownames(MM.f) %in% preleukemic.integration@meta.data$Cell)
#####
rownames(preleukemic.integration@meta.data)<-preleukemic.integration@meta.data$Cell
my.matrix<-merge(preleukemic.integration@meta.data,MM2,by="row.names")

my.matrix$order<-my.matrix$genotype.classification
my.matrix$order[my.matrix$genotype.classification=="preleukemic"]<-3
my.matrix$order[my.matrix$genotype.classification=="MF"]<-2
my.matrix$order[my.matrix$genotype.classification=="WT"]<-1

ggplot(my.matrix,aes(x=order,y=log2(cell_cycle_score+1),colour=order))+
  geom_boxplot(colour="black", outlier.shape = NA,fatten = 6)+
  geom_jitter(size=0.5,alpha=0.5)+
  stat_summary(fun = mean, geom = "point", size = 5,shape=23,fill="white",colour="black")+
  theme_classic()+
  #ylim(values=c(0,6.5))+
  scale_color_manual(values=c("grey60","chartreuse3","blue"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank())+
  theme(legend.position = "none")

ggsave(file="Figure3g.png",height = 4,width = 2.5)


# Assess statistical significance

WT.score<-subset(my.matrix,genotype.classification=="WT")
MF.score<-subset(my.matrix,genotype.classification=="MF")
preLSC.score<-subset(my.matrix,genotype.classification=="preleukemic")

wilcox.test(WT.score$cell_cycle_score,MF.score$cell_cycle_score,paired = FALSE)$p.value
#0.002569984
wilcox.test(WT.score$cell_cycle_score,preLSC.score$cell_cycle_score,paired = FALSE)$p.value
#1.375586e-20
wilcox.test(MF.score$cell_cycle_score,preLSC.score$cell_cycle_score,paired = FALSE)$p.value
#1.107826e-07
