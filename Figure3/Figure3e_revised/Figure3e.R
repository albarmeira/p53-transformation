#Figure3e

################################################
library(SingleCellExperiment)
library(ggplot2)
library(SingCellaR)
load(file = "../../../Sources/preleukemic.integration.revised.rdata")
################################################
object<-preleukemic.integration
table(preleukemic.integration@meta.data$genotype.classification) #730 WT, 1106 MF and 880 preLSC
######################################
#enrichment analysis using GSVA
normalized.umi<-get_normalized_umi(object)
#genes.stem<-c("BEX1","BEX2","MLLT4","CRHBP","AVP","HLF")
genes.stem<-c("MLLT4","CRHBP","AVP","HLF")
genes.HSC<-as.character(genes.stem)

exprs<-as.matrix(normalized.umi[rownames(normalized.umi) %in% genes.HSC,])

#####
MM<-Matrix::colMeans(exprs)
MM.f<-data.frame(score=MM)
colnames(MM.f)<-c("HSC_score")
MM2<-subset(MM.f,rownames(MM.f) %in% preleukemic.integration@meta.data$Cell)
#####
rownames(preleukemic.integration@meta.data)<-preleukemic.integration@meta.data$Cell
my.matrix<-merge(preleukemic.integration@meta.data,MM2,by="row.names")

preLSC<-subset(my.matrix,genotype.classification=="preleukemic")
summary_genotypes<-as.data.frame(table(preLSC$Genotype_curated))
write.table(summary_genotypes,file="summary_genotypes_preLSC.revised.txt",sep="\t")

my.matrix$order<-my.matrix$genotype.classification
my.matrix$order[my.matrix$genotype.classification=="preleukemic"]<-3
my.matrix$order[my.matrix$genotype.classification=="MF"]<-2
my.matrix$order[my.matrix$genotype.classification=="WT"]<-1

table(my.matrix$genotype.classification)

ggplot(my.matrix,aes(x=order,y=log2(HSC_score+1),colour=order))+geom_boxplot(colour="black", outlier.shape = NA,fatten = 6)+
  geom_jitter(size=0.5,alpha=0.5)+
  stat_summary(fun = mean, geom = "point", size = 5,shape=23,fill="white",colour="black")+
  theme_classic()+
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

ggsave(file="Figure3e.png",height = 4,width = 2.5)


# Assess statistical significance

WT.score<-subset(my.matrix,genotype.classification=="WT")
MF.score<-subset(my.matrix,genotype.classification=="MF")
preLSC.score<-subset(my.matrix,genotype.classification=="preleukemic")

wilcox.test(WT.score$HSC_score,MF.score$HSC_score,paired = FALSE)$p.value
#0.002732196 #0.003
wilcox.test(WT.score$HSC_score,preLSC.score$HSC_score,paired = FALSE)$p.value
#6.028342e-06
wilcox.test(MF.score$HSC_score,preLSC.score$HSC_score,paired = FALSE)$p.value
#3.901055e-15

####################################

