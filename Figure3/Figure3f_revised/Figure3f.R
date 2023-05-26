#Figure3f

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

genes.wnt<-c("ADAM17","AXIN1","AXIN2","CCND2","CSNK1E","CTNNB1","CUL1","DKK1",
        "DKK4","DLL1","DVL2","FRAT1","FZD1","FZD8","GNAI1","HDAC11",
        "HDAC2","HDAC5","HEY1","HEY2","JAG1","JAG2","KAT2A","LEF1","MAML1","MYC","NCOR2","NCSTN","NKD1",
        "NOTCH1","NOTCH4","NUMB","PPARD","PSEN2","PTCH1","RBPJ","SKP2","TCF7","TP53","WNT1","WNT5B","WNT6")

genes<-as.character(genes.wnt)

exprs<-as.matrix(normalized.umi[rownames(normalized.umi) %in% genes,])

#####
MM<-Matrix::colMeans(exprs)
MM.f<-data.frame(score=MM)
colnames(MM.f)<-c("Wnt_signalling_score")
MM2<-subset(MM.f,rownames(MM.f) %in% preleukemic.integration@meta.data$Cell)
#####
rownames(preleukemic.integration@meta.data)<-preleukemic.integration@meta.data$Cell
my.matrix<-merge(preleukemic.integration@meta.data,MM2,by="row.names")

my.matrix$order<-my.matrix$genotype.classification
my.matrix$order[my.matrix$genotype.classification=="preleukemic"]<-3
my.matrix$order[my.matrix$genotype.classification=="MF"]<-2
my.matrix$order[my.matrix$genotype.classification=="WT"]<-1

ggplot(my.matrix,aes(x=order,y=log2(Wnt_signalling_score+1),colour=order))+
  geom_boxplot(colour="black", outlier.shape = NA,fatten = 6)+
  geom_jitter(size=0.5,alpha=0.5)+
  stat_summary(fun = mean, geom = "point", size = 5,shape=23,fill="white",colour="black")+
  theme_classic()+
  ylim(values=c(0,6.5))+
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

wilcox.test(WT.score$Wnt_signalling_score,MF.score$Wnt_signalling_score,paired = FALSE)$p.value
#0.2405874
wilcox.test(WT.score$Wnt_signalling_score,preLSC.score$Wnt_signalling_score,paired = FALSE)$p.value
#7.245258e-50
wilcox.test(MF.score$Wnt_signalling_score,preLSC.score$Wnt_signalling_score,paired = FALSE)$p.value
#1.418687e-61
