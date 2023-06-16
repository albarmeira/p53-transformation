#Integration of MF MolCell 2019 samples and p53MPNAML dataset using SingleCellaR
#Data integration
#Author:Alba Rodriguez-Meira
#Date: 24th March 2021
#Modified last: 12th October 2021
####################################################

###################################
library(SingleCellExperiment)
library(Rcpp)
library(Matrix)
library(matrixStats)
library(bigmemory)
library(LinkedMatrix)
library(irlba)
library(Rtsne)
library(ggplot2)
library(gridExtra)
library(cccd)
library(ggpubr)
library(statmod)
library(umap)
library(reticulate)
library(data.table)
library(pbapply)
# library(fgsea)
# library(limma)
# library(RANN)
# library(ComplexHeatmap)
# library(circlize)
# #library(threejs)
# library(RColorBrewer)
# library(igraph)
# library(pheatmap)
# #library(AUCell)
# library(diffusionMap)
# library(destiny)
##################################
##################################
source("SingleCellaR.0.1.6/SingleCellRNASeq.R")
source("SingleCellaR.0.1.6/SingleCellClasses.R")
source("SingleCellaR.0.1.6/SingleCellGenerics.R")
source("SingleCellaR.0.1.6/SingleCellPlots.R")
source("SingleCellaR.0.1.6/SingleCellUtils.R")
source("SingleCellaR.0.1.6/SingleCell_Integration.R")
Rcpp::sourceCpp("SingleCellaR.0.1.4/src/utils.cpp")
##################################
load(file="Robjects/MFp53MPNAML_integration.rdata")
#################
colnames(MPN.integration.all@meta.data)
table(MPN.integration.all@meta.data$donor.type,MPN.integration.all@meta.data$IsPassed)

#MPN.integration.all@meta.data$louvain<-"NA"

MPN.integration.all@meta.data$merged.louvain[is.na(MPN.integration.all@meta.data$merged.louvain)] <- "NA"
table(MPN.integration.all@meta.data$merged.louvain,MPN.integration.all@meta.data$IsPassed)

MPN.integration.all@meta.data$genotype.classification[is.na(MPN.integration.all@meta.data$genotype.classification)] <- "failed_QC"

MPN.integration.all@meta.data$genotype.TP53<-MPN.integration.all@meta.data$genotype.classification
table(MPN.integration.all@meta.data$genotype.TP53,MPN.integration.all@meta.data$donor.type)
MPN.integration.all@meta.data$genotype.TP53[MPN.integration.all@meta.data$donor.type=="normal"] <- "preleukemic"
table(MPN.integration.all@meta.data$genotype.TP53,MPN.integration.all@meta.data$donor.type)
MPN.integration.all@meta.data$genotype.TP53[MPN.integration.all@meta.data$genotype.TP53=="preleukemic"] <- "TP53_WT"
table(MPN.integration.all@meta.data$genotype.TP53,MPN.integration.all@meta.data$donor.type)

table(MPN.integration.all@meta.data$genotype.collapsed,MPN.integration.all@meta.data$genotype.TP53)
table(MPN.integration.all@meta.data$donor.type,MPN.integration.all@meta.data$genotype.TP53)

MPN.integration.all@meta.data$genotype.groups<-"not_analyzed"
MPN.integration.all@meta.data$genotype.groups[MPN.integration.all@meta.data$genotype.TP53=="TP53_multihit_HOM" & MPN.integration.all@meta.data$donor.type=="AP"] <- "TP53_multihit_AP"
MPN.integration.all@meta.data$genotype.groups[MPN.integration.all@meta.data$genotype.TP53=="TP53_multihit_M2" & MPN.integration.all@meta.data$donor.type=="AP"] <- "TP53_multihit_AP"
MPN.integration.all@meta.data$genotype.groups[MPN.integration.all@meta.data$genotype.TP53=="TP53_WT" & MPN.integration.all@meta.data$donor.type=="normal"] <- "normal"
MPN.integration.all@meta.data$genotype.groups[MPN.integration.all@meta.data$genotype.TP53=="TP53_WT" & MPN.integration.all@meta.data$donor.type=="MF"] <- "MF"
MPN.integration.all@meta.data$genotype.groups[MPN.integration.all@meta.data$genotype.TP53=="TP53_WT" & MPN.integration.all@meta.data$donor.type=="AP"] <- "preLSC"
table(MPN.integration.all@meta.data$genotype.groups)
table(MPN.integration.all@meta.data$genotype.groups,MPN.integration.all@meta.data$merged.louvain)
#per LSC and erythroid clusters
MPN.integration.all@meta.data$genotype.groups[MPN.integration.all@meta.data$genotype.groups=="TP53_multihit_AP" & MPN.integration.all@meta.data$merged.louvain=="cl2_LSC"] <- "TP53_multihit_LSC"
MPN.integration.all@meta.data$genotype.groups[MPN.integration.all@meta.data$genotype.groups=="TP53_multihit_AP" & MPN.integration.all@meta.data$merged.louvain=="cl3_LSC_cycling"] <- "TP53_multihit_LSC"
MPN.integration.all@meta.data$genotype.groups[MPN.integration.all@meta.data$genotype.groups=="TP53_multihit_AP" & MPN.integration.all@meta.data$merged.louvain=="cl1_erythroid"] <- "TP53_multihit_ery"
MPN.integration.all@meta.data$genotype.groups[MPN.integration.all@meta.data$genotype.groups=="TP53_multihit_AP" & MPN.integration.all@meta.data$merged.louvain=="cl5_erythroid_cycling"] <- "TP53_multihit_ery"
MPN.integration.all@meta.data$genotype.groups[MPN.integration.all@meta.data$genotype.groups=="TP53_multihit_AP" & MPN.integration.all@meta.data$donor.type=="AP"] <- "not_analyzed"
table(MPN.integration.all@meta.data$genotype.groups)

TP53.HET<-subset(MPN.integration.all@meta.data,genotype.classification=="TP53HET")
table(TP53.HET$stage)

library(ggplot2)
subset<-subset(MPN.integration.all@meta.data,genotype.groups!="not_analyzed")
ggplot(subset,aes(x=genotype.groups,y=log1p(G2M_Core)))+geom_boxplot()
ggplot(subset,aes(x=genotype.groups,y=log1p(S_phase_Core)))+geom_boxplot()
ggplot(subset,aes(x=genotype.groups,y=UMI_count))+geom_boxplot()
ggplot(subset,aes(x=donor,y=UMI_count,fill=donor))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(subset,aes(x=donor,y=log1p(UMI_count),fill=donor))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(subset,aes(x=donor,y=log1p(detectedGenesPerCell),fill=donor))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(subset,aes(x=genotype.groups,y=detectedGenesPerCell))+geom_boxplot()


save(MPN.integration.all,file="Robjects/MFp53MPNAML_integration_20211201.rdata")
