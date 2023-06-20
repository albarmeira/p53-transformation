#Integration of MF MolCell 2019 samples and p53MPNAML dataset using SingleCellaR
#Data integration
#Author:Alba Rodriguez-Meira and Sean Wen
#Date: 24th March 2021
#Modified last: 1st March 2023
####################################################
library(SingCellaR)
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
library(fgsea)
library(limma)
library(RANN)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(igraph)
library(pheatmap)
library(diffusionMap)
#library(destiny)
##################################
source("../../../Sources/SingleCell_Integration_Edit.R") # Disable requirement that object should be SingCellaR

setClass("SingCellaR_Int",representation(dir_path_SingCellR_object_files="character",SingCellR_object_files="character"))

MPN.integration.all <- new("SingCellaR_Int",
                           dir_path_SingCellR_object_files=paste(path.dropbox, "Robjects/", sep=""),
                           SingCellR_object_files=c("HTMPNAML_final.revised.rdata","MF_MolCell2019.revised.rdata"))

###################################
preprocess_integration(MPN.integration.all)
####################################
filter_cells_and_genes(MPN.integration.all,min_UMIs=1000,max_UMIs=10000000,
                       min_detected_genes=500,max_detected_genes=10000,
                       max_percent_mito=80,
                       isRemovedDoublets=FALSE)
#"0/21955 cells will be filtered out from the downstream analyses!."
####################################
normalize_UMIs(MPN.integration.all,use.scaled.factor = F)
#####################################
#Add metadata 
load(file="Robjects/HTMPNAML_final.revised.rdata")
#colnames(HTMPNAML@meta.data)
#head(HTMPNAML@meta.data)
metadata_p53MPNAML<-HTMPNAML@meta.data
rownames(metadata_p53MPNAML)<-metadata_p53MPNAML$Cell
rm(HTMPNAML)

load(file="Robjects/MF_MolCell2019.revised.rdata")
metadata_MF<-MF_MolCell@meta.data
rm(MF_MolCell)

metadata_p53MPNAML_MF<-rbind(metadata_p53MPNAML,metadata_MF)
metadata_p53MPNAML_MF<-metadata_p53MPNAML_MF[MPN.integration.all@meta.data$Cell,]

MPN.integration.all@meta.data<-merge(MPN.integration.all@meta.data,metadata_p53MPNAML_MF,by.x="Cell",by.y="Cell")

MPN.integration.all@meta.data$sampleID.y<-NULL
MPN.integration.all@meta.data$UMI_count.y<-NULL
MPN.integration.all@meta.data$detectedGenesPerCell.y<-NULL
MPN.integration.all@meta.data$percent_mito.y<-NULL
MPN.integration.all@meta.data$IsPassed.y<-NULL
MPN.integration.all@meta.data$data_set[MPN.integration.all@meta.data$data_set==1]<-"p53MPNAML"
MPN.integration.all@meta.data$data_set[MPN.integration.all@meta.data$data_set==2]<-"MF"
names(MPN.integration.all@meta.data)[names(MPN.integration.all@meta.data) == 'sampleID.x'] <- 'sampleID'
names(MPN.integration.all@meta.data)[names(MPN.integration.all@meta.data) == 'UMI_count.x'] <- 'UMI_count'
names(MPN.integration.all@meta.data)[names(MPN.integration.all@meta.data) == 'detectedGenesPerCell.x'] <- 'detectedGenesPerCell'
names(MPN.integration.all@meta.data)[names(MPN.integration.all@meta.data) == 'percent_mito.x'] <- 'percent_mito'
names(MPN.integration.all@meta.data)[names(MPN.integration.all@meta.data) == 'IsPassed.x'] <- 'IsPassed'

#Classify donors per donor type
AP<-c("GH001_003","GR003_AP","IF0392","GST010","JB4211","GR004_AP","IF0131","IF0393","GR001","GR005_AP","GR007_AP","GR006_AP","IF0391","GR002","IF0318","SB5702")
CP_preTP53sAML<-c("GH001_001","GH001_002","GR003_CP","GR005_CP","GR006_CP","GR007_CP")
CP_MPN<-c("IF9118_05HS31","IF9118_0GX632","IF9131","IF9180","IF9038")
normal<-c("HD15_8650","IF0704","IF0902","IF0903","IF0905","IF0907","IF0908")

MPN.integration.all@meta.data$donor.type<-MPN.integration.all@meta.data$stage
MPN.integration.all@meta.data$donor.type[MPN.integration.all@meta.data$stage %in% AP]<-"AP"
MPN.integration.all@meta.data$donor.type[MPN.integration.all@meta.data$stage %in% CP_preTP53sAML]<-"CP_preTP53sAML"
MPN.integration.all@meta.data$donor.type[MPN.integration.all@meta.data$stage %in% CP_MPN]<-"CP_MPN"
MPN.integration.all@meta.data$donor.type[MPN.integration.all@meta.data$stage %in% normal]<-"normal"

#Get summary of donor types for manuscript
donor.type_cells<-as.data.frame(table(MPN.integration.all@meta.data$donor.type,MPN.integration.all@meta.data$IsPassed))
donor.type_cells$Var2<-NULL
colnames(donor.type_cells)<-c("donor.type","n.cells")
write.table(donor.type_cells,file="output/donors.type_cells.counts.revised.txt",row.names = F)

donors_cells<-as.data.frame(table(MPN.integration.all@meta.data$stage,MPN.integration.all@meta.data$IsPassed))
donors_cells$Var2<-NULL
colnames(donors_cells)<-c("donor.stage","n.cells")
write.table(donors_cells,file="output/donors_cells.counts.revised.txt",row.names = F)


#Remove CP_preTP53sAML and CP_MPN donors
MPN.integration.all@meta.data$IsPassed[MPN.integration.all@meta.data$stage %in% CP_preTP53sAML]<-FALSE
MPN.integration.all@meta.data$IsPassed[MPN.integration.all@meta.data$stage %in% CP_MPN]<-FALSE

#Classify genotypes
#table(MPN.integration.all@meta.data$donor.type,MPN.integration.all@meta.data$genotype.classification)

MPN.integration.all@meta.data$genotype.collapsed[MPN.integration.all@meta.data$genotype.classification=="TP53_multihit_HOM"]<-"TP53_multihit"
MPN.integration.all@meta.data$genotype.collapsed[MPN.integration.all@meta.data$genotype.classification=="TP53_multihit_M2"]<-"TP53_multihit"
MPN.integration.all@meta.data$genotype.collapsed[MPN.integration.all@meta.data$genotype.classification=="TP53HET"]<-"TP53HET"
MPN.integration.all@meta.data$genotype.collapsed[MPN.integration.all@meta.data$genotype.classification=="preleukemic" & MPN.integration.all@meta.data$donor.type=="AP"]<-"preLSC"

MPN.integration.all@meta.data$genotype.collapsed[MPN.integration.all@meta.data$genotype.classification=="preleukemic" & MPN.integration.all@meta.data$donor.type=="CP"]<-"TP53_WT"
MPN.integration.all@meta.data$genotype.collapsed[MPN.integration.all@meta.data$genotype.classification=="preleukemic" & MPN.integration.all@meta.data$donor.type=="MAJIC"]<-"TP53_WT"

MPN.integration.all@meta.data$genotype.collapsed[MPN.integration.all@meta.data$donor.type=="MF"]<-"MF"
MPN.integration.all@meta.data$genotype.collapsed[MPN.integration.all@meta.data$donor.type=="normal"]<-"WT"

#Include louvain clusters of integrated dataset (only for APMPNAML patients)
merged_louvain<-read.table(file=paste(path.dropbox,"../classification/louvain/louvain.Harmony.APMPNAML.revised.txt",sep=""),sep="\t",header = T)

merged_louvain_cl1<-subset(merged_louvain,merged_louvain=="cl1")
merged_louvain_cl2<-subset(merged_louvain,merged_louvain=="cl2")
merged_louvain_cl3<-subset(merged_louvain,merged_louvain=="cl3")

MPN.integration.all@meta.data$merged.louvain[MPN.integration.all@meta.data$donor.type=="MF"]<-"MF"
MPN.integration.all@meta.data$merged.louvain[MPN.integration.all@meta.data$donor.type=="normal"]<-"normal"

MPN.integration.all@meta.data$merged.louvain[MPN.integration.all@meta.data$donor.type=="AP" & MPN.integration.all@meta.data$Cell %in% merged_louvain_cl1$Cell]<-"cl1_LSC"
MPN.integration.all@meta.data$merged.louvain[MPN.integration.all@meta.data$donor.type=="AP" & MPN.integration.all@meta.data$Cell %in% merged_louvain_cl2$Cell]<-"cl2_erythroid"
MPN.integration.all@meta.data$merged.louvain[MPN.integration.all@meta.data$donor.type=="AP" & MPN.integration.all@meta.data$Cell %in% merged_louvain_cl3$Cell]<-"cl3_preleukemic"

#####################################
  
get_variable_genes_by_fitting_GLM_model(MPN.integration.all,mean_expr_cutoff = 1,disp_zscore_cutoff = 0.1,
                                        quantile_genes_expr_for_fitting = 0.6,quantile_genes_cv2_for_fitting = 0.2)
#"Identified :4564 variable genes"
######################################################
runPCA(MPN.integration.all,use.components=50,use.regressout.data = F)
#################

#################
save(MPN.integration.all,file="Robjects/MFp53MPNAML_integration.revised.rdata")
#################
