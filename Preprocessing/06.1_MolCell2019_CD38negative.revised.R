###############################################################################
library(Matrix)
library(matrixStats)
##################################
colData.qc<-read.table(file="HT_TARGET_MF_colData.qc.csv",sep=",",header = T)

counts.MolCell2019<-read.table(file="HT_TARGET_MF_counts.qc.tsv",sep="\t",header = T)
rownames(counts.MolCell2019) <- counts.MolCell2019$Gene
###################################
cell.anno<-colData.qc[colData.qc$cell_id %in% colnames(counts.MolCell2019),]
cell.anno_CD38neg_TARGET<-subset(cell.anno,CD38<4000) #no selection for CD90 expression

cell.anno_CD38neg_TARGET$genotype.classification<-"MF"
normal<-c("Aph1","HD85")
cell.anno_CD38neg_TARGET$genotype.classification[cell.anno_CD38neg_TARGET$donor %in% normal]<-"WT"

counts.qc_CD38neg_TARGET<-counts.MolCell2019[,as.character(cell.anno_CD38neg_TARGET$cell_id)] #1431 cells

####################################

save(counts.qc_CD38neg_TARGET,file="../classification/preleukemic_integration/counts.qc_CD38neg_TARGET.revised.rdata")
save(cell.anno_CD38neg_TARGET,file="../classification/preleukemic_integration/cell.anno_CD38neg_TARGET.revised.rdata")
###################################

