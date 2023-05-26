## Setup notes
# Seurat 2.3.0
# Check if java command line installed (javac -version) on terminal, if not, download jdk-13.0.2_osx-x64_bin.dmg for Mac from https://www.oracle.com/java/technologies/javase/jdk13-archive-downloads.html

# Load additional packages
library(SingCellaR)

#Projection of preleukemic MPNAML cells into MF dataset
#LSI integration MF hierarchy
#Author:Alba Rodriguez-Meira and Sean Wen
#Date: 19th November 2019
#Modified last: 13th October 2021
####################################################

library(Matrix)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(tidyverse)
library(uwot)
library(edgeR)
library(FNN)
library(matrixStats)
library(Rcpp)
library(bigmemory)
####################################################
set.seed(1)
source("../../../Sources/LSI/LSI_functions.R")
####################################################
#Input Data
####################################################
#Read in Summarized Experiment
#Please Note Code here has been modified to work with finalized summarized experiment

#SE Healthy Cells
seReference <- readRDS("../../../Sources/LSI/MF_reference_RS.rds")
class(seReference) # SingCellaR_Int
#Reference Summarized Experiment
seDisease <- readRDS("../../../Sources/LSI/HTMPNAML_LSI.revised.rds")
class(seDisease) # SummarizedExperiment

#Subset preleukemic cells from AP donors only
preleukemic.cells<-read.table("../../../Sources/preleukemic_integration_metadata.qc.revised.csv",sep=",",header = T)
AP<-c("GR003_AP","GR001","GR002","GR005_AP","GR006_AP","GR007_AP","IF0131","IF0318","IF0391","SB5702")
preleukemic.cells.AP<-subset(preleukemic.cells,stage %in% AP)
id <- preleukemic.cells.AP$cell_id #880 cells
seDisease <- seDisease[,colData(seDisease)$Cell %in% id]

#Identify Gene Universe
gU <- intersect(row.names(seReference@genes.info), row.names(seDisease))
length(row.names(seReference@genes.info)); length(row.names(seDisease)); length(gU) #19385 common genes
gU <- gU[!grepl("^MT", gU)]
length(gU) #19312

#Set Clustering Parameters
resolution <- c(0.2,0.8,0.8) #clustering resolution
varGenesToUse <- c(1000,1000,1000) #number of variable genes

#Optimize LSI Features
seReference_mat <- assays(seReference)$counts
row.names(seReference_mat) <- row.names(seReference@metadata$optimizeLSI$iter3$matNorm)
names(seReference_mat) <- names(seReference@metadata$optimizeLSI$iter3$matNorm)
matAll <- cbind(seReference_mat[gU,], assay(seDisease[gU,]))
lsiObj <- optimizeLSI(matAll, resolution = resolution, varFeatures = varGenesToUse)

#UMAP
set.seed(1)
umap <- uwot::umap(
  lsiObj[[length(lsiObj)]]$lsiObj$matSVD[,1:25], 
  n_neighbors = 30, 
  min_dist = 0.5, 
  metric = "euclidean", 
  n_threads = 5, 
  verbose = TRUE, 
  ret_model = FALSE
)

#Plot Info
cells <- c(rep("reference", ncol(seReference)),rep("disease",ncol(seDisease)))
splitCells <- split(cells,lsiObj[[length(lsiObj)]]$clusters)
df <- data.frame(
  clusters = names(splitCells),
  proportion = unlist(lapply(seq_along(splitCells), function(x) sum(splitCells[[x]]=="disease") / length(splitCells[[x]])))
)

#Plot UMAP Data Frame
plotDF <- data.frame(umap)
rownames(plotDF) <- c(colnames(seReference), colnames(seDisease))
plotDF$type <- cells
plotDF$clusters <- lsiObj[[length(lsiObj)]]$clusters
plotDF$classification <- 0
plotDF$classification[plotDF$type == "disease" & plotDF$clusters %in% paste0(df$clusters[df[,2] > 0.8])] <- 1
plotDF$classification[plotDF$type == "disease"] <- plotDF$classification[plotDF$type == "disease"] + 1
plotDF <- plotDF[order(plotDF$classification), ]

#Formal Classification
plotDF$classificationSTR <- "reference"
plotDF$classificationSTR[plotDF$classification==1] <- "disease" #previously "healthy-like"
plotDF$classificationSTR[plotDF$classification==2] <- "disease" #previously "disease-like"

####################################################
#Project Into LSI UMAP
####################################################

#Previous Reference Summarized Experiment
se <- seReference

#Load Saved UMAP Manifold
umapManifold <- uwot::load_uwot("../../../Sources/LSI/MF_reference-UMAP-model_RS.uwot")

#LSI Projection Matrix
lsiGenes <- metadata(se)$variableGenes
matProjectLSI <- assay(seDisease[lsiGenes,])

#LSI Project
lsiReference <- metadata(se)$optimizeLSI[[length(metadata(se)$optimizeLSI)]]$lsiObj
lsiProjection <- projectLSI(matProjectLSI, lsiReference)

#UMAP Projection
#Set Seed Prior to umap_transform (see uwot github)
set.seed(1)
umapProjection <- uwot::umap_transform(as.matrix(lsiProjection)[,1:25], umapManifold, verbose = TRUE)

#Plot Projection
refDF <- data.frame(row.names = colnames(se), X1 = umapManifold$embedding[,1], X2 = umapManifold$embedding[,2], Type = "reference")
proDF <- data.frame(row.names = colnames(seDisease), X1 = umapProjection[,1], X2 = umapProjection[,2], Type = seDisease$stage)

projectionDF <- rbind(refDF, proDF)

## Separate reference and p53 samples

projectionDF_reference<-subset(projectionDF,Type=="reference")
projectionDF_NOTreference<-subset(projectionDF,Type!="reference")
projectionDF_NOTreference$Type<-"preLSCs"

p<-ggplot(projectionDF_reference, aes(X1,X2,color=Type)) + 
  geom_point() +
  theme_classic()+
  theme(legend.position = "none", axis.title =element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank())+
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") 

p+stat_density_2d(data=projectionDF_NOTreference,aes(X1,X2,color=Type))+
  scale_color_manual(values=c("blue","lightgrey"))

ggsave(paste(path.local,'Figure3b_right.png',sep=""),height=5,width=5)
