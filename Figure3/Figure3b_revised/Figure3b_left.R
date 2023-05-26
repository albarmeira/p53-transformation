## Setup notes
#Seurat 2.3.0 is required to run LSI projection scripts. To install this version of Seurat, run the following commands:
#install.packages('remotes')
#remotes::install_version(package = 'Seurat', version = package_version('2.3.0'))
#Java is also required. Check if java is installed by typing (javac -version) on the command line terminal. If not, download jdk-13.0.2_osx-x64_bin.dmg for Mac from https://www.oracle.com/java/technologies/javase/jdk13-archive-downloads.html

####################################################
#Mapping preLSCs from AP-HTMPNAML dataset into Greenleaf's LSI human hematopoetic map
#LSI integration
#Author:Alba Rodriguez-Meira and Sean Wen
#Date: 19th November 2019
#Modified last: 12th March 2023
####################################################

library(Matrix)
library(SummarizedExperiment)
library(tidyverse)
library(uwot)
library(edgeR)
library(FNN)
library(matrixStats)
library(Rcpp)
set.seed(1)

source("../../../Sources/LSI/LSI_functions.R")

####################################################
#Input Data
####################################################

#Import Healthy Cells reference (generated in 01_Greenleaf_Healthy_Reference.R)
seReference <- readRDS("../../../Sources/LSI/scRNA-Healthy-Hematopoiesis.revised.rds")

#Read in Summarized Experiment from HTMPNAML dataset (generated in 02_prepare_SummExp.R)
seDisease <- readRDS("../../../Sources/LSI/HTMPNAML_LSI.revised.rds")

#Subset preleukemic cells from AP donors only
preleukemic.cells<-read.table("../../../Sources/preleukemic_integration_metadata.qc.revised.csv",sep=",",header = T)
AP<-c("GR003_AP","GR001","GR002","GR005_AP","GR006_AP","GR007_AP","IF0131","IF0318","IF0391","SB5702")
preleukemic.cells.AP<-subset(preleukemic.cells,stage %in% AP)
id <- preleukemic.cells.AP$cell_id
seDisease <- seDisease[,colData(seDisease)$Cell %in% id]

#Identify Gene Universe
gU <- intersect(rownames(seReference), rownames(seDisease))
gU <- gU[!grepl("^MT", gU)]

#Set Clustering Parameters
resolution <- c(0.2,0.8,0.8) #clustering resolution
varGenesToUse <- c(1000,1000,1000) #number of variable genes

#Optimize LSI Features
matAll <- cbind(assay(seReference[gU,]), assay(seDisease[gU,]))
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

#Plot UMAP Data Frame
plotDF <- data.frame(umap)
rownames(plotDF) <- c(colnames(seReference), colnames(seDisease))
plotDF$type <- cells
plotDF$clusters <- lsiObj[[length(lsiObj)]]$clusters
plotDF$classification <- 0
plotDF$classification[plotDF$type == "disease"]<- 1
plotDF <- plotDF[order(plotDF$classification), ]

#Formal Classification
plotDF$classificationSTR <- "reference"
plotDF$classificationSTR[plotDF$classification==1] <- "disease"
table(plotDF$classificationSTR) #880 preleukemic cells used for plotting; 35582 cells used as reference

####################################################
#Project Into LSI UMAP
####################################################

#Previous Reference Summarized Experiment
se <- seReference

#Load Saved UMAP Manifold
umapManifold <- uwot::load_uwot("../../../Sources/LSI/scRNA-Healthy-Hematopoiesis.revised-UMAP-model.uwot")

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
proDF <- data.frame(row.names = colnames(seDisease), X1 = umapProjection[,1], X2 = umapProjection[,2], Type = seDisease$sampleID)
projectionDF <- rbind(refDF, proDF)

projectionDF$color[projectionDF$Type=="reference"]<-"reference"
projectionDF$color[projectionDF$Type!="reference"]<-"preLSCs"

projectionDF_reference<-subset(projectionDF,Type=="reference")
projectionDF_NOTreference<-subset(projectionDF,Type!="reference")

p<-ggplot(projectionDF_reference, aes(X1,-X2,color=color)) + 
  geom_point(size=0.5,alpha=0.5) +
  theme_classic()+
  theme(legend.position = "none", axis.title =element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank())+
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") 

p+stat_density_2d(data=projectionDF_NOTreference,aes(X1,-X2,color=color),size=0.5)+
  scale_color_manual(values=c("blue","lightgrey"))

ggsave('Figure3b_left.png',height=6,width=6)
