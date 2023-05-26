#Projection of p53-mutant MPNAML cells into MF dataset
#LSI integration MF hierarchy
#Author:Alba Rodriguez-Meira and Sean Wen
#Date: 19th November 2019
#Modified last: 16th March 2023
####################################################
installation of package ‘qqconf’ had non-zero exit status
installation of package ‘Hmisc’ had non-zero exit status
installation of package ‘metap’
## Setup notes
#Seurat 2.3.0 is required to run LSI projection scripts. To install this version of Seurat, run the following commands:
#install.packages('remotes')
#remotes::install_version(package = 'Seurat', version = package_version('2.3.0'))
#Java is also required. Check if java is installed by typing (javac -version) on the command line terminal. If not, download jdk-13.0.2_osx-x64_bin.dmg for Mac from https://www.oracle.com/java/technologies/javase/jdk13-archive-downloads.html

#Set working directory
path.dropbox <- "/Users/seanwen/Dropbox/p53_paper_IAD/Figures_RAW/FigureS5/FigureS5e_revised/"
path.local <- "/Users/seanwen/Documents/Alba/data_check/Figures_RAW/FigureS5/FigureS5e_revised/"
setwd(path.dropbox)

####################################################
library(SingCellaR)
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
#Read in Summarized Experiment
#Please Note Code here has been modified to work with finalized summarized experiment

#SE Healthy Cells
seReference <- readRDS("../../../Sources/LSI/MF_reference_RS.rds")

#Reference Summarized Experiment
seDisease <- readRDS("../../../Sources/LSI/HTMPNAML_LSI.revised.rds")

#Subset p53 mutant cells from AP donors only
TP53.multihit<-subset(colData(seDisease),genotype.classification=="TP53_multihit_HOM"|genotype.classification=="TP53_multihit_M2")
AP<-c("GH001_003","GR003_AP","IF0392","GST010","GR004_AP","IF0131","JB4211",
      "IF0393","GR001","GR005_AP","GR007_AP","GR006_AP","IF0391","GR002","SB5702","IF0318")
TP53.multihit<-subset(TP53.multihit,stage %in% AP)
id <- TP53.multihit$Cell
seDisease <- seDisease[,colData(seDisease)$Cell %in% id]

#Identify Gene Universe
gU <- intersect(row.names(seReference@genes.info), row.names(seDisease))
length(row.names(seReference@genes.info)); length(row.names(seDisease)); length(gU)
gU <- gU[!grepl("^MT", gU)]

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
#If disease cells are clustered with healthy cluster (proportion > 0.8) we will classify these as healthy-like
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
projectionDF_NOTreference$Type<-"p53MPNAML"

p<-ggplot(projectionDF_reference, aes(X1,X2,color=Type)) + 
  geom_point() +
  theme_classic()+
  theme(legend.position = "none", axis.title =element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank())+
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") 

p+stat_density_2d(data=projectionDF_NOTreference,aes(X1,X2,color=Type))+
  scale_color_manual(values=c("red","lightgrey"))

ggsave(paste(path.local,'FigureS5e.png',sep=""),height=5,width=5)

