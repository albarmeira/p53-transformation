## Setup notes
#Seurat 2.3.0 is required to run LSI projection scripts. To install this version of Seurat, run the following commands:
#install.packages('remotes')
#remotes::install_version(package = 'Seurat', version = package_version('2.3.0'))
#Java is also required. Check if java is installed by typing (javac -version) on the command line terminal. If not, download jdk-13.0.2_osx-x64_bin.dmg for Mac from https://www.oracle.com/java/technologies/javase/jdk13-archive-downloads.html

####################################################

#Mapping p53-MPNAML dataset into Greenleaf's LSI projection
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
seReference <- readRDS("../../../Sources/LSI/scRNA-Healthy-Hematopoiesis.rds")

#Read in Summarized Experiment from HTMPNAML dataset (generated in 02_prepare_SummExp.R)
seDisease <- readRDS("../../../Sources/LSI/HTMPNAML_LSI.revised.rds")

#Subset p53 mutant cells from AP donors only
TP53.multihit<-subset(colData(seDisease),genotype.classification=="TP53_multihit_HOM"|genotype.classification=="TP53_multihit_M2")
AP<-c("GH001_003","GR003_AP","IF0392","GST010","GR004_AP","IF0131","JB4211",
      "IF0393","GR001","GR005_AP","GR007_AP","GR006_AP","IF0391","GR002","SB5702","IF0318")
TP53.multihit<-subset(TP53.multihit,stage %in% AP)
id <- TP53.multihit$Cell
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
plotDF$classificationSTR[plotDF$classification==1] <- "disease" #previously healthy-like
plotDF$classificationSTR[plotDF$classification==2] <- "disease" #previously disease-like
table(plotDF$classificationSTR)  #7881 disease and 35582 reference

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
projectionDF$color[projectionDF$Type!="reference"]<-"p53MPNAML"

projectionDF_reference<-subset(projectionDF,Type=="reference")
projectionDF_NOTreference<-subset(projectionDF,Type!="reference")

p<-ggplot(projectionDF_reference, aes(X1,-X2,color=color)) + 
  geom_point(size=0.5,alpha=0.5) +
  theme_classic()+
  theme(legend.position = "none", axis.title =element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank())+
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2") 

p+stat_density_2d(data=projectionDF_NOTreference,aes(X1,-X2,color=color),size=0.5, bins=15)+
  scale_color_manual(values=c("red","lightgrey"))

#ggsave('Figure2d_Panel1.png',height=6,width=6)

##############################################

# utils::sessionInfo()
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Rcpp_1.0.4                  FNN_1.1.3                   edgeR_3.28.1               
# [4] limma_3.42.2                uwot_0.1.8                  forcats_0.5.0              
# [7] stringr_1.4.0               dplyr_1.0.7                 purrr_0.3.4                
# [10] readr_1.3.1                 tidyr_1.0.2                 tibble_3.0.1               
# [13] ggplot2_3.3.5               tidyverse_1.3.0             SummarizedExperiment_1.16.1
# [16] DelayedArray_0.12.3         BiocParallel_1.20.1         matrixStats_0.56.0         
# [19] Biobase_2.46.0              GenomicRanges_1.38.0        GenomeInfoDb_1.22.1        
# [22] IRanges_2.20.2              S4Vectors_0.24.4            BiocGenerics_0.32.0        
# [25] Matrix_1.2-17              
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.1.4             reticulate_1.15        R.utils_2.9.2          tidyselect_1.1.1      
# [5] htmlwidgets_1.5.1      grid_3.6.1             Rtsne_0.15             munsell_0.5.0         
# [9] codetools_0.2-16       mutoss_0.1-12          ica_1.0-2              withr_2.2.0           
# [13] colorspace_1.4-1       knitr_1.28             rstudioapi_0.11        Seurat_2.3.4          
# [17] ROCR_1.0-7             robustbase_0.93-6      dtw_1.21-3             gbRd_0.4-11           
# [21] labeling_0.3           Rdpack_0.11-1          lars_1.2               GenomeInfoDbData_1.2.2
# [25] mnormt_1.5-7           farver_2.0.3           bit64_0.9-7            vctrs_0.3.8           
# [29] generics_0.0.2         TH.data_1.0-10         xfun_0.13              diptest_0.75-7        
# [33] R6_2.4.1               isoband_0.2.1          locfit_1.5-9.4         hdf5r_1.3.2           
# [37] flexmix_2.3-15         bitops_1.0-6           assertthat_0.2.1       SDMTools_1.1-221.1    
# [41] scales_1.1.0           multcomp_1.4-13        nnet_7.3-12            gtable_0.3.0          
# [45] npsurv_0.4-0           sandwich_2.5-1         rlang_0.4.11           splines_3.6.1         
# [49] acepack_1.4.1          broom_0.5.6            checkmate_2.0.0        yaml_2.2.1            
# [53] reshape2_1.4.4         modelr_0.1.7           backports_1.1.6        Hmisc_4.4-0           
# [57] tools_3.6.1            ellipsis_0.3.2         gplots_3.0.3           RColorBrewer_1.1-2    
# [61] proxy_0.4-24           ggridges_0.5.2         TFisher_0.2.0          plyr_1.8.6            
# [65] base64enc_0.1-3        zlibbioc_1.32.0        RCurl_1.98-1.2         rpart_4.1-15          
# [69] pbapply_1.5-0          cowplot_1.0.0          zoo_1.8-8              haven_2.2.0           
# [73] cluster_2.1.0          fs_1.4.1               magrittr_1.5           RSpectra_0.16-0       
# [77] data.table_1.12.8      lmtest_0.9-37          reprex_0.3.0           RANN_2.6.1            
# [81] mvtnorm_1.1-0          fitdistrplus_1.0-14    hms_0.5.3              lsei_1.2-0            
# [85] jpeg_0.1-8.1           mclust_5.4.6           readxl_1.3.1           gridExtra_2.3         
# [89] compiler_3.6.1         KernSmooth_2.23-15     crayon_1.3.4           R.oo_1.23.0           
# [93] htmltools_0.5.2        segmented_1.1-0        Formula_1.2-3          snow_0.4-3            
# [97] lubridate_1.7.8        DBI_1.1.0              dbplyr_1.4.3           MASS_7.3-51.4         
# [101] fpc_2.2-7              cli_2.0.2              R.methodsS3_1.8.0      gdata_2.18.0          
# [105] metap_1.3              igraph_1.2.5           pkgconfig_2.0.3        sn_1.6-1              
# [109] numDeriv_2016.8-1.1    foreign_0.8-71         xml2_1.3.2             foreach_1.5.1         
# [113] multtest_2.42.0        XVector_0.26.0         bibtex_0.4.2.2         rvest_0.3.5           
# [117] digest_0.6.25          RcppAnnoy_0.0.16       tsne_0.1-3             cellranger_1.1.0      
# [121] htmlTable_1.13.3       kernlab_0.9-29         gtools_3.8.2           modeltools_0.2-23     
# [125] lifecycle_1.0.1        nlme_3.1-140           jsonlite_1.6.1         fansi_0.4.1           
# [129] pillar_1.6.3           lattice_0.20-38        fastmap_1.1.0          httr_1.4.1            
# [133] plotrix_3.7-8          DEoptimR_1.0-8         survival_3.1-12        glue_1.4.0            
# [137] png_0.1-7              prabclus_2.3-2         iterators_1.0.13       bit_4.0.4             
# [141] class_7.3-15           stringi_1.4.6          mixtools_1.2.0         doSNOW_1.0.18         
# [145] latticeExtra_0.6-29    caTools_1.18.0         irlba_2.3.3            ape_5.3 
