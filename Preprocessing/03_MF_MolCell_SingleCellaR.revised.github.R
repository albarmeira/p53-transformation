####################################################
#High throughput scRNA-seq of MF MolCell 2019 samples - SingleCellaR processing
#Data integration
#Author:Alba Rodriguez-Meira and Sean Wen
#Date: 24th March 2021
#Modified last: 26th February 2023
####################################################

###################################
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
library(threejs)
library(RColorBrewer)
library(igraph)
library(pheatmap)
library(AUCell)
library(diffusionMap)
library(destiny)
##################################
source("../../../Sources/SingleCellRNASeq_Edit.R")
setClass("TargetSeq",representation(GenesExpressionMatrixFile="character",CellsMetaDataFile="character"))

MF_MolCell<-new("TargetSeq",GenesExpressionMatrixFile=paste(path.dropbox, "HT_TARGET_MF_counts.qc.tsv", sep=""),
               CellsMetaDataFile=paste(path.dropbox,"HT_TARGET_MF_colData.qc.csv",sep=""))

######Load gene expression#######
load_gene_expression_from_a_file(MF_MolCell,isTargetSeq = T,sep="\t")
######Automatically process cell annotation###
TargetSeq_process_cells_annotation(MF_MolCell,mito_genes_start_with = "MT-",ERCC_genes_start_with = "ERCC-")
######Identify good cells########
TargetSeq_plot_cells_annotation(MF_MolCell,type = "boxplot")

# Filter low quality cells (qc performed previously; not needed explicitly)
TargetSeq_filter_cells_and_genes(MF_MolCell,min_Reads = 2000,min_detected_genes = 200,
                                 max_percent_mito = 80,max_percent_ERCC = 50) 
#"0/2734 cells will be filtered out from the downstream analyses!."

######Load external cell metadata (e.g., cell's genotype)
md.1 <- MF_MolCell@meta.data
md.2 <- as.data.frame(fread(paste(path.dropbox,"HT_TARGET_MF_colData.qc.csv",sep=""), sep=",", header=TRUE, stringsAsFactors=FALSE))
names(md.2)[which(names(md.2)=="cell_id")] <- "Cell"
md.1 <- plyr::join(md.1, md.2, by="Cell", type="left")
MF_MolCell@meta.data <- md.1

#####Data normalization##############
normalize_UMIs(MF_MolCell,use.scaled.factor = F)

# Add cell cycle score
add_cell_cycle_genes_score(MF_MolCell,gmt.file = paste(path.dropbox,"MC.human.signature.genes.gmt",sep=""))

#############################################################################
remove_unwanted_confounders(MF_MolCell,residualModelFormulaStr="~detectedGenesPerCell+UMI_count") 
#############################################################################
get_variable_genes_by_fitting_GLM_model(MF_MolCell,mean_expr_cutoff = 1,disp_zscore_cutoff = 0.1,
                                        quantile_genes_expr_for_fitting = 0.6,quantile_genes_cv2_for_fitting = 0.2)
#"Identified :4002 variable genes"

remove_unwanted_genes_from_variable_gene_set(MF_MolCell,gmt.file = paste(path.dropbox,"MC.human.signature.genes.gmt",sep=""),
                                             removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene"))
#"6 genes are removed from the variable gene set."

plot_variable_genes(MF_MolCell,quantile_genes_expr_for_fitting = 0.6,quantile_genes_cv2_for_fitting = 0.2)
#####################################
runPCA(MF_MolCell,use.components=20,use.regressout.data = F)
######################################
plot_PCA_Elbowplot(MF_MolCell)
######################################
runUMAP(MF_MolCell,n.dims.use=15,n.neighbors=5,uwot.min.dist=0.30,uwot.metric = "correlation",n.seed = 1,
        useIntegrativeEmbeddings=FALSE,
        dim_reduction_method="pca")

plot_umap_label_by_a_feature_of_interest(MF_MolCell,feature = "donor")
plot_umap_label_by_a_feature_of_interest(MF_MolCell,feature = "stage")

plot_umap_label_by_multiple_gene_sets(MF_MolCell,gmt.file = paste(path.dropbox,"MC.human.signature.genes.gmt",sep=""),
                                      show_gene_sets = c("Erythroid","Myeloid","Lymphoid","Megakaryocyte","HSPC"),
                                      custom_color = c("red","green","yellow","blue","orange"),
                                      isNormalizedByHouseKeeping = T)

save(MF_MolCell,file="Robjects/MF_MolCell2019.revised.rdata")
######################################
