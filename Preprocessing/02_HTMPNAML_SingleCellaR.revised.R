#High throughput scRNA-seq of MPNAML samples - SingleCellaR processing
#Data preprocessing and QC metrics
#Author:Alba Rodriguez-Meira and sean Wen
#Date: 6th March 2021
#Modified last: 26th February 2023
####################################################

###################################
library(SingCellaR) #1.2.1
library(SingleCellExperiment) # 1.12.0
library(Rcpp) # 1.0.7
library(Matrix) # 1.5-1
library(matrixStats) # 0.58.0
library(bigmemory) # 4.5.36
library(LinkedMatrix) # 1.4.0
library(irlba) # 2.3.3
library(Rtsne) # 0.15
library(ggplot2) # 3.3.6
library(gridExtra) # 2.3
library(cccd) # 1.5
library(ggpubr) # 0.4.0
library(statmod) # 1.4.36
library(umap) # 0.2.7.0
library(reticulate) # 1.22
library(data.table)
library(pbapply) # 1.4-3
library(fgsea) # 1.16.0
library(limma) # 3.46.0
library(RANN) # 2.6.1
library(ComplexHeatmap) # 2.11.1
library(circlize) # 0.4.12
library(threejs) # 0.3.3
library(RColorBrewer) # 1.1-2
library(igraph) # 1.2.6
library(pheatmap) # 1.0.12
library(AUCell) # 1.12.0
library(diffusionMap) # 1.2.0
library(destiny) # 3.9.0
##################################

source("../../../Sources/SingleCellRNASeq_Edit.R")
setClass("TargetSeq",representation(GenesExpressionMatrixFile="character",CellsMetaDataFile="character"))

HTMPNAML<-new("TargetSeq",GenesExpressionMatrixFile="p53_HTMPNAML_counts.qc.all.revised.new.id.tsv",CellsMetaDataFile="p53_HTMPNAML_colData.qc.final.revised.new.id.csv")

######Load gene expression#######
load_gene_expression_from_a_file(HTMPNAML,isTargetSeq = T,sep="\t")

######Automatically process cell annotation###
TargetSeq_process_cells_annotation(HTMPNAML,mito_genes_start_with = "MT-",ERCC_genes_start_with = "ERCC-")
######Identify good cells########
TargetSeq_plot_cells_annotation(HTMPNAML,type = "boxplot")

# Filter low quality cells (qc performed previously; not needed explicitly)
TargetSeq_filter_cells_and_genes(HTMPNAML,min_Reads = 2000,min_detected_genes = 200,
                                 max_percent_mito = 80,max_percent_ERCC = 50) #initially MPNAML13PL72_17B was excluded for the number of genes; this cell should be included

######Load external cell metadata (e.g., cell's genotype)
md.1 <- HTMPNAML@meta.data
md.2 <- as.data.frame(fread("p53_HTMPNAML_colData.qc.final.revised.new.id.csv", sep=",", header=TRUE, stringsAsFactors=FALSE))
names(md.2)[which(names(md.2)=="cell_id")] <- "Cell"
md.1 <- plyr::join(md.1, md.2, by="Cell", type="left")
HTMPNAML@meta.data <- md.1

#####Data normalization##############
normalize_UMIs(HTMPNAML,use.scaled.factor = F)
# Add cell cycle score
add_cell_cycle_genes_score(HTMPNAML,gmt.file = paste(path.dropbox, "MC.human.signature.genes.gmt", sep=""))
#############################################################################
remove_unwanted_confounders(HTMPNAML,residualModelFormulaStr="~detectedGenesPerCell+UMI_count")
#############################################################################
get_variable_genes_by_fitting_GLM_model(HTMPNAML,mean_expr_cutoff = 1,disp_zscore_cutoff = 0.1,
                                        quantile_genes_expr_for_fitting = 0.6,quantile_genes_cv2_for_fitting = 0.2)

#"Identified :4578 variable genes"
remove_unwanted_genes_from_variable_gene_set(HTMPNAML,gmt.file = paste(path.dropbox, "MC.human.signature.genes.gmt",sep=""),
                                             removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene"))
#"25 genes are removed from the variable gene set."

plot_variable_genes(HTMPNAML,quantile_genes_expr_for_fitting = 0.6,quantile_genes_cv2_for_fitting = 0.2)
#####################################
runPCA(HTMPNAML,use.components=20,use.regressout.data = F)
######################################
#plot_PCA_Elbowplot(HTMPNAML)
######################################
runUMAP(HTMPNAML,n.dims.use=15,n.neighbors=5,uwot.min.dist=0.30,uwot.metric = "correlation",n.seed = 1,
        useIntegrativeEmbeddings=FALSE,
        dim_reduction_method="pca")
plot_umap_label_by_a_feature_of_interest(HTMPNAML,feature = "cell_type")
plot_umap_label_by_a_feature_of_interest(HTMPNAML,feature = "donor")
######################################
plot_umap_label_by_multiple_gene_sets(HTMPNAML,gmt.file=paste(path.dropbox,"MC.human.signature.genes.gmt", sep=""),show_gene_sets=c("LSC","Ery_AML"),
                                      custom_color=c("orange","red"))

plot_umap_label_by_multiple_gene_sets(HTMPNAML,gmt.file = paste(path.dropbox, "SingleCellaR.0.1.4/Data/genesets/MC.human.signature.genes.gmt", sep=""),
                                      show_gene_sets = c("Erythroid","Myeloid","Lymphoid","Megakaryocyte","HSPC"),
                                      custom_color = c("red","green","yellow","blue","orange"),
                                      isNormalizedByHouseKeeping = T)

######################################
save(HTMPNAML,file="Robjects/HTMPNAML_final.revised.rdata")
######################################
