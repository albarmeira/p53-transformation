# Set working directories
path.dropbox<-"/Users/miguelcachosoblechero/Dropbox (Partners HealthCare)/p53MPNAML/p53_paper_RawData/01.Preprocessing/HT_transcriptome_revised/06-SinglCellaR/"
#path.dropbox <-"/Users/alba/Dropbox (Partners HealthCare)/p53MPNAML/p53_paper_RawData/01.Preprocessing/HT_transcriptome_revised/06-SinglCellaR/"
#path.dropbox <- "/Users/seanwen/Dropbox/p53_paper_IAD/01.Preprocessing/HT_transcriptome_revised/06-SinglCellaR/"
#path.local <- "/Users/seanwen/Documents/Alba/data_check/01.Preprocessing/HT_transcriptome_revised/06-SinglCellaR/"
#setwd(path.local)

# Load additional packages
library(data.table)

#SinglCellaR processing of preleukemic integration (ND,MF,AML-AP)
#Data integration
#Author:Alba Rodriguez-Meira and Sean Wen
#Date: 16th October 2021 
#Modified last: 4th March 2023
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
library(SingCellaR) #1.2.1
##################################
#source("../../../Sources/SingleCellRNASeq_Edit.R")
setClass("TargetSeq",representation(GenesExpressionMatrixFile="character",CellsMetaDataFile="character"))

preleukemic.integration<-new("TargetSeq",GenesExpressionMatrixFile="../classification/preleukemic_integration/counts_preleukemic.integration.revised.tsv",
               CellsMetaDataFile="../classification/preleukemic_integration/preleukemic_integration_metadata.qc.revised.csv")

######Load gene expression#######
load_gene_expression_from_a_file(preleukemic.integration,isTargetSeq = T,sep="\t")
######Automatically process cell annotation###
#TargetSeq_process_cells_annotation(preleukemic.integration,mitochondiral_genes_start_with = "MT-",ERCC_genes_start_with = "ERCC-")
TargetSeq_process_cells_annotation(preleukemic.integration,mito_genes_start_with = "MT-",ERCC_genes_start_with = "ERCC-")
######Identify good cells########
#TargetSeq_plot_cells_annotation(preleukemic.integration,type = "boxplot")

# Filter low quality cells (qc performed previously; not needed explicitly)
TargetSeq_filter_cells_and_genes(preleukemic.integration,min_Reads = 2000,min_detected_genes = 200,
                                 max_percent_mito = 80,max_percent_ERCC = 50) 

#"0/2716 cells will be filtered out from the downstream analyses!."

######Load external cell metadata (e.g., cell's genotype)
md.1 <- preleukemic.integration@meta.data
md.2 <- read.table("../classification/preleukemic_integration/preleukemic_integration_metadata.qc.revised.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
names(md.2)[which(names(md.2)=="cell_id")] <- "Cell"
md.1 <- plyr::join(md.1, md.2, by="Cell", type="left")
preleukemic.integration@meta.data <- md.1

#####Data normalization##############
normalize_UMIs(preleukemic.integration,use.scaled.factor = F)

#Error in `assays<-`(`*tmp*`, withDimnames = withDimnames, ..., value = `*vtmp*`) : 
  #nb of rows in 'assay' (32830) must equal nb of rows in 'rowData' (1)

#Add cell cycle score
add_cell_cycle_genes_score(preleukemic.integration,gmt.file =paste(path.dropbox, "SingleCellaR.0.1.4/Data/genesets/MC.human.signature.genes.gmt", sep=""))

#############################################################################
remove_unwanted_confounders(preleukemic.integration,residualModelFormulaStr="~detectedGenesPerCell+UMI_count") 
#############################################################################
get_variable_genes_by_fitting_GLM_model(preleukemic.integration,mean_expr_cutoff = 1,disp_zscore_cutoff = 0.1,
                                        quantile_genes_expr_for_fitting = 0.6,quantile_genes_cv2_for_fitting = 0.2)

remove_unwanted_genes_from_variable_gene_set(preleukemic.integration,gmt.file = paste(path.dropbox,"SingleCellaR.0.1.4/Data/genesets/human.signature.genes.gmt",sep=""),
                                             removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene"))

plot_variable_genes(preleukemic.integration,quantile_genes_expr_for_fitting = 0.6,quantile_genes_cv2_for_fitting = 0.2)
#####################################
runPCA(preleukemic.integration,use.components=30,use.regressout.data = F)
######################################
plot_PCA_Elbowplot(preleukemic.integration)
######################################
runUMAP(preleukemic.integration,n.dims.use=10,n.neighbors=5,uwot.min.dist =0.30,uwot.metric  = "correlation",n.seed = 1,
        useIntegrativeEmbeddings=FALSE,
        dim_reduction_method="pca")

plot_umap_label_by_a_feature_of_interest(preleukemic.integration,feature="genotype.classification",point.size=2)
plot_umap_label_by_a_feature_of_interest(preleukemic.integration,feature="donor",point.size=2)
######################################

runFA2_ForceDirectedGraph(preleukemic.integration,n.dims.use = 10)
p<-plot_forceDirectedGraph_label_by_a_feature_of_interest(preleukemic.integration,feature = "genotype.classification")
p+scale_color_manual(values=c("green","blue","grey60"))

plot_forceDirectedGraph_label_by_a_feature_of_interest(preleukemic.integration,feature = "donor")

######################################
save(preleukemic.integration,file="Robjects/preleukemic.integration.revised.rdata")
######################################
