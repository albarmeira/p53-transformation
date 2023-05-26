################################################
library(SingCellaR)
library(ggplot2)
load(file="../../../Sources/HTMPNAML_final.revised.rdata")
######################################

# Select patient IF0131/SB5702
HTMPNAML@meta.data$IsPassed<-FALSE
HTMPNAML@meta.data$IsPassed[HTMPNAML@meta.data$stage=="IF0131"]<-TRUE
HTMPNAML@meta.data$IsPassed[HTMPNAML@meta.data$stage=="SB5702"]<-TRUE
table(HTMPNAML@meta.data$IsPassed) #1081

#############################################################################
normalize_UMIs(HTMPNAML,use.scaled.factor = F)
remove_unwanted_confounders(HTMPNAML,residualModelFormulaStr="~detectedGenesPerCell+UMI_count") 
#############################################################################

get_variable_genes_by_fitting_GLM_model(HTMPNAML,mean_expr_cutoff = 1,disp_zscore_cutoff = 0.1,
                                        quantile_genes_expr_for_fitting = 0.6,quantile_genes_cv2_for_fitting = 0.3)
#"Identified :2945 variable genes"

remove_unwanted_genes_from_variable_gene_set(HTMPNAML,gmt.file = "../../../Sources/MC.human.signature.genes.gmt",
                                             removed_gene_sets=c("Ribosomal_gene","Mitocondrial_gene"))

plot_variable_genes(HTMPNAML,quantile_genes_expr_for_fitting = 0.6,quantile_genes_cv2_for_fitting = 0.3)
#####################################
runPCA(HTMPNAML,use.components=30,use.regressout.data = T)
plot_PCA_Elbowplot(HTMPNAML)
#####################################
#Add CNV/genotype clusters to metadata
IF0131.CNV<-read.table(file="cnv_summary_IF0131_subclusters_genotype.txt",sep="\t")
IF0131.CNV<-subset(IF0131.CNV,genotype=="TP53_multihit_M2")
IF0131.CNV$cell.id<-rownames(IF0131.CNV)
dim(IF0131.CNV) #726

HTMPNAML@meta.data<-merge(HTMPNAML@meta.data,IF0131.CNV,by.x="Cell",by.y="cell.id")
#####################################
runUMAP(HTMPNAML,n.dims.use=10,n.neighbors=20,uwot.min.dist=0.30,uwot.metric = "correlation",n.seed = 1,
        useIntegrativeEmbeddings=FALSE,
        dim_reduction_method="pca")
plot_umap_label_by_a_feature_of_interest(HTMPNAML,feature = "cluster_corrected")
######################################
res.umap<-HTMPNAML@umap.result
res.umap<-subset(res.umap,cluster_corrected!="NA")
point.size<-2
ggplot(res.umap,aes(x=UMAP1,y=-UMAP2,colour=cluster_corrected))+geom_point(size=point.size)+theme_void()+
  scale_colour_manual(values=c("deeppink4","violetred2","cyan3"))+
  theme(legend.position="none")
ggsave("FigureS3m.png",width = 5,height = 5)
