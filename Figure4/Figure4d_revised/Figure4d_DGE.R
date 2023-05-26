#Differential gene expression analysis of p53 heterozygous cells 
#from pre-TP53-sAML and CP-TP53-MPN patients

###################################
library(SingleCellExperiment)
library(SingCellaR)
#library(Matrix)
#library(matrixStats)
library(bigmemory)
library(ggplot2)
##################################
source("../../../Sources/TARGETseq_functions.R")
##################################
#Load counts object
load(file="../../../Sources/MFp53MPNAML_integration.revised.rdata")

#Load metadata and subset p53-heterozygous cells from each patient group
metadata<-read.table(file="../../../Sources/metadata_MPNAMLp53_with_index_genotype.revised.txt",sep="\t",header = T)

preTP53sAML_donors<-c("GH001_001","GH001_002","GR003_CP","GR005_CP","GR007_CP","GR006_CP")
CP_TP53_MPN_donors<-c("IF9118_05HS31","IF9118_0GX632","IF9131","IF9180","IF9038")

preTP53sAML_p53HET<-subset(metadata,stage %in% preTP53sAML_donors)
preTP53sAML_p53HET<-subset(preTP53sAML_p53HET,genotype.classification=="TP53HET")
dim(preTP53sAML_p53HET)
table(preTP53sAML_p53HET$stage)

CP_TP53_MPN_p53HET<-subset(metadata,stage %in% CP_TP53_MPN_donors)
CP_TP53_MPN_p53HET<-subset(CP_TP53_MPN_p53HET,genotype.classification=="TP53HET")
dim(CP_TP53_MPN_p53HET)
table(CP_TP53_MPN_p53HET$stage)

################################################
#Perform DE analysis
################################################
#Extract counts from object
# Obtained normalized counts from the R object MPN.integration.all (which contains
#gene counts from MF and MPNAML patients as well as normal donor controls)
umi.dat<-get_normalized_umi(MPN.integration.all)
#Extracting gene information
genes.info<-get_genes_metadata(MPN.integration.all)
# Subset the counts matrix by expressed genes
genes.info<-subset(genes.info,IsExpress==TRUE)
gene.index<-rownames(umi.dat) %in% rownames(genes.info)
umi.dat<-umi.dat[gene.index,]

# Normalized counts using log2 transformation before doing DGE analysis
umi.dat.log2<-umi.dat
min_expr<-1 
umi.dat.log2[umi.dat.log2 < min_expr]<- 0
umi.dat.log2[umi.dat.log2 >= min_expr]<-log2(umi.dat.log2[umi.dat.log2 >= min_expr])

# Subset names from each genotype group
cols.A<-preTP53sAML_p53HET$cell_id
cols.B<-CP_TP53_MPN_p53HET$cell_id

a.index<-which(colnames(umi.dat.log2) %in% cols.A)
b.index<-which(colnames(umi.dat.log2) %in% cols.B)

# Subset matrix including only cells from selected genotype groups
cellsA.m<-umi.dat.log2[,a.index] #24207 genes 296 cells
cellsB.m<-umi.dat.log2[,b.index] #24207 genes 273 cells

# Now perform differential expression analysis between groups
p53HET_DGE<-dif_expression(groupA="pre-TP53-sAML",groupB="CP-TP53-MPN",groupA.m=cellsA.m,groupB.m=cellsB.m,pvalue_cutoff=0.1)

p53HET_DGE$class<-"NA"
p53HET_DGE$class[p53HET_DGE$log2fc<(-1) & p53HET_DGE$p.adjusted<0.1]<-"up_in_pre-TP53-sAML"
p53HET_DGE$class[p53HET_DGE$log2fc>1 & p53HET_DGE$p.adjusted<0.1]<-"up_in_CP-TP53-MPN"

genes<-c("IFITM1","BANF1","OAS1","IFITM3","BST2","IFITM2","IFI16")
p53HET_DGE$Name<-p53HET_DGE$gene
p53HET_DGE$Name[!(p53HET_DGE$gene %in% genes)]<-NA

library(ggrepel)

ggplot(p53HET_DGE,aes(x=-log2fc,y=-log10(p.adjusted),colour=class))+geom_point()+
  theme_classic()+
  geom_vline(xintercept = -1,colour="grey50")+
  geom_vline(xintercept = 1,colour="grey50")+
  scale_x_continuous(limits = c(-5,5))+
  scale_y_continuous(limits = c(0,100))+
  geom_text_repel(aes(label=Name),hjust=-1.5, vjust=0.3, size=4,colour="black")+
  scale_colour_manual(values = c("grey","forestgreen","darkorange"))+
  theme(axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15))

ggsave(filename="Figure4d_volcano.png",width = 6,height = 5)
write.table(p53HET_DGE,file="p53HET_DGE.txt",sep="\t",row.names = F)
