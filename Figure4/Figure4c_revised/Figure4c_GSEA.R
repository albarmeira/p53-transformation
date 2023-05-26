#Differential gene expression analysis of p53 heterozygous cells 
#from pre-TP53-sAML and CP-TP53-MPN patients

###################################
library(SingleCellExperiment)
library(Matrix)
library(matrixStats)
library(bigmemory)
library(SingCellaR)
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
#Extract counts and normalize data
################################################
##Extract counts from SingCellaR integration object
#This object contains normalized counts from all patients and normal donors included in the study,
#inclusive of Rodriguez-Meira et al, 2019 dataset from myelofibrosis patients

umi.dat<-get_normalized_umi(MPN.integration.all)

#Extracting gene information
genes.info<-get_genes_metadata(MPN.integration.all)

# Subset the counts matrix by expressed genes
genes.info<-subset(genes.info,IsExpress==TRUE)
gene.index<-rownames(umi.dat) %in% rownames(genes.info)
umi.dat<-umi.dat[gene.index,]

# Normalized counts using log2 transformation prior to DGE analysis
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
cellsA.m<-umi.dat.log2[,a.index] #24207 genes 296 cells, preTP53sAML_p53HET group
cellsB.m<-umi.dat.log2[,b.index] #24207 genes 273 cells, CP_TP53_MPN_p53HET group

##################################################################
#### Prepare tables for GSEA analysis ############################
##################################################################
groupA<-as.character(colnames(cellsA.m)) #preTP53sAML_p53HET group
groupB<-as.character(colnames(cellsB.m)) #CP_TP53_MPN_p53HET group

groupA.m<-as.matrix(cellsA.m)
groupB.m<-as.matrix(cellsB.m)
#############################################################

my.phenotype.s1<-rep("preTP53sAML_p53HET",ncol(groupA.m))
my.phenotype.s2<-rep("CP_TP53_MPN_p53HET",ncol(groupB.m))

my.phenotype<-c(my.phenotype.s1,my.phenotype.s2)

my.data<-cbind(groupA.m,groupB.m)

my.info<-data.frame(NAME=rownames(my.data))
my.info$DESCRIPTION<-"na"

my.final<-cbind(my.info,my.data)

h1<-paste(ncol(my.data),"2","1",sep=" ")
h2<-paste("#","preTP53sAML_p53HET","CP_TP53_MPN_p53HET",sep=" ")
h3<-paste(c(rep("preTP53sAML_p53HET",length(groupA)),rep("CP_TP53_MPN_p53HET",length(groupB)),sep=" "))

cat(h1,file=paste("GSEA/preTP53sAML_vs_CP_TP53_MPN_p53HET.logged_GSEA-phenotype",".cls",sep=""),sep="\n")
cat(h2,file=paste("GSEA/preTP53sAML_vs_CP_TP53_MPN_p53HET.logged_GSEA-phenotype",".cls",sep=""),sep="\n",append=TRUE)
cat(h3,file=paste("GSEA/preTP53sAML_vs_CP_TP53_MPN_p53HET.logged_GSEA-phenotype",".cls",sep=""),append=TRUE)

write.table(my.final,file="GSEA/preTP53sAML_vs_CP_TP53_MPN_p53HET.logged_GSEA.txt",
            append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)

#############################################################