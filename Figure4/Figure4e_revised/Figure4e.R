###################################
### Violin plots from IFN-response genes upregulated in p53HET cells from tCP as compared to MAJIC patients
### Author: Alba Rodriguez-Meira
### Date created: 22nd October 2021
### Date modified: 11th March 2023
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

CP_TP53_MPN_p53HET<-subset(metadata,stage %in% CP_TP53_MPN_donors)
CP_TP53_MPN_p53HET<-subset(CP_TP53_MPN_p53HET,genotype.classification=="TP53HET")
dim(CP_TP53_MPN_p53HET) 

################################################
#Extract counts from R object
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
#############################################################

IFITM1<-as.data.frame(subset(umi.dat.log2,rownames(umi.dat.log2)=="IFITM1"))
rownames(IFITM1)<-colnames(umi.dat.log2)
IFITM1$group<-"preTP53sAML"
colnames(IFITM1)<-c("value","group")
names.C<-c(as.character(cols.A),as.character(cols.B))
IFITM1<-IFITM1[rownames(IFITM1) %in% names.C,]
IFITM1$group[rownames(IFITM1) %in% cols.B]<-"CP_TP53_MPN"


ggplot(IFITM1,aes(x=group,y=value,fill=group))+
  geom_violin()+
  geom_jitter(height = 0,size=2,width=0.2,colour="grey30",alpha=0.6)+
  theme_classic()+
  scale_fill_manual(values = c("forestgreen","darkorange"))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position= "none")+
  scale_y_continuous(limits = c(0,9))

ggsave(filename="Figure4e_IFITM1_violin.png",width = 3,height = 5)
################
IFITM2<-as.data.frame(subset(umi.dat.log2,rownames(umi.dat.log2)=="IFITM2"))
rownames(IFITM2)<-colnames(umi.dat.log2)
IFITM2$group<-"preTP53sAML"
colnames(IFITM2)<-c("value","group")
names.C<-c(as.character(cols.A),as.character(cols.B))
IFITM2<-IFITM2[rownames(IFITM2) %in% names.C,]
IFITM2$group[rownames(IFITM2) %in% cols.B]<-"CP_TP53_MPN"


ggplot(IFITM2,aes(x=group,y=value,fill=group))+
  geom_violin()+
  geom_jitter(height = 0,size=2,width=0.2,colour="grey30",alpha=0.6)+
  theme_classic()+
  scale_fill_manual(values = c("forestgreen","darkorange"))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position= "none")+
  scale_y_continuous(limits = c(0,9))

ggsave(filename="Figure4e_IFITM2_violin.png",width = 3,height = 5)
################
IFITM3<-as.data.frame(subset(umi.dat.log2,rownames(umi.dat.log2)=="IFITM3"))
rownames(IFITM3)<-colnames(umi.dat.log2)
IFITM3$group<-"preTP53sAML"
colnames(IFITM3)<-c("value","group")
names.C<-c(as.character(cols.A),as.character(cols.B))
IFITM3<-IFITM3[rownames(IFITM3) %in% names.C,]
IFITM3$group[rownames(IFITM3) %in% cols.B]<-"CP_TP53_MPN"


ggplot(IFITM3,aes(x=group,y=value,fill=group))+
  geom_violin()+
  geom_jitter(height = 0,size=2,width=0.2,colour="grey30",alpha=0.6)+
  theme_classic()+
  scale_fill_manual(values = c("forestgreen","darkorange"))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position= "none")+
  scale_y_continuous(limits = c(0,11))

ggsave(filename="Figure4e_IFITM3_violin.png",width = 3,height = 5)
#########################
OAS1<-as.data.frame(subset(umi.dat.log2,rownames(umi.dat.log2)=="OAS1"))
rownames(OAS1)<-colnames(umi.dat.log2)
OAS1$group<-"preTP53sAML"
colnames(OAS1)<-c("value","group")
names.C<-c(as.character(cols.A),as.character(cols.B))
OAS1<-OAS1[rownames(OAS1) %in% names.C,]
OAS1$group[rownames(OAS1) %in% cols.B]<-"CP_TP53_MPN"


ggplot(OAS1,aes(x=group,y=value,fill=group))+
  geom_violin()+
  geom_jitter(height = 0,size=2,width=0.2,colour="grey30",alpha=0.6)+
  theme_classic()+
  scale_fill_manual(values = c("forestgreen","darkorange"))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position= "none")+
  scale_y_continuous(limits = c(0,9))

ggsave(filename="Figure4e_OAS1_violin.png",width = 3,height = 5)
#########################
BST2<-as.data.frame(subset(umi.dat.log2,rownames(umi.dat.log2)=="BST2"))
rownames(BST2)<-colnames(umi.dat.log2)
BST2$group<-"preTP53sAML"
colnames(BST2)<-c("value","group")
names.C<-c(as.character(cols.A),as.character(cols.B))
BST2<-BST2[rownames(BST2) %in% names.C,]
BST2$group[rownames(BST2) %in% cols.B]<-"CP_TP53_MPN"


ggplot(BST2,aes(x=group,y=value,fill=group))+
  geom_violin()+
  geom_jitter(height = 0,size=2,width=0.2,colour="grey30",alpha=0.6)+
  theme_classic()+
  scale_fill_manual(values = c("forestgreen","darkorange"))+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position= "none")+
  scale_y_continuous(limits = c(0,9))

ggsave(filename="Figure4e_BST2_violin.png",width = 3,height = 5)
#########################
