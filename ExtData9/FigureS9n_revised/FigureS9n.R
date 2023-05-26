###################################
### Violin plots from IFN-response genes upregulated in p53HET cells from tCP as compared to MAJIC patients
### Author: Alba Rodriguez-Meira
### Date created: 22nd October 2021
### Date modified: 11th March 2023
###################################
library(SingleCellExperiment)
library(SingCellaR)
library(bigmemory)
library(ggplot2)
##################################
source("../../../Sources/TARGETseq_functions.R")
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

IFNGR1<-as.data.frame(subset(umi.dat.log2,rownames(umi.dat.log2)=="IFNGR1"))
rownames(IFNGR1)<-colnames(umi.dat.log2)
IFNGR1$group<-"preTP53sAML"
colnames(IFNGR1)<-c("value","group")
names.C<-c(as.character(cols.A),as.character(cols.B))
IFNGR1<-IFNGR1[rownames(IFNGR1) %in% names.C,]
IFNGR1$group[rownames(IFNGR1) %in% cols.B]<-"CP_TP53_MPN"


ggplot(IFNGR1,aes(x=group,y=value,fill=group))+
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

ggsave(filename="FigureS9n_IFNGR1_violin.png",width = 4,height = 5)
################
IFNGR2<-as.data.frame(subset(umi.dat.log2,rownames(umi.dat.log2)=="IFNGR2"))
rownames(IFNGR2)<-colnames(umi.dat.log2)
IFNGR2$group<-"preTP53sAML"
colnames(IFNGR2)<-c("value","group")
names.C<-c(as.character(cols.A),as.character(cols.B))
IFNGR2<-IFNGR2[rownames(IFNGR2) %in% names.C,]
IFNGR2$group[rownames(IFNGR2) %in% cols.B]<-"CP_TP53_MPN"

ggplot(IFNGR2,aes(x=group,y=value,fill=group))+
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

ggsave(filename="FigureS9n_IFNGR2_violin.png",width = 4,height = 5)
################
IFNAR1<-as.data.frame(subset(umi.dat.log2,rownames(umi.dat.log2)=="IFNAR1"))
rownames(IFNAR1)<-colnames(umi.dat.log2)
IFNAR1$group<-"preTP53sAML"
colnames(IFNAR1)<-c("value","group")
names.C<-c(as.character(cols.A),as.character(cols.B))
IFNAR1<-IFNAR1[rownames(IFNAR1) %in% names.C,]
IFNAR1$group[rownames(IFNAR1) %in% cols.B]<-"CP_TP53_MPN"

ggplot(IFNAR1,aes(x=group,y=value,fill=group))+
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

ggsave(filename="FigureS9n_IFNAR1_violin.png",width = 4,height = 5)
#########################
IFNAR2<-as.data.frame(subset(umi.dat.log2,rownames(umi.dat.log2)=="IFNAR2"))
rownames(IFNAR2)<-colnames(umi.dat.log2)
IFNAR2$group<-"preTP53sAML"
colnames(IFNAR2)<-c("value","group")
names.C<-c(as.character(cols.A),as.character(cols.B))
IFNAR2<-IFNAR2[rownames(IFNAR2) %in% names.C,]
IFNAR2$group[rownames(IFNAR2) %in% cols.B]<-"CP_TP53_MPN"

ggplot(IFNAR2,aes(x=group,y=value,fill=group))+
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

ggsave(filename="FigureS9n_IFNAR2_violin.png",width = 4,height = 5)
#########################
