##################################
#DGE and GSEA analysis from cytogenetically distinct subclones in patient IF0131
#Author: Alba Rodriguez-Meira, DPhil
#Date created: 14th October 2021
#Date last modified: 14th October 2021
##################################
source("../../../Sources/TARGETseq_functions.R")
##################################
#Import metadata, counts tables and gene annotations
metadata<-read.table(file="input/cnv_summary_IF0131_subclusters_genotype.txt",header = T,sep="\t")
#counts<-read.table(file="../Sources/p53_HTMPNAML_counts.qc.all.tsv",sep="\t",header = T)
#genes<-read.table(file="../Sources/genes_withERCC_uniqueids.tsv",sep="\t",header = F)
#rownames(counts)<-genes$V3
##Subset counts to include only IF0131-p53 CNV-containing cells
#counts.subset<-counts[,colnames(counts) %in% metadata$cell_id]
#write.table(counts.subset,file="input/counts.subset_IF0131_CNV.txt",sep="\t")
counts.subset<-read.table(file="input/counts.subset_IF0131_CNV.txt",sep="\t",row.names = 1)

### Perform counts normalization
normalize_counts <- function(your.matrix,scale.factor = 1e4,use.scaled.factor=T){
  
  umi.dat<-your.matrix
  if(use.scaled.factor==T){
    totalUMI.lib<-Matrix::colSums(umi.dat)
    normalized.umi<-(t(t(umi.dat)/totalUMI.lib))*scale.factor
    return(normalized.umi)
  }else{
    totalUMI.lib<-Matrix::colSums(umi.dat)
    normalized.umi<-(t(t(umi.dat)/totalUMI.lib))*round(mean(totalUMI.lib))
    return(normalized.umi)
  }
  
}

counts.normalized<-normalize_counts(counts.subset,use.scaled.factor=F)
counts.normalized.IF0131.CNV<-counts.normalized[,as.character(metadata$cell_id)]

#Remove ERCC from analysis
cell.count.no.ERCC <-counts.normalized.IF0131.CNV[-c(grep("ERCC-", rownames(counts.normalized.IF0131.CNV))),]
###################################
#Remove genes expressed in less than 5 cells
genes_with_expressing_cells<-5
cells.per.genes<-as.data.frame(apply(cell.count.no.ERCC>1,1,sum))
colnames(cells.per.genes)<-c("cell.exp.genes")
genes<-rownames(cells.per.genes[cells.per.genes$cell.exp.genes>5,,drop=FALSE])
cell.count.use<-as.matrix(cell.count.no.ERCC[genes,])
#############################################################
#Perform log2 transformation
cell.count.use.log2<-cell.count.use
min_expr<-1
cell.count.use.log2[cell.count.use.log2 < min_expr]<- 0
cell.count.use.log2[cell.count.use.log2 >= min_expr]<-log2(cell.count.use.log2[cell.count.use.log2 >= min_expr])
#############################################################

####### DE analysis of chr3_chr5 TP53-leukemic cells versus chr3_chr5_chr7 subclone

cloneA<-metadata[metadata$cluster_corrected=="chr3_chr5",c("cell_id")]
cloneB<-metadata[metadata$cluster_corrected=="chr3_chr5_chr7",c("cell_id")]

a.index<-which(colnames(cell.count.use.log2) %in% cloneA)
b.index<-which(colnames(cell.count.use.log2) %in% cloneB)

cellsA.m<-cell.count.use.log2[,a.index]
cellsB.m<-cell.count.use.log2[,b.index]

cloneAvscloneB<-dif_expression(groupA="cloneA",groupB="cloneB_del7",groupA.m=cellsA.m,groupB.m=cellsB.m,pvalue_cutoff=0.1)

cloneAvscloneB$class<-"NA"
cloneAvscloneB$class[cloneAvscloneB$log2fc<(-0.5) & cloneAvscloneB$p.adjusted<0.1]<-"up_in_cloneA"
cloneAvscloneB$class[cloneAvscloneB$log2fc>0.5 & cloneAvscloneB$p.adjusted<0.1]<-"up_in_cloneBdel7"

cloneAvscloneB$Name<-cloneAvscloneB$gene
cloneAvscloneB$Name[-log10(cloneAvscloneB$p.adjusted)<10]<-NA

library(ggplot2)
ggplot(cloneAvscloneB,aes(x=log2fc,y=-log10(p.adjusted),colour=class))+geom_point()+
  theme_bw()+
  geom_vline(xintercept = -0.5)+
  geom_vline(xintercept = 0.5)+
  scale_x_continuous(limits = c(-3,3))+
  scale_y_continuous(limits = c(0,50))+
  geom_text(aes(label=Name),hjust=1.2, vjust=0.5, size=3)+
  scale_colour_manual(values = c("grey","green","red"))

write.table(cloneAvscloneB,file="output/DGE_cloneAvscloneB.txt",col.names = TRUE,sep = "\t")

###############################################################################
#### Prepare tables for GSEA analysis #########################################
###############################################################################
my.phenotype.s1<-rep("cloneA",ncol(cellsA.m))
my.phenotype.s2<-rep("cloneB",ncol(cellsB.m))

my.phenotype<-c(my.phenotype.s1,my.phenotype.s2)

my.data<-cbind(cellsA.m,cellsB.m)

my.info<-data.frame(NAME=rownames(my.data))
my.info$DESCRIPTION<-"na"

my.final<-cbind(my.info,my.data)

h1<-paste(ncol(my.data),"2","1",sep=" ")
h2<-paste("#","cloneA","cloneB",sep=" ")
h3<-paste(c(rep("cloneA",length(cloneA)),rep("cloneB",length(cloneB)),sep=" "))

cat(h1,file=paste("output/cloneAvscloneB.logged_GSEA-phenotype",".cls",sep=""),sep="\n")
cat(h2,file=paste("output/cloneAvscloneB.logged_GSEA-phenotype",".cls",sep=""),sep="\n",append=TRUE)
cat(h3,file=paste("output/cloneAvscloneB.logged_GSEA-phenotype",".cls",sep=""),append=TRUE)

write.table(my.final,file="output/cloneAvscloneB.logged_GSEA.txt",
            append=FALSE, sep="\t", quote=FALSE,row.names=FALSE, col.names=TRUE)
