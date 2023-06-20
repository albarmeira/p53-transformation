###################################
# Integration of CD38neg cells from MolCell2019 dataset and p53MPNAML dataset
# Selected donors
###################################
load(file="../classification/preleukemic_integration/counts.qc_CD38neg_TARGET.revised.rdata")
load(file="../classification/preleukemic_integration/cell.anno_CD38neg_TARGET.revised.rdata")

load(file="../classification/preleukemic_integration/counts.qc_CD38neg_p53.revised.rdata")
load(file="../classification/preleukemic_integration/cell.anno_CD38neg_p53.revised.rdata")
###################################

#Merge metadata
table(is.element(colnames(cell.anno_CD38neg_TARGET),colnames(cell.anno_CD38neg_p53)))
table(is.element(colnames(cell.anno_CD38neg_p53),colnames(cell.anno_CD38neg_TARGET)))

cell.anno_CD38neg_p53$sort<-NULL

###################################
## Merge metadata datasets

metadata_combined<-rbind(cell.anno_CD38neg_TARGET,cell.anno_CD38neg_p53) #Column of QC=23
write.table(metadata_combined,file="../classification/preleukemic_integration/preleukemic_integration_metadata.qc.revised.csv",sep=",",row.names = F)

## Merge count tables
merged_counts<-cbind(counts.qc_CD38neg_TARGET,counts.qc_CD38neg_p53)
genes<-read.table("../genes_withERCC_uniqueids.tsv",header = F,sep = "\t")
Gene<-as.matrix(genes$V3)
counts_preleukemic.integration<-cbind(Gene,merged_counts)
colnames(counts_preleukemic.integration)[1] <- "Gene"
write.table(counts_preleukemic.integration,file="../classification/preleukemic_integration/counts_preleukemic.integration.revised.tsv",sep="\t",row.names = F)
###################################
