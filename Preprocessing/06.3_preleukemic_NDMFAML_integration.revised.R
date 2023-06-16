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

# cell.anno_CD38neg<-cell.anno_CD38neg[,c("cell.id","cell.type","donor","donor.type","batch","genotype_class","genotype_class1")]
# preleukemic_MPNAML<-preleukemic_MPNAML[,c("cell_id","cell_type","donor","genotype.corrected")]
# colnames(preleukemic_MPNAML)<-c("cell.id","cell.type","donor","genotype.corrected")
# preleukemic_MPNAML$donor.type<-"MPNAML"
# preleukemic_MPNAML$donor.type[preleukemic_MPNAML$donor=="HD15_8650"]<-"normal"
# preleukemic_MPNAML$batch<-"MPNAML"
# 
# preleukemic_MPNAML$genotype_class1<-NA
# preleukemic_MPNAML$genotype_class1[preleukemic_MPNAML$genotype.corrected=="CALR_Het_TP53_278_Het_TP53_162_Het"]<-"leukemic"
# preleukemic_MPNAML$genotype_class1[preleukemic_MPNAML$genotype.corrected=="CALR_Het_TP53_278_WT_TP53_162_WT"]<-"CALR-MPNAML"
# preleukemic_MPNAML$genotype_class1[preleukemic_MPNAML$genotype.corrected=="JAK2_HET_TP53_175_HOM_TET2_1863_HOM"]<-"leukemic"
# preleukemic_MPNAML$genotype_class1[preleukemic_MPNAML$genotype.corrected=="JAK2_HET_TP53_175_WT_TET2_1863_WT"]<-"JAK2-MPNAML"
# preleukemic_MPNAML$genotype_class1[preleukemic_MPNAML$genotype.corrected=="JAK2_HET_U2AF1_HET_TP53_p272_HET_TP53_c672_HET"]<-"leukemic"
# preleukemic_MPNAML$genotype_class1[preleukemic_MPNAML$genotype.corrected=="JAK2_HET_U2AF1_HET_TP53_p272_WT_TP53_c672_WT"]<-"JAK2-MPNAML"
# preleukemic_MPNAML$genotype_class1[preleukemic_MPNAML$genotype.corrected=="JAK2_WT_TP53_175_WT_TET2_1863_WT"]<-"WT-MPNAML"
# preleukemic_MPNAML$genotype_class1[preleukemic_MPNAML$genotype.corrected=="JAK2_WT_U2AF1_WT_TP53_p272_WT_TP53_c672_WT"]<-"WT-MPNAML"
# preleukemic_MPNAML$genotype_class1[preleukemic_MPNAML$genotype.corrected=="TP53_p241_HET_TP53_p273_HET"]<-"leukemic"
# preleukemic_MPNAML$genotype_class1[preleukemic_MPNAML$genotype.corrected=="TP53_p241_WT_TP53_p273_WT"]<-"CALR-MPNAML"


# preleukemic_MPNAML$genotype_class<-NA
# preleukemic_MPNAML$genotype_class[preleukemic_MPNAML$genotype_class1=="leukemic"]<-"leukemic"
# preleukemic_MPNAML$genotype_class[preleukemic_MPNAML$genotype_class1=="CALR-MPNAML"]<-"preleukemic"
# preleukemic_MPNAML$genotype_class[preleukemic_MPNAML$genotype_class1=="JAK2-MPNAML"]<-"preleukemic"
# preleukemic_MPNAML$genotype_class[preleukemic_MPNAML$genotype_class1=="WT-MPNAML"]<-"WT"

# preleukemic_MPNAML$genotype.corrected<-NULL

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
