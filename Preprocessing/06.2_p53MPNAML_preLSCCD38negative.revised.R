###############################################################################
library(Matrix)
library(matrixStats)
##################################
colData.qc<-read.table(file="p53_HTMPNAML_colData.qc.final.revised.new.id.csv",sep=",",header = T)

counts.p53MPNAML<-read.table(file="p53_HTMPNAML_counts.qc.all.revised.new.id.tsv",sep="\t",header = T)
#rownames(counts.p53MPNAML) <- counts.p53MPNAML$Gene
###################################
#cell.anno<-colData.qc[colData.qc$cell_id %in% colnames(counts.p53MPNAML),]

#Subset metadata for CD38 negative cells and preLSCs (defined genetically)
colData.qc$sort<-sapply(strsplit(as.character(colData.qc$cell_id),"PL"), `[`, 1)

#Each sort has a different CD38 threshold
CD38_MPNAML08<-8000
CD38_MPNAML10<-6000
CD38_MPNAML13<-3500
CD38_MPNAML16<-5000
CD38_MPNAML21<-1000
CD38_MPNAML22<-1400
CD38_MPNAML23<-1500
CD38_MPNAML24<-370
CD38_MPNAML29<-1500
CD38_MPNAML30<-1500
CD38_MPNAML31<-2000

#no selection for CD90 expression
cell.anno_CD38neg_MPNAML08<-subset(colData.qc,CD38<CD38_MPNAML08 & colData.qc$sort=="MPNAML08")
cell.anno_CD38neg_MPNAML10<-subset(colData.qc,CD38<CD38_MPNAML10 & colData.qc$sort=="MPNAML10")
cell.anno_CD38neg_MPNAML13<-subset(colData.qc,CD38<CD38_MPNAML13 & colData.qc$sort=="MPNAML13")
cell.anno_CD38neg_MPNAML16<-subset(colData.qc,CD38<CD38_MPNAML16 & colData.qc$sort=="MPNAML16")
cell.anno_CD38neg_MPNAML21<-subset(colData.qc,CD38<CD38_MPNAML21 & colData.qc$sort=="MPNAML21")
cell.anno_CD38neg_MPNAML22<-subset(colData.qc,CD38<CD38_MPNAML22 & colData.qc$sort=="MPNAML22")
cell.anno_CD38neg_MPNAML23<-subset(colData.qc,CD38<CD38_MPNAML23 & colData.qc$sort=="MPNAML23")
cell.anno_CD38neg_MPNAML24<-subset(colData.qc,CD38<CD38_MPNAML24 & colData.qc$sort=="MPNAML24")
cell.anno_CD38neg_MPNAML29<-subset(colData.qc,CD38<CD38_MPNAML29 & colData.qc$sort=="MPNAML29")
cell.anno_CD38neg_MPNAML30<-subset(colData.qc,CD38<CD38_MPNAML30 & colData.qc$sort=="MPNAML30")
cell.anno_CD38neg_MPNAML31<-subset(colData.qc,CD38<CD38_MPNAML31 & colData.qc$sort=="MPNAML31")

cell.anno_CD38neg_p53<-rbind(cell.anno_CD38neg_MPNAML08,cell.anno_CD38neg_MPNAML10,cell.anno_CD38neg_MPNAML13,
                         cell.anno_CD38neg_MPNAML16,cell.anno_CD38neg_MPNAML21,cell.anno_CD38neg_MPNAML22,
                         cell.anno_CD38neg_MPNAML23,cell.anno_CD38neg_MPNAML24,cell.anno_CD38neg_MPNAML29,
                         cell.anno_CD38neg_MPNAML30,cell.anno_CD38neg_MPNAML31)
#check how many CD38+ cells are still included
table(cell.anno_CD38neg_p53$cell_type)

#Select specific genotypes and donors
table(cell.anno_CD38neg_p53$genotype.classification)
table(cell.anno_CD38neg_p53$stage) 

#Only normal donor HD15_8650 with >100 cells will be selected to minimize batch variation of the rest of donors
#Only AP donors with more than 20 preleukemic cells will be selected

donors<-c("GR003_AP","GR004_AP","IF0131","GR001","GR005_AP","GR007_AP","GR006_AP","IF0391","GR002","IF0318","SB5702","HD15_8650")
cell.anno_CD38neg_p53<-subset(cell.anno_CD38neg_p53,stage %in% donors)

cell.anno_CD38neg_p53$genotype.classification<-as.character(cell.anno_CD38neg_p53$genotype.classification)
cell.anno_CD38neg_p53$genotype.classification[cell.anno_CD38neg_p53$stage=="HD15_8650"]<-"WT"
cell.anno_CD38neg_p53<-subset(cell.anno_CD38neg_p53,genotype.classification=="preleukemic"|genotype.classification=="WT")

counts.qc_CD38neg_p53<-counts.p53MPNAML[,as.character(cell.anno_CD38neg_p53$cell_id)] #1285 cells; 880 preLSC and 405 WT from HD15; same are before revised data

####################################

save(counts.qc_CD38neg_p53,file="../classification/preleukemic_integration/counts.qc_CD38neg_p53.revised.rdata")
save(cell.anno_CD38neg_p53,file="../classification/preleukemic_integration/cell.anno_CD38neg_p53.revised.rdata")
###################################

