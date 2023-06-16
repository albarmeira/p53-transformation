#Merge counts tables and colData from all batches of sequencing; convert patient IDs to paper IDs
#Data preprocessing and QC metrics
#Author:Alba Rodriguez-Meira
#Date: 6th March 2021
#Modified last: 26th February 2023
####################################################

#Import metadata tables from QC'ed cells

colDataB1.qc<-read.table("../02-MPNAMLB1/p53_HTMPNAMLB1_colData.qc.csv",header = T,sep=",") #Read metadata table with QC'ed cells
colDataB2.qc<-read.table("../03-MPNAMLB2/p53_HTMPNAMLB2_colData.qc.csv",header = T,sep=",") #Read metadata table with QC'ed cells
colDataB3.qc<-read.table("../04-MPNAMLB3/p53_HTMPNAMLB3_colData.qc.csv",header = T,sep=",") #Read metadata table with QC'ed cells

colData.qc.final<-rbind(colDataB1.qc,colDataB2.qc,colDataB3.qc)

colnames(colData.qc.final)
dim(colData.qc.final) #19221 cells in 23 columns; QC column=23; previously 19383; 162 QC'ed cells removed from the dataset

table(colData.qc.final$donor)
## Convert donor IDs to new paper IDs
colData.qc.final$donor[colData.qc.final$donor=="JB4211"]<-"IF0308"
colData.qc.final$donor[colData.qc.final$donor=="GST010"]<-"IF0308"
colData.qc.final$donor[colData.qc.final$donor=="VO4409"]<-"IF0318"
colData.qc.final$donor[colData.qc.final$donor=="PM4908"]<-"IF0391"
colData.qc.final$donor[colData.qc.final$donor=="DT5210"]<-"IF0392"
colData.qc.final$donor[colData.qc.final$donor=="KEFA"]<-"GR001"
colData.qc.final$donor[colData.qc.final$donor=="SAJU"]<-"GR002"
colData.qc.final$donor[colData.qc.final$donor=="JH4406"]<-"IF0393"
colData.qc.final$donor[colData.qc.final$donor=="BEMA"]<-"GR003"
colData.qc.final$donor[colData.qc.final$donor=="GUMI"]<-"GR004"
colData.qc.final$donor[colData.qc.final$donor=="LAJA"]<-"GR005"
colData.qc.final$donor[colData.qc.final$donor=="MUGU1"]<-"GR006"
colData.qc.final$donor[colData.qc.final$donor=="LERO"]<-"GR007"
colData.qc.final$donor[colData.qc.final$donor=="AHGH"]<-"GH001"
colData.qc.final$donor[colData.qc.final$donor=="MAJIC38"]<-"IF9038"
colData.qc.final$donor[colData.qc.final$donor=="MAJIC118"]<-"IF9118"
colData.qc.final$donor[colData.qc.final$donor=="MAJIC131"]<-"IF9131"
colData.qc.final$donor[colData.qc.final$donor=="MAJIC180"]<-"IF9180"

table(colData.qc.final$donor)

##
table(colData.qc.final$stage)

colData.qc.final$stage[colData.qc.final$stage=="AHGH_001"]<-"GH001_001"
colData.qc.final$stage[colData.qc.final$stage=="AHGH_002"]<-"GH001_002"
colData.qc.final$stage[colData.qc.final$stage=="AHGH_003"]<-"GH001_003"
colData.qc.final$stage[colData.qc.final$stage=="BEMA_AP"]<-"GR003_AP"
colData.qc.final$stage[colData.qc.final$stage=="BEMA_CP"]<-"GR003_CP"
colData.qc.final$stage[colData.qc.final$stage=="DT5210"]<-"IF0392"
colData.qc.final$stage[colData.qc.final$stage=="GUMI_AP"]<-"GR004_AP"
colData.qc.final$stage[colData.qc.final$stage=="GUMI_CP"]<-"GR004_CP"
colData.qc.final$stage[colData.qc.final$stage=="JH4406"]<-"IF0393"
colData.qc.final$stage[colData.qc.final$stage=="KEFA"]<-"GR001"
colData.qc.final$stage[colData.qc.final$stage=="LAJA_AP"]<-"GR005_AP"
colData.qc.final$stage[colData.qc.final$stage=="LAJA_CP"]<-"GR005_CP"
colData.qc.final$stage[colData.qc.final$stage=="LERO_AP"]<-"GR007_AP"
colData.qc.final$stage[colData.qc.final$stage=="LERO_CP"]<-"GR007_CP"
colData.qc.final$stage[colData.qc.final$stage=="MAJIC118_05HS31"]<-"IF9118_05HS31"
colData.qc.final$stage[colData.qc.final$stage=="MAJIC118_0GX632"]<-"IF9118_0GX632"
colData.qc.final$stage[colData.qc.final$stage=="MAJIC131"]<-"IF9131"
colData.qc.final$stage[colData.qc.final$stage=="MAJIC180"]<-"IF9180"
colData.qc.final$stage[colData.qc.final$stage=="MAJIC38"]<-"IF9038"
colData.qc.final$stage[colData.qc.final$stage=="MUGU1_AP"]<-"GR006_AP"
colData.qc.final$stage[colData.qc.final$stage=="MUGU1_CP"]<-"GR006_CP"
colData.qc.final$stage[colData.qc.final$stage=="PM4908"]<-"IF0391"
colData.qc.final$stage[colData.qc.final$stage=="SAJU"]<-"GR002"
colData.qc.final$stage[colData.qc.final$stage=="VO4409"]<-"IF0318"

table(colData.qc.final$stage)

##
write.table(colData.qc.final,file="p53_HTMPNAML_colData.qc.final.revised.new.id.csv",sep=",",row.names = F)

#Import counts tables from QC'ed cells

counts.B1.qc<-read.table("../02-MPNAMLB1/p53_HTMPNAMLB1_counts.qc.tsv",header = T,sep="\t") #Read metadata table with QC'ed cells
#32830 genes 7123 cells
counts.B2.qc<-read.table("../03-MPNAMLB2/p53_HTMPNAMLB2_counts.qc.tsv",header = T,sep="\t") #Read metadata table with QC'ed cells
#32830 genes 5779 cells
counts.B2.qc[1]<-NULL #Remove Gene column
counts.B3.qc<-read.table("../04-MPNAMLB3/p53_HTMPNAMLB3_counts.qc.tsv",header = T,sep="\t") #Read metadata table with QC'ed cells
#32830 genes 6319 cells; previously 6481 cells
counts.B3.qc[1]<-NULL #Remove Gene column

counts.qc.final<-cbind(counts.B1.qc,counts.B2.qc,counts.B3.qc) #32830 genes 19221 cells
dim(counts.qc.final)

write.table(counts.qc.final,file="p53_HTMPNAML_counts.qc.all.revised.tsv",sep="\t",row.names = F)




