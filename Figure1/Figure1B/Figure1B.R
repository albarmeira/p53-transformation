load(file="metadata_IF0131.qc.rdata")

table(metadata_IF0131.qc$cell_type)

#Exclude sorted HSCs from the count
subset<-c("CD34+CD38-","CD34+CD38+")
metadata_IF0131.qc<-metadata_IF0131.qc[metadata_IF0131.qc$cell_type %in% subset,]

#Percentage of CD34+CD38- sorted cells (considering sorting gates)
#13.59 %
#Percentage of CD34+CD38+ sorted cells (considering sorting gates)
#83.01 %

accepted_genotypes<-c("JAK2_WT_U2AF1_WT_TP53_c672_WT_TP53_p272_WT",
                      "JAK2_WT_U2AF1_HET_TP53_c672_WT_TP53_p272_WT",
                      "JAK2_HET_U2AF1_HET_TP53_c672_WT_TP53_p272_WT",
                      "JAK2_HET_U2AF1_HET_TP53_c672_HET_TP53_p272_WT",
                      "JAK2_HET_U2AF1_HET_TP53_c672_HET_TP53_p272_HET",
                      "JAK2_HOM_U2AF1_HET_TP53_c672_HET_TP53_p272_HET")

metadata_IF0131.qc.accepted.genotypes<-metadata_IF0131.qc[metadata_IF0131.qc$genotype %in% accepted_genotypes,]

metadata_IF0131.qc.CD38neg<-subset(metadata_IF0131.qc.accepted.genotypes,cell_type=="CD34+CD38-")
metadata_IF0131.qc.CD38pos<-subset(metadata_IF0131.qc.accepted.genotypes,cell_type=="CD34+CD38+")
  
table(metadata_IF0131.qc.CD38neg$genotype)
table(metadata_IF0131.qc.CD38pos$genotype)

f_CD38neg<-13.59
f_CD38pos<-83.01
  
#Multiply CD34+CD38- cells by the CD34+CD38- cell ratio 13.59/(83.01+13.59)=0.1406832
metadata_IF0131.qc.CD38neg_df<-as.data.frame(table(metadata_IF0131.qc.CD38neg$genotype))
metadata_IF0131.qc.CD38neg_df$Freq<-metadata_IF0131.qc.CD38neg_df$Freq*(f_CD38neg/(f_CD38neg+f_CD38pos))

#Multiply CD34+CD38+ cells by the CD34+CD38+ cell ratio 83.01/(83.01+13.59)=0.8593168
metadata_IF0131.qc.CD38pos_df<-as.data.frame(table(metadata_IF0131.qc.CD38pos$genotype))
metadata_IF0131.qc.CD38pos_df$Freq<-metadata_IF0131.qc.CD38pos_df$Freq*(f_CD38pos/(f_CD38neg+f_CD38pos))

metadata_IF0131.merged<-merge(metadata_IF0131.qc.CD38neg_df,metadata_IF0131.qc.CD38pos_df,by="Var1")

metadata_IF0131.merged$Freq.total<-metadata_IF0131.merged$Freq.x+metadata_IF0131.merged$Freq.y
metadata_IF0131.merged$fraction<-metadata_IF0131.merged$Freq.total/sum(metadata_IF0131.merged$Freq.total)

metadata_IF0131.merged

