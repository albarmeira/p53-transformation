load(file="metadata_GR001.qc.rdata")

table(metadata_GR001.qc$cell_type)

#Exclude sorted HSCs from the count
subset<-c("CD34+CD38-","CD34+CD38+")
metadata_GR001.qc<-metadata_GR001.qc[metadata_GR001.qc$cell_type %in% subset,]
table(metadata_GR001.qc$genotype)

#Percentage of CD34+CD38- sorted cells (considering sorting gates)
#8.93 %
#Percentage of CD34+CD38+ sorted cells (considering sorting gates)
#80.62 %

accepted_genotypes<-c("CALRtype2_WT_TP53_p278_WT_TP53_pI162_WT",
                      "CALRtype2_HET_TP53_p278_WT_TP53_pI162_WT",
                      "CALRtype2_HET_TP53_p278_WT_TP53_pI162_HET",
                      "CALRtype2_HET_TP53_p278_WT_TP53_pI162_HOM",
                      "CALRtype2_HET_TP53_p278_HET_TP53_pI162_HET")

metadata_GR001.qc.accepted.genotypes<-metadata_GR001.qc[metadata_GR001.qc$genotype %in% accepted_genotypes,]

metadata_GR001.qc.CD38neg<-subset(metadata_GR001.qc.accepted.genotypes,cell_type=="CD34+CD38-")
metadata_GR001.qc.CD38pos<-subset(metadata_GR001.qc.accepted.genotypes,cell_type=="CD34+CD38+")
  
table(metadata_GR001.qc.CD38neg$genotype)
table(metadata_GR001.qc.CD38pos$genotype)

f_CD38neg<-8.93
f_CD38pos<-80.62
  
#Multiply CD34+CD38- cells by the CD34+CD38- cell ratio 13.59/(83.01+13.59)=0.1406832
metadata_GR001.qc.CD38neg_df<-as.data.frame(table(metadata_GR001.qc.CD38neg$genotype))
metadata_GR001.qc.CD38neg_df$Freq<-metadata_GR001.qc.CD38neg_df$Freq*(f_CD38neg/(f_CD38neg+f_CD38pos))

#Multiply CD34+CD38+ cells by the CD34+CD38+ cell ratio 83.01/(83.01+13.59)=0.8593168
metadata_GR001.qc.CD38pos_df<-as.data.frame(table(metadata_GR001.qc.CD38pos$genotype))
metadata_GR001.qc.CD38pos_df$Freq<-metadata_GR001.qc.CD38pos_df$Freq*(f_CD38pos/(f_CD38neg+f_CD38pos))

metadata_GR001.merged<-merge(metadata_GR001.qc.CD38neg_df,metadata_GR001.qc.CD38pos_df,by="Var1")

metadata_GR001.merged$Freq.total<-metadata_GR001.merged$Freq.x+metadata_GR001.merged$Freq.y
metadata_GR001.merged$fraction<-metadata_GR001.merged$Freq.total/sum(metadata_GR001.merged$Freq.total)

metadata_GR001.merged
