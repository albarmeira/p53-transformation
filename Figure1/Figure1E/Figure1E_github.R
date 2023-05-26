load(file="metadata_GR006.qc.rdata")

table(metadata_GR006.qc$cell_type)

#Include only CD34+ sorted cells in the count
subset<-c("CD34+")
metadata_GR006.qc<-metadata_GR006.qc[metadata_GR006.qc$cell_type %in% subset,]
table(metadata_GR006.qc$genotype)

accepted_genotypes<-c("JAK2_V617F_WT_TP53_p141_WT_TP53_p350_WT",
                      "JAK2_V617F_HET_TP53_p141_WT_TP53_p350_WT",
                      "JAK2_V617F_HOM_TP53_p141_WT_TP53_p350_WT",
                      "JAK2_V617F_WT_TP53_p141_WT_TP53_p350_HET",
                      "JAK2_V617F_WT_TP53_p141_HET_TP53_p350_HET")

metadata_GR006.qc.accepted.genotypes<-metadata_GR006.qc[metadata_GR006.qc$genotype %in% accepted_genotypes,]
  
#get frequency from each genotype
metadata_GR006.merged<-as.data.frame(table(metadata_GR006.qc.accepted.genotypes$genotype))

metadata_GR006.merged$fraction<-metadata_GR006.merged$Freq/sum(metadata_GR006.merged$Freq)

metadata_GR006.merged
