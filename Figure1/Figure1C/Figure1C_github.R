load(file="metadata_IF0391.qc.rdata")

table(metadata_IF0391.qc$cell_type)

#Include only CD34+ sorted cells in the count
subset<-c("CD34+")
metadata_IF0391.qc<-metadata_IF0391.qc[metadata_IF0391.qc$cell_type %in% subset,]
table(metadata_IF0391.qc$genotype)

accepted_genotypes<-c("JAK2_WT_TP53_p175_WT_TET2_p1863_WT",
                      "JAK2_HET_TP53_p175_WT_TET2_p1863_WT",
                      "JAK2_HET_TP53_p175_WT_TET2_p1863_HET",
                      "JAK2_HET_TP53_p175_HET_TET2_p1863_HET",
                      "JAK2_HET_TP53_p175_HET_TET2_p1863_HOM",
                      "JAK2_HET_TP53_p175_HOM_TET2_p1863_HOM")

metadata_IF0391.qc.accepted.genotypes<-metadata_IF0391.qc[metadata_IF0391.qc$genotype %in% accepted_genotypes,]
  
#get frequency from each genotype
metadata_IF0391.merged<-as.data.frame(table(metadata_IF0391.qc.accepted.genotypes$genotype))

metadata_IF0391.merged$fraction<-metadata_IF0391.merged$Freq/sum(metadata_IF0391.merged$Freq)

metadata_IF0391.merged

