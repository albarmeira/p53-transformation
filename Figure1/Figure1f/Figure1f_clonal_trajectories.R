#Read IF0131 table

IF0131_CNV_genotype<-read.table(file="cnv_summary_IF0131_subclusters_genotype.txt",sep="\t",header = T)
IF0131_CNV_genotype$CNV_genotype<-paste(IF0131_CNV_genotype$genotype,IF0131_CNV_genotype$cluster_corrected,sep="_")
table(IF0131_CNV_genotype$CNV_genotype)

#Summary of preleukemic clusters
#preleukemic_chr3_chr5_chr7 - 1
#preleukemic_normal - 205
#preleukemic_chr3_chr5 - 42

#Summary of TP53-het clusters
#TP53HET_chr3_chr5 - 8
#TP53HET_chr3_chr5_chr7 - 5
#TP53HET_normal - 1

#Summary of TP53_multi_hit clusters
#TP53_multihit_M2_normal - 14
#TP53_multihit_M2_chr3_chr5 - 486
#TP53_multihit_M2_chr3_chr5_chr7 - 226

#Read GR001 table

GR001_CNV_genotype<-read.table(file="cnv_summary_GR001_subclusters_genotype.txt",sep="\t",header = T)
GR001_CNV_genotype$CNV_genotype<-paste(GR001_CNV_genotype$genotype,GR001_CNV_genotype$cluster_corrected,sep="_")
table(GR001_CNV_genotype$CNV_genotype)

#Summary of preleukemic clusters
#preleukemic_chr5_chr17_chr21 - 5
#preleukemic_chr5_chr7_chr9_chr21 - 1
#preleukemic_normal - 105

#Summary of TP53-het clusters
#TP53HET_chr5_chr17_chr21 - 22
#TP53HET_normal - 1
#TP53HET_normal - 1

#Summary of TP53_multi_hit_HOM clusters
#TP53_multihit_HOM_chr5_chr17_chr21 - 581
#TP53_multihit_HOM_chr5_chr7_chr9_chr21 - 5
#TP53_multihit_HOM_normal - 4

#Summary of TP53_multihit_M2_chr5_chr17_chr21 clusters
#TP53_multihit_M2_chr5_chr17_chr21 - 28
#TP53_multihit_M2_chr5_chr7_chr9_chr21 - 103
#TP53_multihit_M2_normal - 27