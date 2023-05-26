#Deciphering TP53-mutant cancer evolution with single-cell multi-omics
#Author: Alba Rodriguez-Meira
##############################

#Summarize karyotipic abnormalities in TP53-mutant cells from TP53-sAML patients

GH001<-read.table(file="cnv_tables/cnv_summary_GH001_subclusters_genotype.txt",sep="\t",header = T)
GH001<-subset(GH001,stage=="GH001_003" & genotype=="TP53_multihit_HOM")
table(GH001$cluster_corrected) #560 - 100%

GR003<-read.table(file="cnv_tables/cnv_summary_GR003_subclusters_genotype.revised.txt",sep="\t",header = T)
GR003<-subset(GR003,stage=="GR003_AP" & genotype=="TP53_multihit_HOM")
table(GR003$cluster_corrected)
#clone1 normal 
#801     26 


IF0392<-read.table(file="cnv_tables/cnv_summary_IF0392_subclusters_genotype.txt",sep="\t",header = T)
table(IF0392$genotype)
IF0392<-subset(IF0392,genotype=="TP53_multihit_HOM")
table(IF0392$cluster_corrected)
#cluster1 cluster2  #cluster2 seeems to be driven by lower transcriptome quality, but the chromosomal abnormalities are similar
#cluster1 and cluster2 will be merged in a single cluster for this analysis
#665      246

IF0308<-read.table(file="cnv_tables/cnv_summary_IF0308_subclusters_genotype.txt",sep="\t",header = T)
table(IF0308$genotype)
IF0308<-subset(IF0308,genotype=="TP53_multihit_HOM")
table(IF0308$cluster_corrected)
#cluster1 
#684

GR004<-read.table(file="cnv_tables/cnv_summary_GR004_subclusters_genotype.txt",sep="\t",header = T)
table(GR004$genotype)
GR004<-subset(GR004,genotype=="TP53_multihit_M2")
table(GR004$cluster_corrected)
#cluster1 
#677 

IF0131<-read.table(file="cnv_tables/cnv_summary_IF0131_subclusters_genotype.txt",sep="\t",header = T)
table(IF0131$genotype)
IF0131<-subset(IF0131,genotype=="TP53_multihit_M2")
table(IF0131$cluster_corrected)
#chr3_chr5 chr3_chr5_chr7         normal 
#486            226             14 

IF0393<-read.table(file="cnv_tables/cnv_summary_IF0393_subclusters_genotype.txt",sep="\t",header = T)
table(IF0393$genotype)
IF0393<-subset(IF0393,genotype=="TP53_multihit_HOM")
table(IF0393$cluster_corrected)
#del20 
#650 

GR001<-read.table(file="cnv_tables/cnv_summary_GR001_subclusters_genotype.txt",sep="\t",header = T)
table(GR001$genotype)
GR001<-subset(GR001,genotype=="TP53_multihit_HOM"|genotype=="TP53_multihit_M2")
table(GR001$cluster_corrected)
#chr5_chr17_chr21 chr5_chr7_chr9_chr21               normal 
#609                  108                   31


GR005<-read.table(file="cnv_tables/cnv_summary_GR005_subclusters_genotype.txt",sep="\t",header = T)
table(GR005$genotype)
GR005<-subset(GR005,stage=="GR005_AP" & genotype=="TP53_multihit_M2")
table(GR005$cluster_corrected)
#chr1_chr3_chr6_chr9           chr1_chr6              normal 
#246                  55                  50 

GR007<-read.table(file="cnv_tables/cnv_summary_GR007_subclusters_genotype.txt",sep="\t",header = T)
table(GR007$genotype)
GR007<-subset(GR007,stage=="GR007_AP" & genotype=="TP53_multihit_M2")
table(GR007$cluster_corrected)
#cluster1   normal 
#424       54

GR006<-read.table(file="cnv_tables/cnv_summary_GR006_subclusters_genotype.txt",sep="\t",header = T)
table(GR006$genotype)
GR006<-subset(GR006,stage=="GR006_AP" & genotype=="TP53_multihit_M2")
table(GR006$cluster_corrected)

#del5_del7_amp8_del12               normal 
#367                   11 

IF0391<-read.table(file="cnv_tables/cnv_summary_IF0391_subclusters_genotype.txt",sep="\t",header = T)
table(IF0391$genotype)
IF0391<-subset(IF0391,genotype=="TP53_multihit_HOM")
table(IF0391$cluster_corrected)

#del3_abn9_amp21          normal 
#107               6 

GR002<-read.table(file="cnv_tables/cnv_summary_GR002_subclusters_genotype.txt",sep="\t",header = T)
table(GR002$genotype)
GR002<-subset(GR002,genotype=="TP53_multihit_M2")
table(GR002$cluster_corrected)
#abn19_abn21 abn19_abn21_amp2           normal 
#226               34               11 

IF0318<-read.table(file="cnv_tables/cnv_summary_IF0318_subclusters_genotype.txt",sep="\t",header = T)
table(IF0318$genotype)
IF0318<-subset(IF0318,genotype=="TP53_multihit_M2")
table(IF0318$cluster_corrected)

#del5_del9 del5_del9_del13_amp19                normal 
#11                   488                     5

summary<-read.table(file="Figure1g_Summary.txt",sep="\t",header = T)
summary<-summary[,c("Sample","Normal_pc","CNV_pc_SB1","CNV_pc_SB2")]
library(reshape2)
plot<-melt(summary,c("Sample"))


# Stacked
library(ggplot2)
ggplot(plot, aes(fill=variable, y=value, x=Sample)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=c("cyan3","deeppink4","violetred2"))+
  theme_classic()+
  scale_y_continuous(limits = c(0,1.02), expand = c(0, 0))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=10,angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=15),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background = element_blank(),
        legend.position = "none")

ggsave(file="Figure1g.png",width = 3,height = 3)

