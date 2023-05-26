#Annotated healthy reference map for LSI integration using Greenleaf's dataset
#Figure2c
#Author:Alba Rodriguez-Meira
#Date: 19th November 2019
#Modified last: 14th September 2021
####################################################

####################################################
#Import object
#Object generated in 01_Greenleaf_Healthy_Reference.R

load(file="../../../Sources/Greenleaf_reference_annotated.Rdata")
library(ggplot2)
####################################################

#Re-colour BioClassification published by Granja,Klemm,McGinnis et al,
#Nature Biotechnology 2019, DOI: 10.1038/s41587-019-0332-7 
df$colorv2<-df$color
df$colorv2<-as.character(df$colorv2)
df$colorv2[df$color=="01_HSC"]<-"HSC"
df$colorv2[df$color=="02_Early.Eryth"]<-"early_ery"
df$colorv2[df$color=="03_Late.Eryth"]<-"late_ery"
df$colorv2[df$color=="04_Early.Baso"]<-"unk"
df$colorv2[df$color=="05_CMP.LMPP"]<-"CMP_LMPP"
df$colorv2[df$color=="06_CLP.1"]<-"lymph"
df$colorv2[df$color=="07_GMP"]<-"myelo"
df$colorv2[df$color=="08_GMP.Neut"]<-"myelo"
df$colorv2[df$color=="09_pDC"]<-"myelo"
df$colorv2[df$color=="10_cDC"]<-"myelo"
df$colorv2[df$color=="11_CD14.Mono.1"]<-"myelo"
df$colorv2[df$color=="12_CD14.Mono.2"]<-"myelo"
df$colorv2[df$color=="13_CD16.Mono"]<-"myelo"
df$colorv2[df$color=="14_Unk"]<-"myelo"
df$colorv2[df$color=="15_CLP.2"]<-"lymph"
df$colorv2[df$color=="16_Pre.B"]<-"lymph"
df$colorv2[df$color=="17_B"]<-"lymph"
df$colorv2[df$color=="18_Plasma"]<-"lymph"
df$colorv2[df$color=="19_CD8.N"]<-"lymphT"
df$colorv2[df$color=="20_CD4.N1"]<-"lymphT"
df$colorv2[df$color=="21_CD4.N2"]<-"lymphT"
df$colorv2[df$color=="22_CD4.M"]<-"lymphT"
df$colorv2[df$color=="23_CD8.EM"]<-"lymphT"
df$colorv2[df$color=="24_CD8.CM"]<-"lymphT"
df$colorv2[df$color=="25_NK"]<-"lymphT"
df$colorv2[df$color=="26_Unk"]<-"unk"

colours<-c("greenyellow","orange","deepskyblue","red","purple","blue","forestgreen","grey90")

ggplot(df,aes(x,y,color=colorv2)) + 
  geom_point() + 
  theme_bw() + 
  xlab("UMAP Dimension 1") + 
  ylab("UMAP Dimension 2")+
  theme_classic()+
  theme(legend.position = "none", axis.title =element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank())+
  scale_color_manual(values = colours) #all

ggsave('Figure2c.png',height=10,width=10)

