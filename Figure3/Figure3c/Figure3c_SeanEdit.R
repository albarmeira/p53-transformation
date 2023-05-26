####################################################
#Figure 3C
#Author:Alba Rodriguez-Meira
#Date: 12th March 2021
#Modified last: 4th October 2021
####################################################

library(ggplot2)

#Import list of preleukemic cells
metadata<-read.table("../../../Sources/metadata_MPNAMLp53_with_index_genotype.revised.txt",header = T,sep="\t")

#Check that there's no negative values to plot in log-axis
min(metadata$CD38,na.rm=TRUE) #-2256.31
min(metadata$CD90,na.rm=TRUE) #-3823.967
min(metadata$CD123,na.rm=TRUE) #-2566.551
min(metadata$CD45RA,na.rm=TRUE) #-14749.16

#Compute psedovalues for negative ones
CD38_pseudovalue=-(min(metadata$CD38,na.rm=TRUE))+100
CD90_pseudovalue=-(min(metadata$CD90,na.rm=TRUE))+100
CD123_pseudovalue=-(min(metadata$CD123,na.rm=TRUE))+100
CD45RA_pseudovalue=-(min(metadata$CD45RA,na.rm=TRUE))+100

############################################################
################### IF0131 (CD38-) #########################
############################################################

#IF0131 genotype plots
IF0131<-subset(metadata,donor=="IF0131")
CD38_threshold<-10000
IF0131.CD38neg<-subset(IF0131,CD38<=CD38_threshold)
IF0131.CD38pos<-subset(IF0131,CD38>CD38_threshold)
# table(IF0131.CD38neg$genotype.classification)
# include<-c("preleukemic","TP53_multihit_M2")
# IF0131.CD38neg<-subset(IF0131.CD38neg,genotype.classification %in% include)

p<-ggplot(IF0131.CD38neg,aes(y=CD90+CD90_pseudovalue,x=CD45RA+CD45RA_pseudovalue,fill=genotype.classification))+
  geom_point(size=5,shape=21,stroke=1) +
  theme_classic()+
  theme(legend.position = "none",
        axis.line = element_line(colour = 'black', size = 2))+
  #scale_x_continuous(limits = c(10000,100000))+
  #scale_y_log10(limits = c(1000,75000))+
  scale_x_log10()+
  scale_y_log10()+
  geom_hline(yintercept = 10000+CD90_pseudovalue,size=1.5)+
  geom_vline(xintercept = 7000+CD45RA_pseudovalue,size=1.5)+
  scale_fill_manual(values=c("blue","red","orange","grey80"))

p
ggplot_build(p)$layout$panel_scales_x[[1]]$range$range #xmin=10^4.178946 #xmax=10^4.966436
ggplot_build(p)$layout$panel_scales_y[[1]]$range$range #ymin=10^2.911946 #ymax=10^5.392709

#ggsave(file="points/IF0131.CD38neg.png",width = 5,height = 5)

############################################################
#################### GR001 (CD38-) #########################
############################################################

#GR001 genotype plots
GR001<-subset(metadata,stage=="GR001")
CD38_threshold<-5000
GR001.CD38neg<-subset(GR001,CD38<=CD38_threshold)
GR001.CD38pos<-subset(GR001,CD38>CD38_threshold)

CD90_threshold<-10000
CD45RA_threshold<-5000
CD123_threshold<-1000

p<-ggplot(GR001.CD38neg,aes(y=CD90+CD90_pseudovalue,x=CD45RA+CD45RA_pseudovalue,fill=genotype.classification))+
  geom_point(size=5,shape=21,stroke=1) +
  theme_classic()+
  theme(legend.position = "none",axis.line = element_line(colour = 'black', size = 2))+
  scale_x_log10()+
  scale_y_log10()+
  geom_hline(yintercept = CD90_threshold+CD90_pseudovalue,size=1.5)+
  geom_vline(xintercept = CD45RA_threshold+CD45RA_pseudovalue,size=1.5)+
  scale_fill_manual(values=c("blue","red3","red","orange","grey80"))

p
ggplot_build(p)$layout$panel_scales_x[[1]]$range$range #xmin=10^4.184758 #xmax=10^4.577334
ggplot_build(p)$layout$panel_scales_y[[1]]$range$range #ymin=10^3.412812 #ymax=10^5.431015

############################################################
################### IF0318 (CD38-) #########################
############################################################

#IF0318 genotype plots
IF0318<-subset(metadata,stage=="IF0318")
CD38_threshold<-6000
IF0318.CD38neg<-subset(IF0318,CD38<=CD38_threshold)
IF0318.CD38pos<-subset(IF0318,CD38>CD38_threshold)

CD90_threshold<-5000
CD45RA_threshold<-5000
CD123_threshold<-1000

p<-ggplot(IF0318.CD38neg,aes(y=CD90+CD90_pseudovalue,x=CD45RA+CD45RA_pseudovalue,fill=genotype.classification))+
  geom_point(size=5,shape=21,stroke=1) +
  theme_classic()+
  theme(legend.position = "none",axis.line = element_line(colour = 'black', size = 2))+
  scale_x_continuous(limits = c(15000,25000))+
  scale_y_log10()+
  geom_hline(yintercept = CD90_threshold+CD90_pseudovalue,size=1.5)+
  geom_vline(xintercept = CD45RA_threshold+CD45RA_pseudovalue,size=1.5)+
  scale_fill_manual(values=c("blue","red","orange","grey80"))

p
ggplot_build(p)$layout$panel_scales_x[[1]]$range$range #xmin=15282.19 #xmax=19849.16
ggplot_build(p)$layout$panel_scales_y[[1]]$range$range #ymin=10^3.479989 #ymax=10^5.112116

#ggsave(file="points/IF0318.CD38neg.png",width = 5,height = 5)

############################################################
################### IF0391 (CD38-) #########################
############################################################

#IF0391 genotype plots
IF0391<-subset(metadata,stage=="IF0391")
CD38_threshold<-10000
IF0391.CD38neg<-subset(IF0391,CD38<=CD38_threshold)
IF0391.CD38pos<-subset(IF0391,CD38>CD38_threshold)

CD90_threshold<-10000
CD45RA_threshold<-5000
CD123_threshold<-1000

p<-ggplot(IF0391.CD38neg,aes(y=CD90+CD90_pseudovalue,x=CD45RA+CD45RA_pseudovalue,fill=genotype.classification))+
  geom_point(size=5,shape=21,stroke=1) +
  theme_classic()+
  theme(legend.position = "none",axis.line = element_line(colour = 'black', size = 2))+
  scale_x_log10()+
  scale_y_log10()+
  geom_hline(yintercept = CD90_threshold+CD90_pseudovalue,size=1.5)+
  geom_vline(xintercept = CD45RA_threshold+CD45RA_pseudovalue,size=1.5)+
  scale_fill_manual(values=c("blue","red3","grey80"))

p
ggplot_build(p)$layout$panel_scales_x[[1]]$range$range #xmin=10^4.179889 #xmax=10^4.608397
ggplot_build(p)$layout$panel_scales_y[[1]]$range$range #ymin=10^3.569079 #ymax=10^5.040351

#ggsave(file="IF0391.CD38neg.png",width = 5,height = 5)
