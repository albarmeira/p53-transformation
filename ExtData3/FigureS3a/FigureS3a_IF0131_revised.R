################################################################
library(ggplot2)
####

metadata<-read.table('../../../Sources/metadata_MPNAMLp53_with_index_genotype.revised.txt',header = T,sep = '\t')
rownames(metadata)<-metadata$cell_id
colnames(metadata)
table(metadata$donor)

#Subset IF0131 donor
metadata<-subset(metadata,donor=="IF0131")

metadata$Genotype_labels <- factor(metadata$Genotype_labels)
metadata$Genotype_curated <- factor(metadata$Genotype_curated)

table(metadata$Genotype_labels)
table(metadata$Genotype_curated)
table(metadata$genotype.classification)

#Remove rows in which genotype.classification=NA (non-classified genotypes)

metadata<-metadata[!is.na(metadata$genotype.classification), ]
metadata<-subset(metadata,genotype.classification=="TP53_multihit_M2")

################################################################
### Plot index sorting data
################################################################

#Check that there's no negative values to plot in log-axis
min(metadata$CD38,na.rm=TRUE)
min(metadata$CD90,na.rm=TRUE)
min(metadata$CD123,na.rm=TRUE)
min(metadata$CD45RA,na.rm=TRUE)

#Compute psedovalues for negative ones
#CD38_pseudovalue=-(min(metadata$CD38,na.rm=TRUE))+100
CD90_pseudovalue=-(min(metadata$CD90,na.rm=TRUE))+100
CD123_pseudovalue=-(min(metadata$CD123,na.rm=TRUE))+100
#CD45RA_pseudovalue=-(min(metadata.patient.qc$CD45RA,na.rm=TRUE))+100

ggplot(metadata,aes(y=CD38,x=CD34,fill=cell_type))+
  geom_point(size=5,shape=21,stroke=1) + theme_bw()+
  scale_x_continuous(limits = c(1,120000))+
  scale_y_log10(limits = c(1000,500000))+
  geom_hline(yintercept = 10000)

ggplot(metadata,aes(y=CD38,x=CD34,fill=genotype.classification))+
  geom_point(size=3,shape=21,stroke=1,alpha=0.9) + theme_bw() +
  scale_x_continuous(limits = c(1,120000))+
  scale_y_log10(limits = c(1000,500000))+
  geom_hline(yintercept = 10000)+
  #scale_fill_manual(values = c("deepskyblue","brown",'red','orange'))
  scale_fill_manual(values = c("red"))

CD38_threshold<-10000

CD38pos<-subset(metadata,CD38>CD38_threshold)
CD38neg<-subset(metadata,CD38<=CD38_threshold)

ggplot(CD38neg,aes(y=CD90+CD90_pseudovalue,x=CD45RA,fill=genotype.classification))+
  theme_bw()+
  theme(legend.position = "none")+
  stat_density_2d(color="red",size=1) +
  geom_point(size=3,shape=21,stroke=1,alpha=0.2) +
  scale_y_log10(limits = c(500,500000))+
  scale_x_log10(limits = c(500,100000))+
  geom_hline(yintercept = 50000+CD90_pseudovalue)+
  geom_vline(xintercept = 5000)+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x=element_blank())+
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y=element_blank()) +
  theme(panel.border = element_rect(size=4,colour = "black"))+
  theme(text = element_blank())+
  theme(axis.ticks = element_blank())+
  scale_fill_manual(values = c("red"))

ggsave('FigureS3a_IF0131_CD38neg.png',height=5,width=5)

ggplot(CD38pos,aes(y=CD123+CD123_pseudovalue,x=CD45RA,fill=genotype.classification))+
  theme_bw() +
  theme(legend.position = "none")+
  stat_density_2d(color="red",size=1) +
  geom_point(size=3,shape=21,stroke=1,alpha=0.2) +
  scale_y_log10(limits=c(1000,50000))+
  scale_x_log10(limits=c(1000,100000))+
  geom_hline(yintercept = 1200+CD123_pseudovalue)+
  geom_vline(xintercept = 5000)+
  theme(panel.grid.major.x = element_blank(),panel.grid.minor.x=element_blank())+
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y=element_blank()) +
  theme(panel.border = element_rect(size=4,colour = "black"))+
  theme(text = element_blank())+
  theme(axis.ticks = element_blank())+
  scale_fill_manual(values = c("red"))

ggsave('FigureS3a_IF0131_CD38pos.png',height=5,width=5)

#############################################################
