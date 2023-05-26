#FigureS5j - GATA1/CEBPA ratios from isogenic TP53-mutant cell lines from Ebert dataset (Science, 2019)
#cell lines from Ebert dataset (Science, 2019): GEO:GSE131592
#Author: Alba Meira

colData<-read.table("Ebert_GSE131592/colData_Ebert.txt",header = T,sep="\t") 
counts.normalized<- read.table("Ebert_GSE131592/counts.normalized.Ebert.txt",header = T,sep="\t")

#Re-order counts.normalized table
counts.normalized<-counts.normalized[,colData$Sample_ID]
log2.counts<-log2(counts.normalized+1)

colData$status<-"TP53Mutant"
colData$status[colData$Genotype=="WT"]<-"TP53_WT"

colData$CEBPA<-t(log2.counts["CEBPA",]+0.001)
colData$GATA1<-t(log2.counts["GATA1",]+0.001)
colData$ratio<-colData$CEBPA/colData$GATA1

library(ggplot2)
ggplot(colData, aes(x=status, y=ratio,fill=status)) + geom_boxplot()+
   theme_classic()+
   scale_fill_manual(values = c("grey","red"))+
   #geom_jitter(colData, mapping=aes(x=status, y=ratio), position=position_jitter(width=0.1, height=0), size=2,alpha=0.5) +
   theme(axis.text.y = element_text(size=20),legend.position = "none",
        axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_blank())

ggsave('FigureS5m_right.png',height=4,width=2)

t.test(colData[colData$status=="TP53Mutant",c("ratio")],colData[colData$status=="TP53_WT",c("ratio")])
#p-value = 0.05