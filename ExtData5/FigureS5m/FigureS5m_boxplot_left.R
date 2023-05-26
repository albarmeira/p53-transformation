#FigureS5j - erythroid-signature heatmap expression from isogenic TP53-mutant 
#cell lines from Ebert dataset (Science, 2019): GEO:GSE131592
#Author: Alba Meira

colData<-read.table("Ebert_GSE131592/colData_Ebert.txt",header = T,sep="\t") 
counts.normalized<- read.table("Ebert_GSE131592/counts.normalized.Ebert.txt",header = T,sep="\t")

#Re-order counts.normalized table
counts.normalized<-counts.normalized[,colData$Sample_ID]

#Subset erythoid.signature genes
erythroid_genes<-c("KLF1", "GATA1", "ZFPM1", "GATA2", "GYPA", "TFRC")
counts.normalized.erythroid<-counts.normalized[erythroid_genes,]

#Log-normalize and calculate z-scores
counts.log.erythroid<-log1p(counts.normalized.erythroid)

#Calculate erythroid score
df.small <- apply(counts.log.erythroid, 2, mean)
df.small <- data.frame("sample.id"=names(df.small), "score"=df.small, stringsAsFactors=FALSE)
row.names(df.small) <- NULL

names(df.small)[which(names(df.small)=="score")] <- "erythroid.score"
colData <- merge(colData, df.small, by.x="Sample_ID",by.y="sample.id", all=T)

colData$status<-"TP53Mutant"
colData$status[colData$Genotype=="WT"]<-"TP53_WT"

library(ggplot2)
ggplot(colData, aes(x=status, y=erythroid.score,fill=status)) + geom_boxplot(outlier.shape = NA)+
  theme_classic()+
  geom_jitter(colData, mapping=aes(x=status, y=erythroid.score), position=position_jitter(width=0.1, height=0), size=2,alpha=0.5) +
  scale_fill_manual(values = c("grey","red"))+
  theme(axis.text.y = element_text(size=20),legend.position = "none",
        axis.text.x=element_blank(),axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_blank())

ggsave('FigureS5m_left.png',height=4,width=2)

t.test(colData[colData$status=="TP53Mutant",c("erythroid.score")],colData[colData$status=="TP53_WT",c("erythroid.score")])
#p-value = 0.01963
