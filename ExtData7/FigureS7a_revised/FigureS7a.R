#FigureS7a - summarize preleukemic donors and genotypes  

preleukemic.integration<-read.table("preleukemic_integration_metadata.qc.revised.csv",sep=",",header = T)
preleukemic.integration<-subset(preleukemic.integration,genotype.classification=="preleukemic")
dim(preleukemic.integration) #880 cells

#Pie-chart 1
table(preleukemic.integration$stage)

preleukemic.integration$evolution<-"not_classified"
preleukemic.integration$evolution[preleukemic.integration$stage %in% c("GR003_AP","IF0391")]<-"hemizygous"
preleukemic.integration$evolution[preleukemic.integration$stage=="GR006_AP"]<-"JAK2_negative"
preleukemic.integration$evolution[preleukemic.integration$stage %in% c("GR001","GR007_AP")]<-"parallel_evol"
preleukemic.integration$evolution[preleukemic.integration$stage %in% c("GR002","GR005_AP","IF0131","IF0318","SB5702")]<-"biallelic"

#Pie-chart 2
table(preleukemic.integration$evolution)

preleukemic.integration$WTvsmutant<-"mutant"
preleukemic.integration$WTvsmutant[preleukemic.integration$Genotype_curated=="WT"]<-"WT"
table(preleukemic.integration$WTvsmutant)

#Pie-chart 3
table(preleukemic.integration$evolution,preleukemic.integration$WTvsmutant)


