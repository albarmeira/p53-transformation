#FigureS5f

################################################
library(diffusionMap)
library(SingCellaR)
library("gg3D")
library(ggplot2)
################################################
load(file = "../../../Sources/AP_HTMPNAML_Harmony.revised.rdata")
################################################

AP.HTMPNAML@meta.data
table(AP.HTMPNAML@meta.data$stage)

GR003_AP.MEP<-subset(AP.HTMPNAML@meta.data,stage=="GR003_AP" & CD38>1500 & CD45RA<3000 & CD123<800)
IF0131.MEP<-subset(AP.HTMPNAML@meta.data,stage=="IF0131" & CD38>10000 & CD45RA<7000 & CD123<1000)
SB5702.MEP<-subset(AP.HTMPNAML@meta.data,stage=="SB5702" & CD38>10000 & CD45RA<7000 & CD123<1000)
GR001.MEP<-subset(AP.HTMPNAML@meta.data,stage=="GR001" & CD38>5000 & CD45RA<5000 & CD123<1000)
GR005_AP.MEP<-subset(AP.HTMPNAML@meta.data,stage=="GR005_AP" & CD38>370 & CD45RA<3500 & CD123<1000)
GR007_AP.MEP<-subset(AP.HTMPNAML@meta.data,stage=="GR007_AP" & CD38>2000 & CD45RA<4000 & CD123<1000)
GR006_AP.MEP<-subset(AP.HTMPNAML@meta.data,stage=="GR006_AP" & CD38>1500 & CD45RA<5000 & CD123<1200)
IF0391.MEP<-subset(AP.HTMPNAML@meta.data,stage=="IF0391" & CD38>10000 & CD45RA<5000 & CD123<1000)
GR002.MEP<-subset(AP.HTMPNAML@meta.data,stage=="GR002" & CD38>1400 & CD45RA<3000 & CD123<230)
IF0318.MEP<-subset(AP.HTMPNAML@meta.data,stage=="IF0318" & CD38>6000 & CD45RA<5000 & CD123<1000)
IF0392.MEP<-subset(AP.HTMPNAML@meta.data,stage=="IF0392" & CD38>5000 & CD45RA<4000 & CD123<800)
GST010.MEP<-subset(AP.HTMPNAML@meta.data,stage=="GST010" & CD38>5000 & CD45RA<4000 & CD123<800)
JB4211.MEP<-subset(AP.HTMPNAML@meta.data,stage=="JB4211" & CD38>5000 & CD45RA<4000 & CD123<800)
IF0393.MEP<-subset(AP.HTMPNAML@meta.data,stage=="IF0393" & CD38>6000 & CD45RA<5000 & CD123<1000)


#AHGH003 doesn't have any MEPs by FACS

MEPs<-c(as.character(GR003_AP.MEP$Cell),
        as.character(IF0131.MEP$Cell),
        as.character(SB5702.MEP$Cell),
        as.character(GR001.MEP$Cell),
        as.character(GR005_AP.MEP$Cell),
        as.character(GR007_AP.MEP$Cell),
        as.character(GR006_AP.MEP$Cell),
        as.character(IF0391.MEP$Cell),
        as.character(GR002.MEP$Cell),
        as.character(IF0318.MEP$Cell))

AP.HTMPNAML@diffusionmap.result$MEP<-"FALSE"
AP.HTMPNAML@diffusionmap.result$MEP[AP.HTMPNAML@diffusionmap.result$Cell %in% MEPs]<-"MEP"
table(AP.HTMPNAML@diffusionmap.result$MEP) #904 cells classified as MEPs included

######################################
DiffMap<-AP.HTMPNAML@diffusionmap.result

# An empty plot with 3 axes
qplot(x=0, y=0, z=0, geom="blank") +
  theme_void() +
  axes_3D()

# Bring MEPs to the front
#levels <- c("notMEP","MEP")
#DiffMap$MEP <- factor(DiffMap$MEP, levels=levels)
#DiffMap <- DiffMap[order(DiffMap$MEP), ]


ggplot(DiffMap, aes(x=-DC1, y=-DC2, z=-DC4,color=MEP, alpha=MEP)) +
  theme_void() +
  axes_3D() +
  stat_3D() +
  scale_colour_manual(values=c("grey80","red"))+
  geom_point(size=3)+
  theme(legend.position="none")

ggsave(filename="FigureS5f_MEPs.png",width = 10,height = 10)


