################################################
#Figure2a
################################################
library(SingCellaR)
load(file = "../../../Sources/AP_HTMPNAML_Harmony.revised.rdata")
library(diffusionMap)
library("gg3D")
################################################
object<-AP.HTMPNAML
######################################
DiffMap<-object@diffusionmap.result

DiffMap<-subset(DiffMap,genotype.collapsed=="TP53_multihit"|genotype.collapsed=="preleukemic") #cells with TP53-het and unassigned genotypes are not show in this plot for clarity of representation
dim(DiffMap) #8988
table(DiffMap$donor) #14 donors

# Bring preLSC to the front
levels <- c("TP53_multihit", "preleukemic")
DiffMap$genotype.collapsed <- factor(DiffMap$genotype.collapsed, levels=levels)
DiffMap <- DiffMap[order(DiffMap$genotype.collapsed), ]

# An empty plot with 3 axes
qplot(x=0, y=0, z=0, geom="blank") +
  theme_void() +
  axes_3D()

ggplot(DiffMap, aes(x=-DC1, y=-DC2, z=-DC4, color=genotype.collapsed)) +
  theme_void() +
  axes_3D() +
  stat_3D() +
    scale_colour_manual(values=c("red3", "deepskyblue"))+
  theme(legend.position="none")

ggsave(filename="Figure2a_Panel1.png",width = 10,height = 10)
