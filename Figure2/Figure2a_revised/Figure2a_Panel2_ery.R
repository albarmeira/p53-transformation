################################################
#Figure2a; erythroid score panel
################################################

# Load packages
library(SingCellaR)
library(plyr)
library(diffusionMap)
library("gg3D")
library(SingleCellExperiment)

# Load modified function to compute gene score
source("../../../Sources/SingleCellPlots_Edit.R")

# Load SingCellaR object
load(file = "../../../Sources/AP_HTMPNAML_Harmony.revised.rdata")

# Retrieve diffusion map coordinates
DiffMap<-AP.HTMPNAML@diffusionmap.result

# Subset selected genotypes to plot
DiffMap<-subset(DiffMap,genotype.collapsed=="TP53_multihit"|genotype.collapsed=="preleukemic")

# Compute gene scores
scores <- plot_diffusionmap_label_by_multiple_gene_sets(
                AP.HTMPNAML,
                gmt.file="../../../Sources/genesets/MC.human.signature.genes.gmt",
                show_gene_sets=c("Ery_AML"),
                return.values=TRUE
                )

#Ery_AML<-c("KLF1","GATA1","ITAG2B","FCER1A","EPOR","HBD")

scores$Cell <- row.names(scores)

DiffMap <- join(DiffMap, scores, by="Cell", type="left")

DiffMap <- DiffMap[order(DiffMap$Ery_AML), ]

# Plot
qplot(x=0, y=0, z=0, geom="blank") +
  axes_3D() +
  theme_void()

ggplot(DiffMap, aes(x=-DC1, y=-DC2, z=-DC4, color=Ery_AML)) +
  geom_point(DiffMap, mapping=aes(x=-DC1, y=-DC2, z=-DC4, color=Ery_AML))+
  theme_void() +
  axes_3D() +
  scale_colour_gradientn(colours = c("gray85","red","red"),values=c(0,0.2,1)) +
  # geom_point(size=0.03) +
  stat_3D()

ggsave(filename="Figure2a_Panel2_ery.png",width = 8,height = 8)
