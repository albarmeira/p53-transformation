################################################
#Figure2g - CEBPA expression
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

# Retrieve diffusion map coordintaes
DiffMap<-AP.HTMPNAML@diffusionmap.result

# Subset selected genotypes to plot
DiffMap<-subset(DiffMap,genotype.collapsed=="TP53_multihit"|genotype.collapsed=="preleukemic")

# Compute gene scores
scores <- plot_diffusionmap_label_by_genes(
                AP.HTMPNAML,
                gene_list=c("CEBPA"),
                return.values=TRUE
                )

DiffMap <- join(DiffMap, scores[,c("Cell", "CEBPA")], by="Cell", type="left")

DiffMap <- DiffMap[order(DiffMap$CEBPA), ]

# Plot
qplot(x=0, y=0, z=0, geom="blank") +
  axes_3D() +
  theme_void()

ggplot(DiffMap, aes(x=-DC1, y=-DC2, z=-DC4, color=CEBPA)) +
  geom_point(DiffMap, mapping=aes(x=-DC1, y=-DC2, z=-DC4, color=CEBPA))+
  theme_void() +
  axes_3D() +
  scale_colour_gradientn(colours = c("gray87","blue","red"),values=c(0,0.01,1)) +
  # geom_point(size=0.03) +
  stat_3D()

ggsave(filename="Figure2g_CEBPA.png", width=8, height=8)

#####
p.x[[j]] <- ggplot(my.dat, aes_string(x="-DC1", y="DC2", z="DC4",colour = genes.x,alpha=genes.x)) + 
  theme_void()+axes_3D()+stat_3D()+ 
  scale_colour_gradientn(colours = c("gray87",point.color1,point.color2),values=c(0,0.01,1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())+
  guides(alpha = FALSE)
