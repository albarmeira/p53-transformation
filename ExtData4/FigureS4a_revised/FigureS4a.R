#Integration of MF MolCell 2019 samples and p53MPNAML dataset using SingleCellaR
#Data integration
#Author:Alba Rodriguez-Meira
#Date: 24th March 2021
#Modified last: 12th October 2021
####################################################

# FDG python setup
#library(reticulate)
#conda_create("r-reticulate", python_version="3.6")
#py_install("fa2", envname="r-reticulate")
#py_install("networkx==2.1", envname="r-reticulate")
#py_module_available("fa2")
#py_module_available("networkx")
#py_module_available("numpy")

# Load  packages
library(SingCellaR)
library(ggplot2)
######################################################
# Load  R object
load(file="../../../Sources/MFp53MPNAML_integration.revised.rdata")
######################################################

runFA2_ForceDirectedGraph(MPN.integration.all,dim_reduction_method="pca",n.dims.use = 30,useIntegrativeEmbeddings = F)

p<-plot_forceDirectedGraph_label_by_a_feature_of_interest(MPN.integration.all,feature = "donor")
p
ggsave(file="FigureS4a.donor.png",width=6,height = 5)

p<-plot_forceDirectedGraph_label_by_a_feature_of_interest(MPN.integration.all,feature = "genotype.collapsed")
p+scale_color_manual(values = c("greenyellow","blue","red","orange","grey40"),na.value="white")
ggsave(file="FigureS4a.genotype.png",width=6,height = 5)

#################