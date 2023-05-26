###################################
library(SingCellaR)
library(SingleCellExperiment)
library(Rcpp)
library(Matrix)
library(matrixStats)
library(bigmemory)
library(ggplot2)
library(pbapply)
library(ComplexHeatmap)
library(circlize)
##################################
load(file = "../../../Sources/MFp53MPNAML_integration.revised.rdata")
##################################

passed.QC<-MPN.integration.all@meta.data[MPN.integration.all@meta.data$IsPassed==TRUE,]
table(passed.QC$donor.type) # 10459 AP 2056 MF 5002 normal
dim(passed.QC)

passed.QC$genotype.classification[is.na(passed.QC$genotype.classification)] <- "not_genotyped"
table(passed.QC$donor.type,passed.QC$genotype.classification)

table(passed.QC$genotype.classification)
dim(passed.QC)

#normal donor cells
normal<-passed.QC[passed.QC$donor.type=="normal",]
table(normal$donor)
dim(normal) #(Revised: 5002 cells from 9 donors)

#MF cells
MF<-passed.QC[passed.QC$donor.type=="MF",]
table(MF$donor)
dim(MF) #(Revised: 2056 cells from 8 donors)

#preleukemic cells from AP patients (TP53-WT)
preLSC<-passed.QC[passed.QC$donor.type=="AP" & passed.QC$genotype.classification=="preleukemic",]
table(preLSC$donor)
dim(preLSC) #(Revised: 1107 cells from 12 patients)

#TP53-mutant cells
TP53<-passed.QC[passed.QC$donor.type=="AP" & passed.QC$genotype.classification=="TP53_multihit_HOM"|
                  passed.QC$donor.type=="AP" & passed.QC$genotype.classification=="TP53_multihit_M2",]
table(TP53$donor)
dim(TP53) #(Revised: 7881 cells from 14 donors)


combined <- rbind(preLSC,normal,MF)


TP53_mutant <- TP53$Cell
combined.cells <- combined$Cell


length(TP53_mutant) # 7881
length(combined.cells) # 8165
 

#Remove genes targeted for transcript-specific amplification from analysis
dim(MPN.integration.all)
#(Revised: 32830, 21955)
MPN.integration.all@lsi.result <- "TP53"
MPN.integration.all <- MPN.integration.all[!rownames(MPN.integration.all) %in% 
                                              c("TP53","JAK2","CALR","TET2","U2AF1","EZH2","DNMT3A"),]
dim(MPN.integration.all)
#(Revised: 32823 21955)

#Extract gene counts
umi.dat<-get_normalized_umi(MPN.integration.all)

# Normalized counts using log2 transformation before doing DGE analysis
umi.dat.log2<-umi.dat
min_expr<-1 
umi.dat.log2[umi.dat.log2 < min_expr]<- 0
umi.dat.log2[umi.dat.log2 >= min_expr]<-log2(umi.dat.log2[umi.dat.log2 >= min_expr])

class(umi.dat.log2)
dim(umi.dat.log2) #32823 21955

MPN.integration.all@assays@data$normalized.umi <- umi.dat.log2

DEG.AP.combined <- identifyDifferentialGenes(
  objectA = MPN.integration.all,
  objectB = MPN.integration.all,
  cellsA = TP53_mutant,
  cellsB = combined.cells,
  min.log2FC = 0.3,
  min.expFraction = 0.2,
  write.to.file = "DEG.TP53mutant.vs.TP53WT.txt")

head(DEG.AP.combined)
write.table(DEG.AP.combined[1:100,],file="DEG.TP53mutant.vs.TP53WT_top100.txt",sep="\t",row.names = F)


plot_heatmap_for_differential_genes(objectA = MPN.integration.all,
                                    objectB = MPN.integration.all,
                                    cellsA = TP53_mutant,
                                    cellsB = combined.cells,
                                    groupA.name = "TP53_mutant",
                                    groupB.name = "TP53_WT",
                                    isClusterByCol = F,
                                    gene_list = DEG.AP.combined$Gene[1:100],
                                    col.low="blue",
                                    col.mid = "white",
                                    col.high = "brown"
                                    )




normal<-passed.QC[passed.QC$donor.type=="normal",]
table(normal$donor)
dim(normal) #(Revised: 5002 cells from 9 donors)

#MF cells
MF<-passed.QC[passed.QC$donor.type=="MF",]
table(MF$donor)
dim(MF) #(Revised: 2056 cells from 8 donors)

#preleukemic cells from AP patients (TP53-WT)
preLSC<-passed.QC[passed.QC$donor.type=="AP" & passed.QC$genotype.classification=="preleukemic",]
table(preLSC$donor)
dim(preLSC) #(Revised 1107 cells from 12 donors)

#TP53-mutant cells
TP53<-passed.QC[passed.QC$donor.type=="AP" & passed.QC$genotype.classification=="TP53_multihit_HOM"|
                  passed.QC$donor.type=="AP" & passed.QC$genotype.classification=="TP53_multihit_M2",]





full.meta.data <- rbind(TP53,combined)
dim(full.meta.data)
full.meta.data$new.genotype[full.meta.data$genotype.classification %in% c("not_genotyped","preleukemic")] <- "TP53_WT"
full.meta.data$new.genotype[full.meta.data$genotype.classification %in% c("TP53_multihit_HOM","TP53_multihit_M2")] <- "TP53_multihit"

table(full.meta.data$new.genotype)

####################################

expr <- get_normalized_umi(MPN.integration.all)
expr <- expr[DEG.AP.combined$Gene[1:100],full.meta.data$Cell]
mat_scaled = t(apply(expr, 1, scale))

###########################################

identical(as.character(full.meta.data$Cell),colnames(expr)) # TRUE
  
donor.type <- full.meta.data$donor.type
donor <- full.meta.data$donor
genotype <- full.meta.data$new.genotype

ha <- HeatmapAnnotation(donor.id = donor,
                             donor.type = donor.type,
                             genotype = genotype,
                        col = list(
                                   donor.type = c("AP" = "red","MF" = "forestgreen", "normal" = "grey60"),
                                   genotype = c(
                                                "TP53_WT" = "gray80",
                                                "TP53_multihit" = "#DC143C"),
                                   donor.id = c("IF0131" = "#f2bd80",
                                             "IF0391" = "#1d6d1f",
                                             "IF0318" = "#8c3ba0",
                                             "IF0393" = "#6533ed",
                                             "IF0392" = "#83e3f0",
                                             "IF0308" = "#fd5917",
                                             "GR001" = "#4f8c9d",
                                             "GR002" = "#eb1fcb",
                                             "GR006" = "#f5cdaf",
                                             "GR005" = "#9698dc",
                                             "GR003" = "#20f53d",
                                             "GR004" = "#f283e3",
                                             "GR007" = "#ffb2be",
                                             "GH001" = "#f3d426",
                                             "IF0704" = "#5ebf72",
                                             "IF0902" = "#a67649",
                                             "IF0903" = "#2f5bb1",
                                             "IF0905" = "#90a479",
                                             "HD15_8650" = "#f6932e",
                                             "IF0907" = "#d59e9a",
                                             "IF0908" = "#caf243",
                                             "Aph1" = "#38b5fc",
                                             "HD85" = "#c82565",
                                             "IF0137" = "#d6061a",
                                             "IF0157" = "#e36f6f",
                                             "IF0155" = "#1dfee1",
                                             "IF0123" = "#506356",
                                             "IF0101" = "black",
                                             "IF0602" = "#b35757",
                                             "IF0140" = "#e68c8c",
                                             "IF0138" = "#364261") ),
                        annotation_name_side = "left"
)

###############

tiff(file = "FigureS4b_heatmap.tiff",units = "in",width =8, height = 10, res=300,compression = "lzw")
Heatmap(mat_scaled,
        name = "expression",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "brown")),
        top_annotation = ha,
        show_column_names = F,
        show_row_names = T,
        cluster_columns = F,
        row_dend_reorder=TRUE,
        row_names_gp = gpar(fontsize = 6,fontface = "italic"))
dev.off()
