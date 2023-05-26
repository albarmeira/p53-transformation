################################################

#Creating an InferCNV object based on your three required inputs:
#the read count matrix, cell type annotations, and the gene ordering file:

# create the infercnv object

library("infercnv")
# mainDir<-getwd()
# subDir<-"inferCNV_IF0131_subclones"
# dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
# 
# infercnv_obj.IF0131 = CreateInfercnvObject(raw_counts_matrix="IF0131_tables/counts.matrix.txt",
#                                            annotations_file="IF0131_tables/annotation_infercnv.txt",
#                                            delim="\t",
#                                            gene_order_file="genes_subset.txt",
#                                            ref_group_names=c("normal"))
# 
# # perform infercnv operations to reveal cnv signal
# infercnv_obj.IF0131 = infercnv::run(infercnv_obj.IF0131,
#                                     cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
#                                     out_dir="inferCNV_IF0131_subclones",  # dir is auto-created for storing outputs
#                                     cluster_by_groups=T,   # cluster
#                                     denoise=T,
#                                     analysis_mode='subclusters',
#                                     #HMM_report_by="cell",
#                                     HMM=T,
# 				                            num_threads = 8)
# 
# save(infercnv_obj.IF0131,file="infercnv_obj_IF0131_subclusters.rdata")

load(file="RawData/objects/infercnv_obj_IF0131_subclusters.rdata")


infercnv_obj.IF0131@options$HMM_report_by == "cell"
by_cells = TRUE
cnv_summary<-infercnv_obj.IF0131@tumor_subclusters$subclusters$patient
cnv_summary<-data.frame(unlist(cnv_summary))
cnv_summary$subcluster<-sapply(strsplit(as.character(rownames(cnv_summary)),".MPNAML",fixed = TRUE),"[[",1)
#cnv_summary$subcluster<-gsub("patient.","IF0131.",cnv_summary$subcluster)
cnv_summary$cell_id<-sapply(strsplit(as.character(rownames(cnv_summary)),".MPN",fixed = TRUE),"[[",2)
cnv_summary$cell_id<-gsub("AML","MPNAML",cnv_summary$cell_id)

cnv_summary<-cnv_summary[,c("subcluster","cell_id")]

#write.table(cnv_summary,file="RawData/CNVSummary_subclusters/cnv_summary_IF0131_subclusters.txt",sep="\t",row.names=F)

##################################################################################
############## Heatmap plot ################ 
# to reproduce the heatmap output by InferCNV 

cnv_summary$cluster_corrected<-NA
cnv_summary$cluster_corrected[cnv_summary$subcluster=="patient.1.1.1.1"]<-"chr3_chr5"
cnv_summary$cluster_corrected[cnv_summary$subcluster=="patient.1.1.1.2"]<-"chr3_chr5_chr7"

cnv_summary$cluster_corrected[cnv_summary$subcluster=="patient.1.1.2.1"]<-"chr3_chr5"
cnv_summary$cluster_corrected[cnv_summary$subcluster=="patient.1.1.2.2"]<-"chr3_chr5"

cnv_summary$cluster_corrected[cnv_summary$subcluster=="patient.1.2.1.1"]<-"normal"
cnv_summary$cluster_corrected[cnv_summary$subcluster=="patient.1.2.1.2"]<-"normal"

cnv_summary$cluster_corrected[cnv_summary$subcluster=="patient.1.2.2.1"]<-"normal"
cnv_summary$cluster_corrected[cnv_summary$subcluster=="patient.1.2.2.2"]<-"normal"

# Generate two identically sorted lists #
#infercnv_observations_step21 <- as.data.frame(read.table(file = "RawData/subclones_folders/inferCNV_IF0131_subclones/infercnv.observations.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE))
infercnv_observations_step21 <- as.data.frame(read.table(file = "RawData/subclones_folder/inferCNV_IF0131_subclones/infercnv.observations.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE))
#infercnv_observations_step19 <- as.data.frame(read.table(file = "RawData/subclones_folders/inferCNV_IF0131_subclones/infercnv.19_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE))
infercnv_observations_step19 <- as.data.frame(read.table(file = "RawData/subclones_folder/inferCNV_IF0131_subclones/infercnv.19_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE))

infercnv_observations_step21_sorted <- data.frame(row.names = row.names(infercnv_observations_step19))
for (cell in colnames(infercnv_observations_step19))
  infercnv_observations_step21_sorted[cell] <- infercnv_observations_step21[, c(cell)]

infercnv_observations<-infercnv_observations_step21_sorted[,order(ncol(infercnv_observations_step21_sorted):1)]

rm(infercnv_observations_step19)
rm(infercnv_observations_step21)
rm(infercnv_observations_step21_sorted)
# 
gene_order_file <- read.table(file = "RawData/genes_subset.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE) 
colnames(gene_order_file) <- c("Gene", "Chr", "Start", "End") 
gene_order_file_used <- gene_order_file[gene_order_file$Gene %in% rownames(infercnv_observations),] 

infercnv_observations_order <- infercnv_observations[match(rownames(infercnv_observations), gene_order_file_used$Gene), ] 
infercnv_observations_order_t <- t(as.matrix(infercnv_observations_order)) 

# make an annotation data frame 

#annotations_file <- read.table(file = "annotations_file8.txt", sep = "\t", header = F, stringsAsFactors = FALSE) 
#annotations_file <- annotations_file[!annotations_file$V2 %in% c("Normal_BM_1","Normal_BM_2", "Normal_BM_3"),] 

annotations_file<-cnv_summary

annotations_file<-annotations_file[match(rownames(infercnv_observations_order_t),annotations_file$cell_id),]
metadata<-read.table(file="../../../Sources/metadata_MPNAMLp53_with_index_genotype.revised.txt",sep="\t",header = T)
metadata<-subset(metadata,cell_id %in% annotations_file$cell_id)
metadata<-metadata[,c("cell_id","genotype.classification")]

merged<-metadata

#Re-order metadata matrix according to annotation file
merged<-merged[match(annotations_file$cell_id, merged$cell_id),]
annotations_file$genotype.corrected<-merged$genotype.classification


df_table_rows <- data.frame(annotations_file$subcluster) 
rownames(df_table_rows) <- annotations_file$cell_id 

colnames(df_table_rows)<- "subcluster" 

df_table_rows$genotype<- annotations_file$genotype.corrected
df_table_rows$cluster_corrected<- annotations_file$cluster_corrected

table(df_table_rows$cluster_corrected)
table(df_table_rows$genotype)

df_table_rows$cluster <- factor(df_table_rows$subcluster, levels = c("patient.1.1.1.1", "patient.1.1.1.2","patient.1.1.2.1", "patient.1.1.2.2","patient.1.2.1.1","patient.1.2.1.2","patient.1.2.2.1","patient.1.2.2.2")) 
df_table_rows$genotype <- factor(df_table_rows$genotype, levels = c("preleukemic", "TP53_multihit_M2","TP53HET","unassigned","NA")) 
df_table_rows$cluster_corrected<- factor(df_table_rows$cluster_corrected, levels = c("chr3_chr5", "chr3_chr5_chr7","normal")) 

df_table_rows$subcluster<-NULL

df_table_cols <- data.frame(gene_order_file_used$Chr) 
rownames(df_table_cols) <- gene_order_file_used$Gene 
colnames(df_table_cols)<- "Chr" 

library(RColorBrewer) 

col_vector <- c("seagreen2", "lightpink4", "gainsboro", "coral1", "thistle1", "brown3", "floralwhite", "midnightblue", "lavenderblush", "salmon2", "firebrick2", 
                "darkorange1", "dodgerblue1", "linen", "aquamarine4", "lightsalmon", "rosybrown", "magenta1", "sienna", "lightcyan", "darkorange2", "darkorange4") 


ann_colors <- list(cluster = c("patient.1.1.1.1" = "red", "patient.1.1.1.2" = "yellow", "patient.1.1.2.1" = "blue", "patient.1.1.2.2" = "green", "patient.1.2.1.1" = "deepskyblue" , "patient.1.2.1.2" = "brown", "patient.1.2.2.1" = "orange", "patient.1.2.2.2"= "grey60"),
                   genotype=c("NA"="white","preleukemic"="blue","TP53HET"="orange","TP53_multihit_M2"="red","unassigned"="grey90"),
                   cluster_corrected=c("chr3_chr5"="red","normal"="grey60","chr3_chr5_chr7"="brown"),
                   Chr = c("chr1" = col_vector[1], "chr2" = col_vector[2], "chr3" = col_vector[3], "chr4" = col_vector[4], "chr5" = col_vector[5], "chr6" = col_vector[6], 
                           "chr7" = col_vector[7], "chr8" = col_vector[8], "chr9" = col_vector[9], 
                           "chr10" = col_vector[10], "chr11" = col_vector[11], "chr12" = col_vector[12], 
                           "chr13" = col_vector[13], "chr14" = col_vector[14], "chr15" = col_vector[15], 
                           "chr16" = col_vector[16], "chr17" = col_vector[17], "chr18" = col_vector[18], 
                           "chr19" = col_vector[19], "chr20" = col_vector[20], "chr21" = col_vector[21], 
                           "chr22" = col_vector[22])) 

library(pheatmap)
png(file ="inferCNV_heatmap_IF0131_Sean.png",width=800,height = 800)
pheatmap(infercnv_observations_order_t, show_rownames=F, show_colnames=F, cluster_cols=F,  cluster_rows=F,
         color=colorRampPalette(rev(brewer.pal(n = 12, name = "RdBu")))(100),cexRow=0.75,cexCol=0.75, border_color =NA,  
         annotation_col = df_table_cols, annotation_row = df_table_rows, annotation_colors=ann_colors,  
         main = "InferCNV IF0131", 
         fontsize_row=12, key=TRUE,symkey=FALSE,keysize = 0.8) 

dev.off()



