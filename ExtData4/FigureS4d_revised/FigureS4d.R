# Load packages
library(SingCellaR)
library(mGSZ)
library(fgsea)
library(plyr)
library(pheatmap)

############################################################
############################ DE ############################
############################################################

# Load R object
load("../../../Sources/MFp53MPNAML_integration.revised.rdata")

# Retrieve metadata
df.pheno <- MPN.integration.all@meta.data

# Retrieve expression matrix
df <- get_normalized_umi(MPN.integration.all)

# Merge genotype w/ cluster assignment for TP53_multihit
table(df.pheno$genotype.collapsed, df.pheno$merged.louvain)
df.pheno$genotype.groups <- NA

index <- df.pheno$genotype.collapsed=="TP53_multihit" & df.pheno$merged.louvain=="cl1_LSC"
df.pheno$genotype.groups[index] <- "TP53_multihit_LSC"

index <- df.pheno$genotype.collapsed=="TP53_multihit" & df.pheno$merged.louvain=="cl2_erythroid"
df.pheno$genotype.groups[index] <- "TP53_multihit_ery"

index <- df.pheno$genotype.collapsed=="MF"
df.pheno$genotype.groups[index] <- "MF"

index <- df.pheno$genotype.collapsed=="WT"
df.pheno$genotype.groups[index] <- "normal"

index <- df.pheno$genotype.collapsed=="preLSC"
df.pheno$genotype.groups[index] <- "preLSC"

table(df.pheno$genotype.groups) # 4546 TP53_multihit_LSC
                                # 3130 TP53_multihit_ery
                                # 2056 MF
                                # 5002 normal
                                # 1107 preLSC

# Define cell groups
non.ref <- c("TP53_multihit_LSC", "TP53_multihit_LSC", "TP53_multihit_LSC",
             "TP53_multihit_ery", "TP53_multihit_ery", "TP53_multihit_ery"
             )
ref <- c("MF", "normal", "preLSC",
         "MF", "normal", "preLSC"
         )

# DE
for(i in 1:length(non.ref)) {

    # Retrieve cell groups
        # group A (non-reference)
        cellsA <- df.pheno[which(df.pheno$genotype.groups==non.ref[i]), "Cell"]
        length(cellsA)
        
        # group B (reference)
        cellsB <- df.pheno[which(df.pheno$genotype.groups==ref[i]), "Cell"]
        length(cellsB)
        
    # Retrieve overlaps with exp matrix
        # group A (non-reference)
        cellsA <- intersect(cellsA, colnames(df))
        length(cellsA)
        
        # group B (reference)
        cellsB <- intersect(cellsB, colnames(df))
        length(cellsB)

    # Downsample: Both groups with equal n cells
        # Upper limit
        cells.max <- min(c(length(cellsA), length(cellsB)))

        # Downsample
        set.seed(1)
        cellsA <- sample(x=cellsA, size=cells.max, replace=FALSE)
        cellsB <- sample(x=cellsB, size=cells.max, replace=FALSE)

    # Generate, write pre-ranked genes
    file <- paste("pre-ranked_", non.ref[i], " vs ", ref[i], ".txt", sep="")

    identifyGSEAPrerankedGenes(
            method = c("wilcoxon"),
            objectA = MPN.integration.all,
            objectB = MPN.integration.all,
            cellsA = cellsA,
            cellsB = cellsB,
            write.to.file=paste("pre-ranked_", non.ref[i], " vs ", ref[i], ".txt", sep=""),
            min.log2FC = 0.1,
            min.expFraction = 0.01
            )

}

############################################################
##################### GSEA: HALLMARK #######################
############################################################

# Read reference gene list
gene.set <- geneSetsList("h.all.v7.4.symbols.gmt")

# Define DE files
files <- list.files()
files <- files[grep("pre-ranked", files, fixed=TRUE)]

for(i in 1:length(files)) {

    # Read R object
    file <- files[i]
    df <- read.table(file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

    # Recode Inf,-Inf values
    df$log2FC[which(df$FoldChange=="Inf")] <- log2(df$ExpA[which(df$FoldChange=="Inf")] + 1)
    df$log2FC[which(df$FoldChange==0)] <- log2(df$ExpB[which(df$FoldChange==0)] + 1) * -1
    
    # Recode 0 p-values
    pval <- min(df$adjusted.pval[which(df$adjusted.pval != 0)])
    df$adjusted.pval[which(df$adjusted.pval == 0)] <- pval
    
    # Transform p-values
    df$adjusted.pval.log <- -log10(df$adjusted.pval)
    df$adjusted.pval.log <- ifelse(df$log2FC < 0, df$adjusted.pval.log* -1, df$adjusted.pval.log * 1)
     
    # Rank list (from lowest to highest! Opposite to GSEA GUI)
    df <- df[order(df$adjusted.pval.log), ]

    # Retrieve ranked genes
    genes.ranked <- df$adjusted.pval.log
    names(genes.ranked) <- df$Gene

    # Capitalize gene names
    names(genes.ranked) <- toupper(names(genes.ranked))

    # GSEA
        # Lastest version do not have eps option
        set.seed(123)
        fgseaRes <- fgsea(pathways=gene.set, stats=genes.ranked, minSize=0, maxSize=5000, nperm=1000)

        # Convert to data frame
        fgseaRes <- as.data.frame(fgseaRes)
        
        # Adjust for mutliple testing (default: fdr)
        #fgseaRes$padj.fdr  <-  p.adjust(fgseaRes$pval, method="fdr", n=length(fgseaRes$padj))

        # Order by p-values
        fgseaRes <- fgseaRes[order(fgseaRes$pval), ]

        # Reformat leading edge column
        . <- as.character(fgseaRes$leadingEdge)
        . <- gsub("c(", "", ., fixed=TRUE)
        . <- gsub(")", "", ., fixed=TRUE)
        . <- gsub("\"", "", ., fixed=TRUE)
        . <- gsub(" ", "" ,.)
        . <- gsub(",", "|", ., fixed=TRUE)
        fgseaRes$leadingEdge <- .

        # Save file
        file <- gsub("pre-ranked", "gsea", file, fixed=TRUE)
        write.table(fgseaRes, file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
        

}

############################################################
######################### GSEA #############################
############################################################

# Read reference gene list
gene.set <- geneSetsList("h.all.v7.4.symbols_GENTLES_LEUKEMIC_STEM_CELL_UP.gmt")

# Define DE files
files <- list.files()
files <- files[grep("pre-ranked", files, fixed=TRUE)]

for(i in 1:length(files)) {

    # Read R object
    file <- files[i]
    df <- read.table(file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

    # Recode Inf,-Inf values
    df$log2FC[which(df$FoldChange=="Inf")] <- log2(df$ExpA[which(df$FoldChange=="Inf")] + 1)
    df$log2FC[which(df$FoldChange==0)] <- log2(df$ExpB[which(df$FoldChange==0)] + 1) * -1
    
    # Recode 0 p-values
    pval <- min(df$adjusted.pval[which(df$adjusted.pval != 0)])
    df$adjusted.pval[which(df$adjusted.pval == 0)] <- pval
    
    # Transform p-values
    df$adjusted.pval.log <- -log10(df$adjusted.pval)
    df$adjusted.pval.log <- ifelse(df$log2FC < 0, df$adjusted.pval.log* -1, df$adjusted.pval.log * 1)
     
    # Rank list (from lowest to highest! Opposite to GSEA GUI)
    df <- df[order(df$adjusted.pval.log), ]

    # Retrieve ranked genes
    genes.ranked <- df$adjusted.pval.log
    names(genes.ranked) <- df$Gene

    # Capitalize gene names
    names(genes.ranked) <- toupper(names(genes.ranked))

    # GSEA
        # Lastest version do not have eps option
        set.seed(123)
        fgseaRes <- fgsea(pathways=gene.set, stats=genes.ranked, minSize=0, maxSize=5000, nperm=1000)

        # Convert to data frame
        fgseaRes <- as.data.frame(fgseaRes)
        
        # Adjust for mutliple testing (default: fdr)
        #fgseaRes$padj.fdr  <-  p.adjust(fgseaRes$pval, method="fdr", n=length(fgseaRes$padj))

        # Order by p-values
        fgseaRes <- fgseaRes[order(fgseaRes$pval), ]

        # Reformat leading edge column
        . <- as.character(fgseaRes$leadingEdge)
        . <- gsub("c(", "", ., fixed=TRUE)
        . <- gsub(")", "", ., fixed=TRUE)
        . <- gsub("\"", "", ., fixed=TRUE)
        . <- gsub(" ", "" ,.)
        . <- gsub(",", "|", ., fixed=TRUE)
        fgseaRes$leadingEdge <- .

        # Save file
        file <- gsub("pre-ranked", "gsea", file, fixed=TRUE)
        write.table(fgseaRes, file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
        

}

############################################################
####################### GSEA HEATMAP #######################
############################################################

# Define GSEA result files
files <- c("gsea_TP53_multihit_LSC vs normal.txt",
           "gsea_TP53_multihit_LSC vs MF.txt",
           "gsea_TP53_multihit_LSC vs preLSC.txt",
           "gsea_TP53_multihit_ery vs normal.txt",
           "gsea_TP53_multihit_ery vs MF.txt",
           "gsea_TP53_multihit_ery vs preLSC.txt"
           )
labels <- c("LSC vs HD",
            "LSC vs MF",
            "LSC vs preLSC",
            "Ery vs HD ",
            "Ery vs MF ",
            "Ery vs preLSC "
            )

df.list <- list()

for(i in 1:length(files)) {
  
  file <- files[i]
  df <- read.table(file, sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
  
  . <- grep("|", df$pathway, fixed=TRUE, value=TRUE)
  
  if(length(.) != 0) {
    
    df <- df[-grep("|", df$pathway, fixed=TRUE), ]
    
  }
  
  
  df.list[[i]] <- df
  
}

length(df.list)

# Tabulate padj
. <- lapply(df.list, function(x) {(x[,c("pathway", "padj")])})
names(.[[1]])[2] <- labels[1]
names(.[[2]])[2] <- labels[2]
names(.[[3]])[2] <- labels[3]
names(.[[4]])[2] <- labels[4]
names(.[[5]])[2] <- labels[5]
names(.[[6]])[2] <- labels[6]
. <- Reduce(function(x,y) join(x=x, y=y, by="pathway", type="full"), .)
row.names(.) <- .$pathway
.$pathway <- NULL
padj <- .

# Tabulate nes
. <- lapply(df.list, function(x) {(x[,c("pathway", "NES")])})
names(.[[1]])[2] <- labels[1]
names(.[[2]])[2] <- labels[2]
names(.[[3]])[2] <- labels[3]
names(.[[4]])[2] <- labels[4]
names(.[[5]])[2] <- labels[5]
names(.[[6]])[2] <- labels[6]
. <- Reduce(function(x,y) join(x=x, y=y, by="pathway", type="full"), .)
row.names(.) <- .$pathway
.$pathway <- NULL
nes <- .

# Censor non-sig NES
threshold <- 0.25
nes[padj > threshold] <- NA

# Remove non-sig. gene sets
. <- apply(nes, 1, function(x) {sum(!is.na(x)) != 0})
nes <- nes[.,]

# Reorder for aesthetic purpose
#nes$up <- apply(nes, 1, function(x) {sum(x > 0, na.rm=TRUE)})
#nes$down <- apply(nes, 1, function(x) {sum(x < 0, na.rm=TRUE)})
#nes$net <- nes$up - nes$down
#nes <- nes[order(nes$net, decreasing=TRUE), ]
#nes <- nes[, -which(names(nes) %in% c("up", "down", "net"))]
nes$sum <- rowSums(nes, na.rm=TRUE)
nes <- nes[order(nes$sum, decreasing=TRUE), ]
nes$sum <- NULL

# Tidy gene set names
row.names(nes) <- gsub("HALLMARK_", "", row.names(nes))
row.names(nes) <- gsub("_", " ", row.names(nes))

# Subset selected pathways
pathways <- c("G2M CHECKPOINT",
              "MTORC1 SIGNALING",
              "E2F TARGETS",
              "PI3K AKT MTOR SIGNALING",
              "IL6 JAK STAT3 SIGNALING",
              "KRAS SIGNALING UP",
              "IL2 STAT5 SIGNALING",
              "APOPTOSIS",
              "INTERFERON ALPHA RESPONSE",
              "INTERFERON GAMMA RESPONSE",
              "GENTLES LEUKEMIC STEM CELL UP"
              )
nes <- nes[pathways, ]

# Definitions
data <- nes

# Plot and save: NES
pdf("FigureS4d.pdf", width=4, height=3)

pheatmap(data, cluster_rows=FALSE, cluster_cols=FALSE, scale="none", fontsize_row=9, border_color="white", legend=TRUE, angle_col="90")

dev.off()


