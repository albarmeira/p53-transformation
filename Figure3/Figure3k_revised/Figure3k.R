#Figure 3k

# Load packages
library(SingCellaR)
library(mGSZ)
library(fgsea)
library(plyr)
library(ggplot2)
library(pheatmap)

##############################################################
################ PREPARE PRE-RANKED GENE LIST ################
##############################################################

# Load R object
load(file = "../../../Sources/preleukemic.integration.revised.rdata")

# Retrieve metadata
df.pheno <- preleukemic.integration@meta.data

# Retrieve expression matrix
df <- get_normalized_umi(preleukemic.integration)

# Define non-ref groups, ref groups
table(df.pheno$genotype.classification)
non.ref <- "preleukemic"
ref <- "WT"

# Retrieve cell groups
    # group A (non-reference)
    cellsA <- df.pheno[which(df.pheno$genotype.classification==non.ref), "Cell"]
    length(cellsA) # 880
    
    # group B (reference)
    cellsB <- df.pheno[which(df.pheno$genotype.classification==ref), "Cell"]
    length(cellsB) # 730
    
# Retrieve overlaps with exp matrix
    # group A (non-reference)
    cellsA <- intersect(cellsA, colnames(df))
    length(cellsA) # 880
    
    # group B (reference)
    cellsB <- intersect(cellsB, colnames(df))
    length(cellsB) # 730

# Downsample: Both groups with equal n cells
    # Upper limit
    cells.max <- min(c(length(cellsA), length(cellsB)))

    # Downsample
    set.seed(1)
    cellsA <- sample(x=cellsA, size=cells.max, replace=FALSE)
    cellsB <- sample(x=cellsB, size=cells.max, replace=FALSE)

# Generate, write pre-ranked genes
identifyGSEAPrerankedGenes(
        method = c("wilcoxon"),
        objectA = preleukemic.integration,
        objectB = preleukemic.integration,
        cellsA = cellsA,
        cellsB = cellsB,
        write.to.file="pre-ranked.txt",
        min.log2FC = 0.1,
        min.expFraction = 0.01
        )

##############################################################
########################## GSEA ##############################
##############################################################

# Read files
    # reference gene list
    gene.set <- geneSetsList("h.all.v7.4.symbols.gmt")

    # Read DE file
    df <- read.table("pre-ranked.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

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
    write.table(fgseaRes, "gsea_hallmark.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

##############################################################
################ PLOT GSEA RESULTS ###########################
##############################################################

# Read GSEA results file
df <- read.table("gsea_hallmark.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)

. <- grep("|", df$pathway, fixed=TRUE, value=TRUE)

if(length(.) != 0) {

df <- df[-grep("|", df$pathway, fixed=TRUE), ]

}

# Censor non-sig NES
df$NES[which(df$padj > 0.25)] <- NA

# Subset selected pathways for plotting
pathways <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
            "HALLMARK_INTERFERON_GAMMA_RESPONSE",
            "HALLMARK_INFLAMMATORY_RESPONSE",
            "HALLMARK_TGF_BETA_SIGNALING",
            "HALLMARK_IL2_STAT5_SIGNALING",
            "HALLMARK_IL6_JAK_STAT3_SIGNALING",
            "HALLMARK_INTERFERON_ALPHA_RESPONSE"
            )
  
df <- df[which(df$pathway %in% pathways), ]

# Tabulate nes
nes <- df[,"NES",drop=FALSE]
row.names(nes) <- df$pathway
nes <- nes[pathways, ,drop=FALSE]

# Tidy gene set names
row.names(nes) <- gsub("HALLMARK_", "", row.names(nes))
row.names(nes) <- gsub("_", " ", row.names(nes))

# Definitions
data <- nes

# Plot and save: NES
pdf("Figure3k.pdf", width=4, height=4)

range <- unlist(nes)
range <- nes[!is.na(nes)]
myBreaks <- c(seq(max(abs(range))*-1, max(abs(range)), length.out=100))

pheatmap(data, cluster_rows=FALSE, cluster_cols=FALSE, scale="none", fontsize_row=9, border_color="white", legend=TRUE, breaks=myBreaks,
         display_numbers=TRUE,fontsize_number=20)

dev.off()

