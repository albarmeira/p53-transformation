#FigureS9o - GSEA in Tp53-homozygous versus Tp53-heterozygous cells from pre-sAML-TP53 patients (Related to FigureS9i-m)
#Authors: Sean Wen and Alba Meira

# Load packages
library(SingCellaR)
library(mGSZ)
library(fgsea)
library(ggplot2)

#################################################################
################## CREATE PRE-RANKED GENE LIST ##################
#################################################################

# Read files
    # SingCellaR object
    load("../../../Sources/MFp53MPNAML_integration.revised.rdata")

    # Sample metadata
    df.pheno <- read.table("../../../Sources/metadata_MPNAMLp53_with_index_genotype.revised.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Subset overlapping cells with gex
df.pheno.2 <- MPN.integration.all@meta.data

cell_ids <- intersect(df.pheno$cell_id, df.pheno.2$Cell)
length(cell_ids) # 19221

index <- which(df.pheno$cell_id %in% cell_ids)
df.pheno <- df.pheno[index, ]

# Define cell groups
    # GH001
        # Group A (non-reference)
        index.1 <- which(df.pheno$stage %in% c("GH001_001", "GH001_002"))
        index.2 <- which(df.pheno$genotype.classification=="TP53_multihit_HOM")
        index <- intersect(index.1, index.2)
        length(index) # 593 cells
        cellsA.1 <- df.pheno[index, "cell_id"]

        # Group B (reference)
        index.1 <- which(df.pheno$stage %in% c("GH001_001", "GH001_002"))
        index.2 <- which(df.pheno$genotype.classification=="TP53HET")
        index <- intersect(index.1, index.2)
        length(index) # 249 cells
        cellsB.1 <- df.pheno[index, "cell_id"]
    
    # GR005
        # Group A (non-reference)
        index.1 <- which(df.pheno$stage %in% c("GR005_CP"))
        index.2 <- which(df.pheno$genotype.classification=="TP53_multihit_M2")
        index <- intersect(index.1, index.2)
        length(index) # 29 cells
        cellsA.2 <- df.pheno[index, "cell_id"]

        # Group B (reference)
        index.1 <- which(df.pheno$stage %in% c("GR005_CP"))
        index.2 <- which(df.pheno$genotype.classification=="TP53HET")
        index <- intersect(index.1, index.2)
        length(index) # 35 cells
        cellsB.2 <- df.pheno[index, "cell_id"]
        
    # Merge
        # Group A (non-reference)
        cellsA <- c(cellsA.1, cellsA.2)
        length(cellsA) # 622
        
        # Group B (non-reference)
        cellsB <- c(cellsB.1, cellsB.2)
        length(cellsB) # 284

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
        objectA = MPN.integration.all,
        objectB = MPN.integration.all,
        cellsA = cellsA,
        cellsB = cellsB,
        write.to.file=paste(path.local,"pre-ranked.txt",sep=""),
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
################## PLOT GSEA RESULTS #########################
##############################################################

# Read GSEA results
p53HETvsHOM_GSEA<-read.table(file="gsea_hallmark.txt",sep="\t",header = T)

pathways <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
              "HALLMARK_INTERFERON_GAMMA_RESPONSE",
              "HALLMARK_IL6_JAK_STAT3_SIGNALING"
              )

p53HETvsHOM_GSEA <- p53HETvsHOM_GSEA[which(p53HETvsHOM_GSEA$pathway %in% pathways), ]

p53HETvsHOM_GSEA$pathway<-gsub("HALLMARK_","",p53HETvsHOM_GSEA$pathway)
rownames(p53HETvsHOM_GSEA)<-p53HETvsHOM_GSEA$pathway
p53HETvsHOM_GSEA$Group<-"HOM"
p53HETvsHOM_GSEA$Size<-as.numeric(abs(p53HETvsHOM_GSEA$NES))
p53HETvsHOM_GSEA$Order<-length(p53HETvsHOM_GSEA$NES):1

#Re-order for plotting
p53HETvsHOM_GSEA$pathway <- factor(p53HETvsHOM_GSEA$pathway, levels = p53HETvsHOM_GSEA$pathway[order(p53HETvsHOM_GSEA$Order)])

#####################################
ggplot(p53HETvsHOM_GSEA,aes(x=-(log10(padj+0.001)),y=pathway,colour=Group))+
  theme_classic()+
  geom_segment(aes(x=0.55, xend=-(log10(padj+0.001)), y=pathway, yend=pathway), colour="black",size=0.5)+
  geom_point(aes(size=Size))+
  scale_color_manual(values=c("red3"))+
  scale_x_continuous(limits = c(0.5,2.0), expand = c(0, 0))+
  theme(axis.text.x = element_text(size=12))

ggsave(filename="FigureS9o.pdf", width=5, height=2)
