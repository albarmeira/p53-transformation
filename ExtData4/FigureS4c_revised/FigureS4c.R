# Load packages
library(mGSZ)
library(SingCellaR)
library(GeneOverlap)

# Read files
    # GMT
    gene.set <- geneSetsList("Master.gmt")
    
    # Convert mouse to human names, if any
    gene.set <- lapply(gene.set, toupper)

    # DE output (from FigS4b)
    df <- read.table("../FigureS4b_revised/DEG.TP53mutant.vs.TP53WT.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
    # SingCellaR object
    load("../../../Sources/MFp53MPNAML_integration.revised.rdata")

############################################################
################### HYPERGEOMETRIC TEST ####################
############################################################

# Retrieve n expressed genes
df.feature <- MPN.integration.all@genes.info
index <- which(df.feature$IsExpress==TRUE)
n.expr <- length(index)

# Find overlaps: Up
    # Retrieve reference genes
    gene_short_names.ref <- unlist(gene.set["Fischer_et_al_p53_up"])
        
    # Retrieve HSPC genes to compare against
    index <- which(df$log2FC < 0)
    length(index)
    df.small <- df[index, ]
    
    # Tabulate shared/unique
    go.obj <- newGeneOverlap(gene_short_names.ref,
                             df.small$Gene,
                             genome.size=n.expr
                             )
    
    # Hyper-geometric test
    go.obj <- testGeneOverlap(go.obj)
    
    # Save as new object
    go.obj_up <- go.obj
    
# Find overlaps: Down
    # Retrieve reference genes
    gene_short_names.ref <- unlist(gene.set["Fischer_et_al_p53_down"])
        
    # Retrieve HSPC genes to compare against
    index <- which(df$log2FC > 0)
    length(index)
    df.small <- df[index, ]
    
    # Tabulate shared/unique
    go.obj <- newGeneOverlap(gene_short_names.ref,
                             df.small$Gene,
                             genome.size=n.expr
                             )
    
    # Hyper-geometric test
    go.obj <- testGeneOverlap(go.obj)

    # Save as new object
    go.obj_down <- go.obj
    
# Find overlaps: Up + down
    # Retrieve reference genes
    gene_short_names.ref_up <- unlist(gene.set["Fischer_et_al_p53_up"])
    gene_short_names.ref_down <- unlist(gene.set["Fischer_et_al_p53_down"])
    gene_short_names.ref <- unique(c(gene_short_names.ref_up, gene_short_names.ref_down))
    
    # Retrieve HSPC genes to compare against
    df.small <- df
    
    # Tabulate shared/unique
    go.obj <- newGeneOverlap(gene_short_names.ref,
                             df.small$Gene,
                             genome.size=n.expr
                             )
    go.obj
    
    # Hyper-geometric test
    go.obj_merged <- testGeneOverlap(go.obj)
    go.obj_merged
    
######################################################

# Save merged results
results <- data.frame("unique_ref"=length(setdiff(go.obj_merged@listA, go.obj_merged@intersection)),
                      "unique_hspc"=length(setdiff(go.obj_merged@listB, go.obj_merged@intersection)),
                      "overlap"=length(go.obj_merged@intersection),
                      "jaccard"=go.obj_merged@Jaccard,
                      "pval"=go.obj_merged@pval
                      )

write.table(results, "HyperGeometric Test.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Save overlap, consistent direction
    # Up in WT, down in MT
    gene_short_names <- go.obj_up@intersection
    
    results.1 <- data.frame("gene_short_name"=gene_short_names,
                            "direction"="Overlap: Fischer up_HSPC MUT down",
                            stringsAsFactors=FALSE
                            )

    # Down in WT, up in MUT
    gene_short_names <- go.obj_down@intersection
    
    results.2 <- data.frame("gene_short_name"=gene_short_names,
                            "direction"="Overlap: Fischer down_HSPC MUT up",
                            stringsAsFactors=FALSE
                            )
    
    # Up in WT, up in MT
    gene_short_names <- intersect(go.obj_up@listA,
                                c(go.obj_down@intersection, go.obj_down@listB)
                                )
    results.3 <- data.frame("gene_short_name"=gene_short_names,
                            "direction"="Overlap: Fischer up_HSPC MUT up",
                            stringsAsFactors=FALSE
                            )
   
    # Down in WT, down in MT
    gene_short_names <- intersect(go.obj_down@listA,
                                c(go.obj_up@intersection, go.obj_up@listB)
                                )
    results.4 <- data.frame("gene_short_name"=gene_short_names,
                            "direction"="Overlap: Fischer down_HSPC MUT down",
                            stringsAsFactors=FALSE
                            )
                            
    # No overlap: unique to Fischer
    gene_short_names <- setdiff(go.obj_merged@listA, go.obj_merged@intersection)
    results.5 <- data.frame("gene_short_name"=gene_short_names,
                            "direction"="No overlap: Fischer unique",
                            stringsAsFactors=FALSE
                            )
                            
   
    # No overlap: unique to HSPCs
    gene_short_names <- setdiff(go.obj_merged@listB, go.obj_merged@intersection)
    results.6 <- data.frame("gene_short_name"=gene_short_names,
                            "direction"="No overlap: HSPC MUT unique",
                            stringsAsFactors=FALSE
                            )
                            
    # Save file
    results <- rbind.data.frame(results.1, results.2, results.3, results.4, results.5, results.6)

    write.table(results, "Fischer-HSPC MUT Direction.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

############################################################
################## GENE ONTOLOGY ANALYSIS ##################
############################################################

# Load packages
library(MARVEL)

# Read direction file
df <- read.table("Fischer-HSPC MUT Direction.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
# Define genes
gene_short_names <- df[grep("^Overlap", df$direction), "gene_short_name"]

# GO
marvel <- list()
marvel <- BioPathways(MarvelObject=marvel,
                      custom.genes=gene_short_names,
                      species="human"
                      )

# Save file
write.table(marvel$DE$BioPathways$Table, "GO_All Overlaps.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

############################################################
############## GENE ONTOLOGY ANALYSIS: PLOT ###############
############################################################

# Load packages
library(ggplot2)

# Read all overlaps file
df <- read.table("GO_All Overlaps.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
    
# Define terms to plot
marvel <- list()
marvel$DE$BioPathways$Table <- df
go.terms <- c("signal transduction by p53 class mediator",
              "cell cycle arrest",
              "DNA damage response, signal transduction by p53 class mediator resulting in cell cycle arrest",
              "signal transduction involved in mitotic G1 DNA damage checkpoint",
              "regulation of intrinsic apoptotic signaling pathway",
              "cellular response to external stimulus",
              "response to radiation",
              "positive regulation of catabolic process",
              "response to oxidative stress",
              "apoptotic mitochondrial changes",
              "cellular response to UV",
              "cellular response to hypoxia",
              "positive regulation of intracellular transport",
              "positive regulation of protein transport",
              "cellular response to chemical stress"
              )

# Barplot
    # Subset terms to plot
    df <- df[which(df$Description %in% go.terms), ]
    
    # Set factor levels
    levels <- rev(df$Description)
    df$Description <- factor(df$Description, levels=levels)
    df <- df[order(df$Description), ]

    # Definition
    data <- df
    x <- data$Description
    y <- -log10(data$p.adjust)
    maintitle <- ""
    ytitle <- "-log10(FDR)"
    xtitle <- ""
    
    # Plot
    plot <- ggplot() +
        geom_bar(data, mapping=aes(x=x, y=y),stat="identity", fill="red3") +
        labs(title=maintitle, x=xtitle, y=ytitle) +
        theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border=element_blank(),
            plot.title=element_text(hjust = 0.5, size=12),
            plot.subtitle=element_text(hjust = 0.5, size=12),
            axis.line.y.left = element_line(color="black"),
            axis.line.x = element_line(color="black"),
            axis.title=element_text(size=12),
            axis.text.x=element_text(size=10, colour="black"),
            axis.text.y=element_blank(),
            legend.position="none"
            ) +
        coord_flip()

    # Save plot
    ggsave("GO_All Overlaps_Top Pathways.pdf", plot, width=2.5, height=3)
