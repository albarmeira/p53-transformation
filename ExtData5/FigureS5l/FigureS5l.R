#FigureS5l - CEBPA and GATA1 levels in BeatAML dataset stratified by CEBPA mutational status
#Authors: Alba Meira

# Load packages
library(data.table)
library(plyr)
library(ggplot2)

# Read file
df <- read.table("../../../Sources/BeatAML/Master.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Set factor level
table(df$genotype.TP53)
df$genotype.TP53 <- factor(df$genotype.TP53,
                           levels=c("Wildtype", "Mutant (1 hit)", "Mutant (2 hits)"),
                           labels=c("WT", "Mut", "Mut")
)
###########################################################

# Read expression files
    # Matrix
    exp <- as.data.frame(fread("../../../01.Preprocessing/BeatAML/FPKM/FPKM.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))

    row.names(exp) <- exp$gene_id
    exp$gene_id <- NULL

    # Gene metadata
    exp.feature <- read.table("../../../01.Preprocessing/BeatAML/FPKM/FPKM_featureData.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Transform and censor values
exp <- log2(exp + 1)
#exp[exp < 1] <- 0

# Subset genes
genes<-c("CEBPA","GATA1")
gene_ids <- exp.feature[which(exp.feature$gene_short_name %in% genes), "gene_id"]
gene_short_names <- exp.feature[which(exp.feature$gene_short_name %in% genes), "gene_short_name"]
exp <- exp[gene_ids, ]
row.names(exp) <- gene_short_names
exp <- as.data.frame(t(exp))
exp$patient.id <- row.names(exp)

# Annotate
df <- join(df, exp, by="patient.id", type="left")
###########################################################
## GATA1

# Boxplot
    # Definition
    data <- df
    x <- data$genotype.TP53
    y <- data[["GATA1"]]
    z <- data$genotype.TP53
    maintitle <- ""
    ytitle <- "log2(FPKM + 1)"
    xtitle <- ""
    #legendtitle <- "CEBPA Genotype"
    #xlabels <- n.cells$label
    #fivenum(y) ; ymin <- -0.5 ; ymax <- 1.0 ; yinterval <- 0.25

    # Plot
    plot <- ggplot() +
      geom_boxplot(data, mapping=aes(x=x, y=y, fill=z), outlier.shape=NA, show.legend = FALSE) +
      scale_fill_manual(values=c("grey60","red3"))+
      #ylim(c(0,10))+
      #scale_y_log10(expand=c(0,0))+
      geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.1,alpha=0.5) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border=element_blank(),
            axis.line.y.left = element_line(color="black"),
            axis.line.x = element_line(color="black"),
            axis.title=element_blank(),
            axis.text=element_text(size=12),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=15, colour="black")
      )

    # Save file
    ggsave(file="FigureS5l_GATA1.png", width=2, height=4)

#Summary of numbers for figure
table(data$genotype.TP53)
#WT; 329
#Mut; 31

# Statistical test
wilcox.test(y~x)$p.value # 0.0006489584

###########################################################
## CEBPA

# Boxplot
    # Definition
    data <- df
    x <- data$genotype.TP53
    y <- data[["CEBPA"]]
    z <- data$genotype.TP53
    maintitle <- ""
    ytitle <- "log2(FPKM + 1)"
    xtitle <- ""
    #legendtitle <- "CEBPA Genotype"
    #xlabels <- n.cells$label
    #fivenum(y) ; ymin <- -0.5 ; ymax <- 1.0 ; yinterval <- 0.25

    # Plot
    plot <- ggplot() +
      geom_boxplot(data, mapping=aes(x=x, y=y, fill=z), outlier.shape=NA, show.legend = FALSE) +
      scale_fill_manual(values=c("grey60","red3"))+
      #ylim(c(0,10))+
      #scale_y_log10(expand=c(0,0))+
      geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.1,alpha=0.5) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border=element_blank(),
            axis.line.y.left = element_line(color="black"),
            axis.line.x = element_line(color="black"),
            axis.title=element_blank(),
            axis.text=element_text(size=12),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=15, colour="black")
      )

    # Save file
    ggsave(file="FigureS5l_CEBPA.png", width=2, height=4)

#Summary of numbers for figure
table(data$genotype.TP53)
#WT; 329
#Mut; 31

# Statistical test
wilcox.test(y~x)$p.value # 0.04765247

