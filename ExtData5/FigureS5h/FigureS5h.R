#FigureS5h - Greenleaf erythroid signatures from TCGA dataset
#Authors: Sean Wen and Alba Meira

# Load packages
library(ggplot2)

# Read file
df <- read.table("../../../Sources/TCGA/Master.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# # Remove entries with missing expression data
df <- df[which(!is.na(df$erythroid.greenleaf.score)), ]
dim(df) #151 patients
table(df$genotype.TP53.final) #140 WT and 11 Mut

# Set factor level
table(df$genotype.TP53.final)
df$genotype.TP53.final <- factor(df$genotype.TP53.final, levels=c("WT", "Mut"))

###########################################################

# Boxplot
    # Definition
    data <- df
    x <- data$genotype.TP53.final
    y <- data$erythroid.greenleaf.score

    #Plot
    plot <- ggplot() +
      geom_boxplot(data, mapping=aes(x=x, y=y, fill=x), outlier.shape=NA, show.legend = FALSE) +
      scale_fill_manual(values=c("grey60","red3"))+
      geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.1,alpha=0.5) +
      #ylim(0,6)+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border=element_blank(),
            axis.line.y.left = element_line(color="black"),
            axis.line.x = element_line(color="black"),
            axis.title=element_text(size=12),
            axis.text=element_text(size=12),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=20, colour="black"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()
      )

    # Save file
    ggsave(file="FigureS5h.png", width=2, height=4)

# Statistical test
wilcox.test(y~x)$p.value #0.0009767291
