
# Load packages
library(survival)
library(survminer)
library(ggplot2)
library(ggridges)
library(plyr)

# Read patient-level file
df <- read.table("../../../Sources/TCGA/Master.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Remove entries with missing survival data
df <- df[which(!is.na(df$time)), ]
df <- df[which(!is.na(df$status)), ]

# Remove entries with missing expression data
df <- df[which(!is.na(df$upP53_LSConly.score)), ]

# Remove entries with missing genotype
df <- df[which(!is.na(df$genotype.TP53.final)), ]

# Binarize scorse
df$score.bi <- ifelse(df$upP53_LSConly.score < median(df$upP53_LSConly.score), "low", "high")
df$score.bi <- factor(df$score.bi, levels=c("low", "high"))
table(df$score.bi)

# Set factor levels
table(df$genotype.TP53.final)
df$genotype.TP53.final <- factor(df$genotype.TP53.final, levels=c("WT", "Mut"), labels=c("WT", "MUT"))

###########################################################

# Density plot (x-axis: score)
    # Definition
    data <- df
    x <- data$upP53_LSConly.score
    y <- data$genotype.TP53
    maintitle <- ""
    xtitle <- "p53LSC score"
    ytitle <- ""
    xintercept <- median(x)

    # Plot
    plot <- ggplot() +
        geom_density_ridges(data, mapping=aes(x=x, y=y, fill=y), alpha=0.75, color="black", scale = 1.0) +
        geom_vline(xintercept=xintercept, linetype="dashed", size=0.25, color="black") +
        scale_fill_manual(values=c("grey90", "red3")) +
        labs(title=maintitle, x=xtitle, y=ytitle) +
        theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border=element_blank(),
            plot.title=element_text(hjust = 0.5, size=12),
            plot.subtitle=element_text(hjust = 0.5, size=15),
            axis.line.y.left = element_line(color="black"),
            axis.line.x = element_line(color="black"),
            axis.title=element_text(size=14),
            axis.text.x=element_blank(),
            axis.text.y=element_text(size=12, colour="black"),
            axis.ticks.x=element_blank(),
            legend.position="none",
            legend.title=element_blank(),
            legend.text=element_text(size=10)
            )
            
    # Save plot
    ggsave("FigureS6b.png", plot, width=3, height=2)

