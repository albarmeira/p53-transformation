################################################
#Figure2E - erythroid signatures from BeatAML dataset
################################################

# Load packages
library(ggplot2)

# Read file
path <- "../../../Sources/BeatAML/"
file <- "Master.txt"
df <- read.table(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)
dim(df)

# Set factor level
table(df$genotype.TP53)
df$genotype.TP53 <- factor(df$genotype.TP53,
                           levels=c("Wildtype", "Mutant (1 hit)", "Mutant (2 hits)"),
                           labels=c("WT", "Mut", "Mut")
)

###########################################################

# Boxplot
    # Definition
    data <- df
    x <- data$genotype.TP53
    y <- data$erythroid.score
    maintitle <- ""
    ytitle <- "Erythroid Score"
    xtitle <- ""

    # Plot
    plot <- ggplot() +
      geom_boxplot(data, mapping=aes(x=x, y=y, fill=x), outlier.shape=NA, show.legend = FALSE) +
      scale_fill_manual(values=c("grey60","red3"))+
      geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.1,alpha=0.5) +
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
    ggsave(file="Figure_2E.png", width=2, height=4)

#Summary of numbers for figure
table(data$genotype.TP53)
#WT; 329
#Mut; 31

# Statistical test
wilcox.test(y~x)$p.value
#0.002663848

