#FigureS6f - LSC17, p53all and p53LSC signature in all BeatAML patients
#Authors: Sean Wen and Alba Meira

# Load packages
library(survival)
library(survminer)
library(ggplot2)

# Read patient-level file
df <- read.table("../../../Sources/BeatAML/Master.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
dim(df) #360 patients; only including de novo and sAML

# Remove entries with missing survival data
df <- df[which(!is.na(df$time)), ]
df <- df[which(!is.na(df$status)), ]
dim(df) #322 with conconmitant survival data
table(df$genotype.TP53) #28 TP53-mutant and 294 TP53-WT

# Define signatures to analyse
sigs <- c("LSC17.score", "upP53_LSCery.score", "upP53_LSConly.score")

# Cox proportional-hazards model
.list <- list()

for(i in 1:length(sigs)) {
    
    # Define sig
    sig <- sigs[i]
    
    # Subset sig
    cols <- c("time", "status", sig)
    df.small <- df[,cols]
    
    # Binarize scorse
    df.small$score.bi <- ifelse(df.small[,sig] < median(df.small[,sig]), "low", "high")
    df.small$score.bi <- factor(df.small$score.bi, levels=c("low", "high"))
    print(table(df.small$score.bi))
    
    # Build model
    model <- coxph(Surv(time, status) ~ score.bi, data=df.small)

    # Retrieve coefficients
    coef <- as.data.frame(summary(model)$coefficient)
    coef.int <- as.data.frame(summary(model)$conf.int[,c(3,4)])
    results <- data.frame("coef"=coef, "coef.95.lower"=coef.int[1,1], "coef.95.upper"=coef.int[2,1])

    # Summarise table
    cols <- c("coef.exp.coef.", "coef.95.lower", "coef.95.upper", "coef.Pr...z..")
    cols.new <- c("hazard_ratio", "hazard_ratio_95_lower",  "hazard_ratio_95_upper", "pval")
    results <- results[,cols]
    names(results) <- cols.new

    # Indicate sig
    . <- data.frame("gene_signature"=sig)
    results <- cbind.data.frame(., results)
    
    # Save into list
    .list[[i]] <- results
    
}

results <- do.call(rbind.data.frame, .list)

# Save file
write.table(results, "FigureS6f.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Set factor levels (for ordering y-axis)
levels <- rev(sigs)
results$gene_signature <- factor(results$gene_signature, levels=levels)

# Forest plot
    # Definition
    data <- results
    y <- data$gene_signature
    x <- data$hazard_ratio
    x.lower <- data$hazard_ratio_95_lower
    x.upper <- data$hazard_ratio_95_upper
    z <- data$gene_signature
    color.breaks <- rev(c("royalblue", "orangered", "brown"))

    # Plot
    plot <- ggplot() +
            geom_errorbarh(data, mapping=aes(y=y, xmin=x.lower, xmax=x.upper), height=0.25) +
            geom_point(data, mapping=aes(x=x, y=y, fill=z), shape=22, size=5) +
            geom_vline(xintercept=1, color="black", linetype="dashed", size=0.5, alpha=1) +
            scale_fill_manual(values=color.breaks) +
            scale_x_continuous(breaks=seq(from=0,to=ceiling(max(x.upper)),by=1),
                               limits=c(0, ceiling(max(x.upper)))
                               ) +
            xlab("Hazard ratio") +
            ylab(" ") +
            theme_bw() +
            theme(panel.border = element_blank(),
                  panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black"),
                  axis.text.y = element_text(size = 12, colour = "black"),
                  axis.text.x.bottom = element_text(size = 12, colour = "black"),
                  axis.title.x = element_text(size = 12, colour = "black"),
                  legend.position="none"
                  )

    # Save file
    ggsave("FigureS6f.png", plot, width=4.5, height=2.5)
