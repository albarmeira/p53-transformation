#FigureS6c - p53LSC signature on total TCGA dataset
#Authors: Sean Wen and Alba Meira

# Load packages
library(survival)
library(survminer)
library(ggplot2)

# Read patient-level file
df <- read.table("../../../Sources/TCGA/Master.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Remove entries with missing survival data
df <- df[which(!is.na(df$time)), ]
df <- df[which(!is.na(df$status)), ]

# Remove entries with missing expression data
df <- df[which(!is.na(df$upP53_LSConly.score)), ]

# Binarize scorse: LSConly signature
df$score.bi <- ifelse(df$upP53_LSConly.score < median(df$upP53_LSConly.score), "low", "high")
df$score.bi <- factor(df$score.bi, levels=c("low", "high"))
table(df$score.bi)

###########################################################

# Cox proportional-hazards model
model <- coxph(Surv(time, status) ~ score.bi, data=df)
summary(model) # HR: 2.30 (95% CI: 1.45-3.65), P: 3.9e-04

# Plot
    # Model
    model <- survfit(Surv(time, status) ~ score.bi, data=df)

    # Plot
    xmin <- 0 ; xmax <- ceiling(max(df$time, na.rm=TRUE))
    plot <- ggsurvplot(model,
                       data=df,
                       risk.table=TRUE,
                       palette=c("grey60","red3"),
                       xlim=c(xmin, xmax),
                       break.x.by=1,
                       fontsize=5,
                       ggtheme = theme_classic2(base_size=20))

    # Save file
    pdf("FigureS6c.pdf")
    print(plot, newpage=FALSE)
    dev.off()
