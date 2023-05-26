#FigureS6d - p53LSC signature on TP53-WT BeatAML dataset; plot figure
#Authors: Wei Wen and Alba Meira

# Load packages
library(survival)
library(survminer)
library(ggplot2)

# Read patient-level file
df <- read.table("../../../Sources/BeatAML/Master.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Remove entries with missing survival data
df <- df[which(!is.na(df$time)), ]
df <- df[which(!is.na(df$status)), ]

# Binarize scorse
df$score.bi <- ifelse(df$upP53_LSConly.score < median(df$upP53_LSConly.score), "low", "high")
df$score.bi <- factor(df$score.bi, levels=c("low", "high"))
table(df$score.bi)

# Select TP53-WT samples
table(df$genotype.TP53)
df <- subset(df,genotype.TP53=="Wildtype")
dim(df)

###########################################################

# Cox proportional-hazards model
model <- coxph(Surv(time, status) ~ score.bi, data=df)
summary(model) # HR: 3.10 (95% CI: 2.22-4.34), P:3.35e-11

# Plot
    # Retrieve sample size
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

    pdf("FigureS6d.pdf")
    print(plot, newpage=FALSE)
    dev.off()
