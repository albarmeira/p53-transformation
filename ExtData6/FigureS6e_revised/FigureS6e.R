#FigureS6e - p53LSC signature on TCGA TP53-WT samples
#Authors: Sean Wen and Alba Meira

# Load packages
library(survival)
library(survminer)
library(ggplot2)

# Read patient-level file
df <- read.table("../../../Sources/TCGA/Master.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Remove entries with missing expression data
df <- df[which(!is.na(df$upP53_LSConly.score)), ]
dim(df) #151
table(df$genotype.TP53.final) #140 WT #11 Mut

# Remove entries with missing survival data
df <- df[which(!is.na(df$time)), ]
df <- df[which(!is.na(df$status)), ]
dim(df) #132
table(df$genotype.TP53.final) #124 WT #8 Mut

# Binarize score
df$score.bi <- ifelse(df$upP53_LSConly.score < median(df$upP53_LSConly.score), "low", "high")
df$score.bi <- factor(df$score.bi, levels=c("low", "high"))
table(df$score.bi)

# Remove entries which are TP53-mutant
df <- subset(df,genotype.TP53.final=="WT")
dim(df) #124

###########################################################
####################### ALL AML ###########################
###########################################################

# Cox proportional-hazards model
model <- coxph(Surv(time, status) ~ score.bi, data=df)
summary(model) # HR: 2.06 (95% CI: 1.28-3.32), P: 0.0028

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

    pdf("FigureS6e.pdf")
    print(plot, newpage=FALSE)
    dev.off()
