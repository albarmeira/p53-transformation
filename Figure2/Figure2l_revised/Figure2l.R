# Set working directory

#Figure2L - p53LSC signature on total BeatAML dataset

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
dim(df) #322 with concomitant survival data
table(df$genotype.TP53) #28 TP53-mutant and 294 TP53-WT

# Binarize scorse: LSConly signature
df$score.bi <- ifelse(df$upP53_LSConly.score < median(df$upP53_LSConly.score), "low", "high")
df$score.bi <- factor(df$score.bi, levels=c("low", "high"))
table(df$score.bi) # 161 low 161 high

###########################################################

# Survival analysis
    # Cox proportional-hazards model
    model <- coxph(Surv(time, status) ~ score.bi, data=df)
    summary(model)

    # Retrieve coefficients
    coef <- as.data.frame(summary(model)$coefficient)
    coef.int <- as.data.frame(summary(model)$conf.int[,c(3,4)])
    results <- data.frame("coef"=coef, "coef.95.lower"=coef.int[1,1], "coef.95.upper"=coef.int[2,1])

    # Summarise table
    cols <- c("coef.exp.coef.", "coef.95.lower", "coef.95.upper", "coef.Pr...z..")
    cols.new <- c("hazard_ratio", "hazard_ratio_95_lower",  "hazard_ratio_95_upper", "pval")
    results <- results[,cols]
    names(results) <- cols.new

    print(results) # HR=3.42 (95% Cl:2.48-4.71)
                   # P=5.13-14

# KM curve
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
    pdf("Figure2l.pdf")
    print(plot, newpage=FALSE)
    dev.off()
