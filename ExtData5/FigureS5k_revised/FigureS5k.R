
###############################################################
########### COMPUTE P53 TARGET GENE SCORES: BEATAML ###########
###############################################################

# Load packages
library(data.table)
library(plyr)
library(ggplot2)

# Read patient-level file
df <- read.table("../../../Sources/BeatAML/Master.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Remove entries with missing expression data
df <- df[which(!is.na(df$erythroid.score)), ]

# Read expression files
    # Matrix
    path <- "../../../Sources/BeatAML/FPKM/"
    file <- "FPKM.txt"
    df.exp <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE))

    row.names(df.exp) <- df.exp$gene_id
    df.exp$gene_id <- NULL
    
    # featureData
    path <- "../../../Sources/BeatAML/FPKM/"
    file <- "FPKM_featureData.txt"
    df.exp.feature <- as.data.frame(fread(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE))

# Transform and censor values
df.exp <- log2(df.exp + 1)
df.exp[df.exp < 1] <- 0

# Define gene set
path <- "../../FigureS4/FigureS4c_revised/"
file <- "Fischer-HSPC MUT Direction.txt"
gene.set <- read.table(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)

gene_short_names <- gene.set[which(gene.set$direction=="Overlap: Fischer up_HSPC MUT down"), "gene_short_name"]
length(gene_short_names)

# Check for non-matches
setdiff(gene_short_names, df.exp.feature$gene_short_name)
# [1] "RARRES3"

grep("PLAAT4", df.exp.feature$gene_short_name, value=TRUE)

gene_short_names <- gsub("RARRES3", "PLAAT4", gene_short_names, fixed=TRUE)

setdiff(gene_short_names, df.exp.feature$gene_short_name) # All genes found

# Subset gene set
df.exp.feature <- df.exp.feature[which(df.exp.feature$gene_short_name %in% gene_short_names), ]
df.exp <- df.exp[df.exp.feature$gene_id, ]

# Compute score
. <- apply(df.exp, 2, function(x) {mean(x)})
. <- data.frame("patient.id"=names(.), "score"=as.numeric(.), stringsAsFactors=FALSE)

# Annotate
df <- join(df, ., by="patient.id", type="left")
sum(is.na(df$score))

names(df)[which(names(df)=="score")] <- "p53.score"

# Save file
write.table(df[,c("patient.id", "p53.score")], "p53 Score_Fisher up_HSPC MUT down.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

###############################################################
########### P53 TARGET GENE SCORES VS ERY-HIGH/LOW ############
###############################################################

# Load packages
library(plyr)
library(ggplot2)

# Read patient-level file
path <- "../../../Sources/BeatAML/"
file <- "Master.txt"
df <- read.table(paste(path, file, sep=""), sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Remove entries with missing expression data
df <- df[which(!is.na(df$erythroid.score)), ]

# Annotate p53 score
    # Read file
    score <- read.table("p53 Score_Fisher up_HSPC MUT down.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

    # Annotate
    df <- join(df, score, by="patient.id", type="left")

#################################################################
######################## erythroid.score ########################
#################################################################

# Binarize ery score: erythroid.greenleaf.score
df$score.bi <- ifelse(df$erythroid.score < median(df$erythroid.score), "low", "high")
df$score.bi <- factor(df$score.bi, levels=c("low", "high"))
table(df$score.bi)
    
# Set factor levels
table(df$score.bi)
df$score.bi <- factor(df$score.bi,
                      levels=c("low", "high")
                      )
                      
# Boxplot
    # Definition
    data <- df
    x <- data$score.bi
    y <- data$p53.score
    z <- data$score.bi
    maintitle <- ""
    xtitle <- ""
    ytitle <- "p53 target gene scores"
    legendtitle <- "Ery"
    
    # Plot
    plot <- ggplot() +
            geom_boxplot(data, mapping=aes(x=x, y=y, fill=z), outlier.shape=NA, show.legend = FALSE) +
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
    ggsave("FigureS5k.png", plot, width=2, height=4)

# Wilcox
wilcox.test(y ~ x)$p.value # 0.0135901

