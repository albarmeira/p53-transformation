#Figure2K - Before and after lasso for p53-LSConly

# Read file
    file <- "upP53_LSConly.txt"
    df.1 <- read.table(file, sep=" ", header=TRUE, row.names=1, stringsAsFactors=FALSE)

    # After lasso
    file <- "upP53_LSConly_Coefficient.txt"
    df.2 <- read.table(file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Tabulate n genes
    # Before lasso
    gene_short_names.before <- df.1$x

    # After lasso
    gene_short_names.after <- df.2$gene_short_name

# Save files
before <- paste(length(gene_short_names.before), " genes before lasso regression", sep="")
after <- paste(length(gene_short_names.after), " genes after lasso regression", sep="")
results <- data.frame("V1"=c(before, after))
write.table(results, "Figure2k.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
