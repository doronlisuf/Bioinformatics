library(matrixStats)

df <- read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt")

matrixDF <- matrix(as.numeric(unlist(df)),nrow=nrow(df))
matrixDF
x <- rowRanges(matrixDF, rows = NULL, cols = NULL, na.rm = FALSE, dim. = dim(matrixDF), useNames = NA)

x

hist(matrixDF[1,])


for(i in 1:26000) {
  hist(matrixDF[i,])
}

