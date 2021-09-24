library(matrixStats)
library(ggplot2)
library(DESeq2)

df <- read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt")

matrix_of_data <- matrix(as.numeric(unlist(df)),nrow=nrow(df))

matrix_of_ranges <- rowRanges(matrix_of_data, rows = NULL, cols = NULL, na.rm = FALSE, dim. = dim(matrix_of_data), useNames = NA)

vector_of_ranges <- 0
for(i in 1:26364) {
  vector_of_ranges[i] <- matrix_of_ranges[i,2]-matrix_of_ranges[i,1]
}

plot(density(vector_of_ranges), log='x')
