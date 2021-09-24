library("matrixStats")
library(ggplot2)
install.packages("matrixStats")

df_counts <- read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt")
df_metadata <- read.table("GSE119290_series_matrix.txt")
matrix_df_counts <- matrix(as.numeric(unlist(df_counts)),nrow=nrow(df_counts))
matrix_df_counts
matrix_df_metadata <- matrix(as.character(unlist(df_metadata)),nrow=nrow(df_metadata))
matrix_df_metadata 

r <- rowRanges(matrixDF, rows = NULL, cols = NULL, na.rm = FALSE, dim. = dim(matrixDF), useNames = NA)

ranges <- 0
for(i in 1:26364) {
  ranges[i] <- r[i,2]-r[i,1]
}
hist(ranges)
max(ranges)
plot(density(ranges))

