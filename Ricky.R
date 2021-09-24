library(matrixStats)
library(ggplot2)

df <- read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt")

matrixDF <- matrix(as.numeric(unlist(df)),nrow=nrow(df))
matrixDF
r <- rowRanges(matrixDF, rows = NULL, cols = NULL, na.rm = FALSE, dim. = dim(matrixDF), useNames = NA)

ranges <- 0
for(i in 1:26364) {
  ranges[i] <- r[i,2]-r[i,1]
}
hist(ranges)
max(ranges)
plot(density(ranges))

