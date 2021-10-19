library("matrixStats")
library(ggplot2)
library(DESeq2)
library(DOSE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#load counts expression data (normalize it)
matrix_of_data <- as.matrix(read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt"))
coldata<-read.csv("coldata.csv",header = T,row.names=1,stringsAsFactors=T)
dds <- DESeqDataSetFromMatrix(countData = matrix_of_data,colData = coldata,design = ~ dex)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
#Subset matrix of data into 5000 most variable genes
#find ranges of each row .3722
matrix_of_ranges <- rowRanges(normalized_counts, rows = NULL, cols = NULL, na.rm = FALSE, dim. = dim(normalized_counts), useNames = NA)
vector_of_ranges <- 0
for(i in 1:26364) {
  vector_of_ranges[i] <- matrix_of_ranges[i,2]-matrix_of_ranges[i,1]
}

#create matrix that stores the indices of the top 5000 variable genes
matrix_of_variable_genes <- 0
count <- 2
for(i in 1:26364) {
  if(vector_of_ranges[i] > 2822){
    matrix_of_variable_genes <- rbind(matrix_of_variable_genes, normalized_counts[i, ])
    row.names(matrix_of_variable_genes)[count] <- row.names(normalized_counts)[i]
    count <- count + 1
  }
}

#remove the first row because I set it as 0
matrix_of_variable_genes <- matrix_of_variable_genes[-1, ]


#create distance matrix with the top 5000 variable genes
d <- dist(matrix_of_variable_genes)

#perform hclustering
hc <- hclust(d, method="ward")
#plot dendrogram
plot(hc)
val <- 0
val <- cutree(hc, k=2)
#tapply(coldata$dex, val, mean)




