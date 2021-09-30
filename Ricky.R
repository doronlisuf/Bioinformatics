library("matrixStats")
library(ggplot2)
library(DESeq2)
library(DOSE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install('org.Hs.eg.db')


df <- read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt")
df_counts <- read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt")


matrix_of_data <- matrix(as.numeric(unlist(df)),nrow=nrow(df))

matrix_of_ranges <- rowRanges(matrix_of_data, rows = NULL, cols = NULL, na.rm = FALSE, dim. = dim(matrix_of_data), useNames = NA)

vector_of_ranges <- 0
for(i in 1:26364) {
  vector_of_ranges[i] <- matrix_of_ranges[i,2]-matrix_of_ranges[i,1]
}

plot(density(vector_of_ranges), log='x')

cts <- as.matrix(read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt"))
coldata<-read.csv("coldata.csv",header = T,row.names=1,stringsAsFactors=T)
coldata

dds <- DESeqDataSetFromMatrix(countData = cts,colData = coldata,design = ~ dex)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("dex"))

#get differentially expressed genes
res <- results(DESeq(dds))
#store values
geneList = res[,2]
#create list of genes
names_list = rownames(res)

library('org.Hs.eg.db')
library(data.table)
#lead gene names into symbols variable
symbols <- rownames(res)
#map symbols to IDs
gene_ids <- mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')

#update ID's in list of gene names
for(i in 1:26364){
  names_list[i] = gene_ids[names_list[i]]
}
names(geneList) = names_list
geneList = sort(geneList, decreasing = TRUE)
gene <- names(geneList)[abs(geneList) > 1.5]
x <- enrichDO(gene          = gene,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = names(geneList),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = TRUE)
output_table <- x@result
setDT(output_table)
output_table
