if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  BiocManager::install("apeglm", update = FALSE)
}

library(magrittr)
library(matrixStats)
library(ggplot2)
library(DESeq2)

cts <- as.matrix(read.table("C:/Users/doron/source/semesterE-UF/Bioinformatics/bioinformatics/GSE119290_Readhead_2018_RNAseq_gene_counts.txt"))
coldata<-read.csv("C:/Users/doron/source/semesterE-UF/Bioinformatics/bioinformatics/coldata.csv",header = T,row.names=1,stringsAsFactors=T)
coldata

dds <- DESeqDataSetFromMatrix(countData = cts,colData = coldata,design = ~ dex)
dds        
head(assay(dds))

nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1,]
nrow(dds)

deseq_object <- DESeq(dds)
deseq_results <- results(deseq_object)
deseq_results <- lfcShrink(deseq_object,coef = 2, res = deseq_results )
head(deseq_results)

deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))
head(deseq_df)

plotCounts(dds, gene = "DDX11L1", intgroup = "dex")

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,lab = deseq_df$Gene,x = "log2FoldChange",
  y = "padj",pCutoff = 0.01 )

volcano_plot

vst <- vst(dds,blind = FALSE)
head(assay(vsd), 3)

sampleDists <- dist(t(assay(vst)))
sampleDists
install.packages("pheatmap")
library ("pheatmap")
library ("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vst$dex, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,col = colors)

topVarGenes <- head(order(rowVars(assay(vst)), decreasing = TRUE), 20)
topVarGenes
mat  <- assay(vst)[topVarGenes, ]
mat<-mat-rowMeans(mat)
anno <- as.data.frame(colData(vst)[c("dex")])
pheatmap(mat, annotation_col = anno)

#df <- read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt")
#df_counts <- read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt")
#df_metadata <- read.table("GSE119290_series_matrix.txt")
#matrix_df_counts <- matrix(as.numeric(unlist(df_counts)),nrow=nrow(df_counts))
#matrix_df_counts
#matrix_df_metadata <- matrix(as.character(unlist(df_metadata)),nrow=nrow(df_metadata))
#matrix_df_metadata 

#matrix_of_data <- matrix(as.numeric(unlist(df)),nrow=nrow(df))
matrix_of_data <- cts
matrix_of_ranges <- rowRanges(matrix_of_data, rows = NULL, cols = NULL, na.rm = FALSE, dim. = dim(matrix_of_data), useNames = NA)

vector_of_ranges <- 0
for(i in 1:26364) {
  vector_of_ranges[i] <- matrix_of_ranges[i,2]-matrix_of_ranges[i,1]
}
vector_of_ranges
plot(density(vector_of_ranges), log='x')

