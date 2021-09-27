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
if (!("pheatmap" %in% installed.packages())) {
  BiocManager::install("apeglm", update = FALSE)
}
if (!("gprofiler2" %in% installed.packages())) {
  BiocManager::install("apeglm", update = FALSE)
}
if (!("clusterProfiler" %in% installed.packages())) {
  BiocManager::install("apeglm", update = FALSE)
}
library(clusterProfiler)
library(gprofiler2)
library ("pheatmap")
library ("RColorBrewer")
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
head(assay(vst), 3)

sampleDists <- dist(t(assay(vst)))
sampleDists


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

topVarGenesGO <- head(order(rowVars(assay(vst)), decreasing = TRUE), 100)
topVarGenesGO
mat_100  <- assay(vst)[topVarGenesGO, ]

gostres <- gost(query = rownames(mat_100), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

names(gostres)
head(gostres$result, 6)
names(gostres$meta)
gostplot(gostres, capped = TRUE, interactive = TRUE)
plot <- gostplot(gostres, capped = FALSE, interactive = FALSE)
plot
publish_gosttable(gostres, highlight_terms = gostres$result[c(1:2,10,120),],
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)

matrix_of_data <- cts
matrix_of_ranges <- rowRanges(matrix_of_data, rows = NULL, cols = NULL, na.rm = FALSE, dim. = dim(matrix_of_data), useNames = NA)

vector_of_ranges <- 0
for(i in 1:26364) {
  vector_of_ranges[i] <- matrix_of_ranges[i,2]-matrix_of_ranges[i,1]
}
vector_of_ranges
plot(density(vector_of_ranges), log='x')

