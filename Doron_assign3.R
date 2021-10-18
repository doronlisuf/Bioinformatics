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
if (!("M3C" %in% installed.packages())) {
  BiocManager::install("M3C", update = FALSE)
}
if (!("Seurat" %in% installed.packages())) {
  BiocManager::install("Seurat", update = FALSE)
}
#library(M3C)
#library(umap)
#library(data.table)
library(Seurat)
library(clusterProfiler)
library(gprofiler2)
library (pheatmap)
library (RColorBrewer)
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


vst <- vst(dds,blind = FALSE)
#plotPCA(vst, intgroup=c("dex"))
head(assay(vst), 3)

topVarGenes <- head(order(rowVars(assay(vst)), decreasing = TRUE), 5000)
topVarGenes
mat_5000  <- assay(vst)[topVarGenes, ]

pbmc <- Seurat::CreateSeuratObject(counts = cts, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
head(pbmc@meta.data, 5)

pbmc <- Seurat::NormalizeData(pbmc)
pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 5000)
head(pbmc@assays, 5)
#top10 <- head(Seurat::VariableFeatures(pbmc), 10)
#plot1 <- Seurat::VariableFeaturePlot(pbmc)
#plot1
#plot2 <- Seurat::LabelPoints(plot = plot1, points = top10)
#plot2
all.genes <- rownames(pbmc)
pbmc <- Seurat::ScaleData(pbmc, features = all.genes)
pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(pbmc), npcs =43)
Seurat::DimPlot(pbmc, reduction = "pca")
pbmc
Seurat::ElbowPlot(pbmc)

pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:20)

pbmc <- Seurat::FindClusters(pbmc, resolution = 1.)
head(Seurat::Idents(pbmc),42)
pbmc <- Seurat::RunUMAP(pbmc, dims = 1:20)
Seurat::DimPlot(pbmc, reduction = "umap")



x <- c(100,1000,10000)
for (val in x) {
  pbmc <- Seurat::CreateSeuratObject(counts = cts, project = "pbmc3k", min.cells = 3, min.features = 200)
  pbmc <- Seurat::NormalizeData(pbmc)
  pbmc <- Seurat::FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = val)
  all.genes <- rownames(pbmc)
  pbmc <- Seurat::ScaleData(pbmc, features = all.genes)
  pbmc <- Seurat::RunPCA(pbmc, features = Seurat::VariableFeatures(pbmc), npcs =43)
  Seurat::DimPlot(pbmc, reduction = "pca")
  Seurat::ElbowPlot(pbmc)
  pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:20)
  pbmc <- Seurat::FindClusters(pbmc, resolution = 1.0)
  head(Seurat::Idents(pbmc),42)
  pbmc <- Seurat::RunUMAP(pbmc, dims = 1:20)
  Seurat::DimPlot(pbmc, reduction = "umap")
}
print(count)

pbmc.markers <- Seurat::FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)
top <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
Seurat::DoHeatmap(pbmc, features = top$gene) + NoLegend()


Seurat::Idents(pbmc)
coldata
chisq <- chisq.test(table(Seurat::Idents(pbmc), coldata$dex))
chisq$p.value
p.adjust(chisq$p.value,method="BH")

#mat<-mat-rowMeans(mat)
#anno <- as.data.frame(colData(vst)[c("dex")])
#pheatmap(mat, annotation_col = anno)

