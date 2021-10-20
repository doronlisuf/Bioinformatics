library("matrixStats")
library(ggplot2)
library(DESeq2)
library(DOSE)
library(gplots)
library(Biobase)
library(pheatmap)
library(mclust)
library(ConsensusClusterPlus)
library(cba)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
#load counts expression data (normalize it)
matrix_of_data <- as.matrix(read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt"))
coldata<-read.csv("coldata.csv",header = T,row.names=1,stringsAsFactors=T)
dds <- DESeqDataSetFromMatrix(countData = matrix_of_data,colData = coldata,design = ~ dex)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#create distance matrix with the top 10 variable genes
mads=apply(matrix_of_data,1,mad)
d=matrix_of_data[rev(order(mads))[1:10],]
d=sweep(d,1, apply(d,1,median,na.rm=T))
d = scale(d)
di <- dist(d)

#perform hclustering and plot dendrogram
hc <- hclust(di, method="ward.D")
plot(hc)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#create distance matrix with the top 100 variable genes
mads=apply(matrix_of_data,1,mad)
d=matrix_of_data[rev(order(mads))[1:100],]
d=sweep(d,1, apply(d,1,median,na.rm=T))
d = scale(d)
di <- dist(d)

#perform hclustering and plot dendrogram
hc <- hclust(di, method="ward.D")
plot(hc)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#create distance matrix with the top 1000 variable genes
mads=apply(matrix_of_data,1,mad)
d=matrix_of_data[rev(order(mads))[1:1000],]
d=sweep(d,1, apply(d,1,median,na.rm=T))
d = scale(d)
di <- dist(d)

#perform hclustering and plot dendrogram
hc <- hclust(di, method="ward.D")
plot(hc)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#create distance matrix with the top 10000 variable genes
mads=apply(matrix_of_data,1,mad)
d=matrix_of_data[rev(order(mads))[1:10000],]
d=sweep(d,1, apply(d,1,median,na.rm=T))
d = scale(d)
di <- dist(d)

#perform hclustering and plot dendrogram
hc <- hclust(di, method="ward.D")
plot(hc)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#create distance matrix with the top 5000 variable genes
mads=apply(matrix_of_data,1,mad)
d=matrix_of_data[rev(order(mads))[1:5000],]
d=sweep(d,1, apply(d,1,median,na.rm=T))
d = scale(d)
di <- dist(d)

#perform hclustering and plot dendrogram
hc <- hclust(di, method="ward.D")
plot(hc)

#draw heatmap with clusters and dendrograms

library(dendextend)
library(ComplexHeatmap)
row_dend = as.dendrogram(hclust(dist(d), method="ward.D"))# row clustering
row_dend = color_branches(row_dend, k = 2)
col_dend = as.dendrogram(hclust(dist(t(d)), method="ward.D")) # column clustering
col_dend = color_branches(col_dend, k = 2)


Heatmap(
          d, 
          name = "expression", 
          column_title = c("Cluster 1", "Cluster 2"),
          row_title = c("Cluster 1", "Cluster 2"),
          row_title_rot = 0,
          row_names_gp = gpar(fontsize = 7),
          cluster_rows = row_dend,
          cluster_columns = col_dend,
          row_split = 2,
          column_split = 2
        )







