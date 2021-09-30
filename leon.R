if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!("org.Hs.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

library(umap)
library(data.table)

df <- transpose(read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt"))
df.umap = umap(df)
df.umap
coldata<-read.csv("coldata.csv",header = T,row.names=1,stringsAsFactors=T)
coldata

xylim <- range(df.umap$layout)
plot(xylim, xylim, type="n")
points(df.umap$layout[,1], df.umap$layout[,2], col=as.integer(coldata[,"dex"]), cex=6, pch=20)

tsne(df, labels=as.factor(coldata[,"dex"]))
as.factor(coldata[,"dex"])

if (!("org.Hs.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}
df <- read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt")
library(DOSE)
library(org.Hs.eg.db)
library(clusterProfiler)
data("geneList")
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

             