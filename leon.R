if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!("M3C" %in% installed.packages())) {
  BiocManager::install("M3C", update = FALSE)
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
