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
res <- results(DESeq(dds))
#store values
geneList = res[,2]
#create list of genes
names_list = rownames(res)

library('org.Hs.eg.db')

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
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
symbols = rownames(res)
mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
             