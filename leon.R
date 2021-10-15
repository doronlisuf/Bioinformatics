if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ConsensusClusterPlus")

library(ConsensusClusterPlus)

matrix_of_data <- as.matrix(read.table("GSE119290_Readhead_2018_RNAseq_gene_counts.txt"))
coldata<-read.csv("coldata.csv",header = T,row.names=1,stringsAsFactors=T)

mads=apply(matrix_of_data,1,mad)
d=matrix_of_data[rev(order(mads))[1:5000],]
d=sweep(d,1, apply(d,1,median,na.rm=T))

title="ConsensusClusterPlus"
results = ConsensusClusterPlus(d,maxK=6,reps=50,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson")
results[[2]][["consensusTree"]]
icl = calcICL(results,title=title,plot="png")
