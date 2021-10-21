library(BS831)
library(Biobase)
library(mclust)
#source("http://bioconductor.org/biocLite.R")
#biocLite("ConsensusClusterPlus")
library(ConsensusClusterPlus)
library(cba)
library(pheatmap)

data(AEDAT.collapsed.mad4k)

eset <- AEDAT.collapsed.mad4k
table(eset$Characteristics.DiseaseState)

## restrict to 3 classes
samples.keep <- eset$Characteristics.DiseaseState %in%
  c("non-basal-like","sporadic basal-like","normal")
cancerSet <- eset[,samples.keep]
pheno <- cancerSet$Characteristics.DiseaseState
table(pheno)

## hclustering + heatmap
assayData <- Biobase::exprs(cancerSet)

## use Eculidean distance for columns/samples
## use ward as agglomeration rule
hc01.col <- hcopt(dist(t(assayData)),method="ward.D")

## use 1-correlation as distance for for rows/genes
## use ward as agglomeration rule
hc01.row <- hcopt(as.dist(1-cor(t(assayData))),method="ward.D")

## making heatmap
annot <- data.frame(as.factor(cancerSet$Characteristics.DiseaseState)) 
rownames(annot) <- colnames(assayData)
colnames(annot) <- c("DiseaseState")

annotCol <- list(DiseaseState = c("green", "orange", "purple"))
names(annotCol$DiseaseState) <- levels(annot$DiseaseState)

CCout <- ConsensusClusterPlus(Biobase::exprs(cancerSet),maxK=6,reps=50,pItem=0.8,pFeature=1,
                              innerLinkage="ward.D", finalLinkage="ward.D",
                              title="cc",clusterAlg="hc",distance="euclidean",seed=1262118388.71279,
                              plot=NULL)

nc <- 3

## remake heatmap, include both subtype and cluster assignments for visual comparison
annot1 <- data.frame(annot,cluster=CCout[[nc]]$consensusClass)
annotCol$cluster <- rainbow(n=nc)
names(annotCol$cluster) <- unique(annot1$cluster)

## use cluster tree from consensus clustering for column ordering in heatmap
clust.col <- CCout[[nc]]$consensusTree
## determine row ordering based on de-novo clustering
clust.row <- hcopt(as.dist(1-cor(t(exprs(cancerSet)))),method="ward.D")

##featureNames(cancerSet) <- fData(cancerSet)$hgnc_symbol
heatmaptitle <- "Heatmap consensus clustering assignment vs. subtype\n top 4k z-score normalized"
pheatmap(d1,
         title="consensus cluster plus with PAM for 5000 genes",
         cluster_rows=hclust(dist(t(d1)), method="ward.D"),
         cluster_cols=results1[[2]]$consensusTree,
         show_rownames = FALSE,
         show_colnames = FALSE,
         scale = "row")
