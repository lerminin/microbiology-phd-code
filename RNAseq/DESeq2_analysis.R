###############################
##
## RNA-seq analysis with DESeq2
##
###############################

# Loading packages and reads

pkgs <- c("calibrate", "DESeq2", "here", "gplots", "caTools", "hexbin", "pheatmap", "vsn", "RColorBrewer", "ggplot2")
lapply(pkgs, library, character = TRUE)
sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252    LC_MONETARY=English_Canada.1252
# [4] LC_NUMERIC=C                    LC_TIME=English_Canada.1252    
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ggplot2_3.3.5               RColorBrewer_1.1-2          vsn_3.62.0                  pheatmap_1.0.12            
# [5] hexbin_1.28.2               caTools_1.18.2              gplots_3.1.1                here_1.0.1                 
# [9] DESeq2_1.34.0               SummarizedExperiment_1.24.0 Biobase_2.54.0              MatrixGenerics_1.6.0       
# [13] matrixStats_0.61.0          GenomicRanges_1.46.1        GenomeInfoDb_1.30.0         IRanges_2.28.0             
# [17] S4Vectors_0.32.3            BiocGenerics_0.40.0         calibrate_1.7.7             MASS_7.3-54     

# Import data from count table
countdata <- read.table(here("polymicrobial_RNAseq_2020", "salmon_quant", "quant_W3110_reademption.txt"), header=TRUE, row.names=1)
head(countdata)

# Assign condition 
(condition <- factor(c(rep("WT_WT", 2), rep("cib_WT", 2), rep("WT_cirA", 2))))
model.matrix(~ condition)

# To change which condition is the reference condition (default: alphabetically)
condition <- relevel(condition, "WT_WT")
model.matrix(~ condition)

# Analysis with DESeq2 ----------------------------------------------------

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds
all(rownames(coldata) == colnames(condition))

# Run the DESeq pipeline
dds <- DESeq(dds)
dds2 <- estimateSizeFactors(dds)
dds3 <- estimateDispersions(dds2)
dds4 <- nbinomWaldTest(dds3)
res <- results(dds)
res

# To pull out pairwise comparisons
res <- results(dds, contrast = c("condition","cib_WT","WT_cirA"))
table(res$padj<0.05)
mcols(res)$description
metadata(res)$filterThreshold

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld),10)
hist(assay(rld))

# To plot standard deviation of rld assay against normalized count mean using hexbin and vsn:
meanSdPlot(assay(rld))
ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

# Colors for plots below - Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Alternative sample distance heatmap with pheatmap
pheatmap(sampleDists)
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
pheatmap(assay(rld)[select,], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE)

# To see heatmap of most differentially expressed genes across samples:
de = (res$padj < 0.05)
de[is.na(de)]=FALSE
pheatmap(t(assay(rld)[de,]), dendrogram = "colum", xlab = "Genes", ylab = "Samples",
         symkey = FALSE, lhei = c(1,3))

# Principal components analysis and scree plot
## Could do with built-in DESeq2 function:
DESeq2::plotPCA(rld, intgroup="condition")
rv <- rowVars(assay(rld))
select <- order(rv,decreasing=TRUE)[seq_len(min(500,length(rv)))]
pca <- prcomp(t(assay(rld)[select,]))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)
scree_plot = data.frame(percentVar)
scree_plot[,2]<-c(1:6)
colnames(scree_plot)<-c("variance","component_number")
ggplot(scree_plot,mapping=aes(x=component_number, y=variance))+geom_bar(stat="identity")

## Jordyn likes her PCA better:
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()

# Get differential expression results
res <- results(dds)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results_cibWTvsWTcirA_reademption.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering - can't get Jordyn's version to work on my machine
#attr(res, "filterThreshold")
#plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")
metadata(res)$filterThreshold
plot(metadata(res)$filterNumRej, type ="b", xlab="quantiles of baseMean", ylab="number of rejections")
abline(v=metadata(res)$filterTheta)


## MA plot
## Could do with built-in DESeq2 function:
DESeq2::plotMA(dds, ylim=c(-2,2))
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot_cibWTvsWTcirA.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot_cibWTvsWTcirA.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()


