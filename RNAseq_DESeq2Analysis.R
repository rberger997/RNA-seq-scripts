######################################################
##              RNA seq analysis pt II 
#       Exploratory Analysis of the DESeq dataset
#               Ryan Berger 1-15-18
# 
# This script will use the DESeq2 pipeline to explore RNA seq results.
#
# Before beginning:
#     Complete RNAseq_tsvtoDEseq.R 
#     Or have a DESeq2 dataset named 'dds'
#    
######################################################

## Add annotation - symbol and entrezID
library(AnnotationDbi)
library(org.Hs.eg.db)
columns(org.Hs.eg.db)

# Add column for gene name (symbol)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = row.names(res),
                     column = 'SYMBOL',
                     keytype = 'ENSEMBL',
                     multiVals = 'first')
# Add column for gene Entrez ID
res$entrez <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = 'ENTREZID',
                     keytype = 'ENSEMBL',
                     multiVals = 'first')
# Put results in order of the gene with the lowest p value
resOrdered <- res[order(res$pvalue),]
head(resOrdered)


# Filtering of rows without any counts
# There are lots of rows with zero counts
nrow(dds)
# Filter out rows with no counts
dds <- dds[rowSums(counts(dds)) > 1, ]
nrow(dds)

## Examine sample distances and PCA plot

## Perform variance stabilizing transformations using rlog
# Should use argument 'blind = FALSE' but it gives an error code if there aren't multiple replicates.
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

# Assess overall similarity between samples using Sample Distances
sampleDists <- dist(t(assay(rld)))
sampleDists
# visualize distances
library(pheatmap)
library(RColorBrewer)

sampleDistMatrix <- as.matrix(sampleDists)
head(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# PCA plot
plotPCA(rld, intgroup = 'condition')


## Differential expression analysis
dds <- DESeq(dds)

# Build results table. Automatically compares fold change between samples (?)
res <- results(dds)
res
summary(res)

## Be more strict with results: lower false discovery rate threshold (padj) from 10% to 5%
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

# Raise log2 fold change threshold
resLFC1 <- results(dds, lfcThreshold = 1)
table(resLFC1$padj < 0.1)

# Results of padj < 0.05, log2 fold change > 1
res.sig <- results(dds, lfcThreshold = 1, alpha = 0.05)
# using results() doesn't seem to be working. Do it like this:
res.sig <- res.sig[res.sig$log2FoldChange > 1 & res.sig$pvalue < 0.05, ]


## Plotting results
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup = 'condition')

# MA plot
plotMA(res, ylim = c(-5,5))

## Exporting results as csv file
resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file = 'diff_expression_results.csv')


# Optional clean up
# rm(list= ls())