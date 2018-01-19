###########################
## RNA seq analysis starting with aligned reads
# Ryan Berger 1-10-18
# Work through RNA seq workflow using HIO RNA seq files that were aligned by David Hill using kallisto
###########################

# Need to import the .tsv output files using tximport

# First set up a gene reference (tx2gene) for annotation from ENSEMBL transcipt IDs to gene IDs
source("https://bioconductor.org/biocLite.R")
biocLite("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
edb
# Create a dataframe of transcripts from ensembl human genome
Tx <- transcripts(edb,
                  columns = c(listColumns(edb , "tx"), "gene_name"),
                  return.type = "DataFrame")
head(Tx)
# assign columns 1 (transcript ID) and 7 (gene ID) to dataframe tx2gene
tx2gene <- Tx[,c(1,7)]
head(tx2gene)

# 
# Assign directory to RNA seq folder with alignment files
dir <- '~/Desktop/RNAseq stuff/RNAseq'

# Load in sample key 
setwd('~/Desktop/RNAseq stuff/RNAseq')
sample_key <- read.csv('sample_key.csv')
sample_key$file_name

# Set up path to read files into tximport
files <- file.path(dir, sample_key$file_name, 'abundance.tsv')
# Add sample IDs to files
names(files) <- sample_key$short_name
files

library(tximport)
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
tail(txi$counts, 100)

# Prep for DESeq2 object
library(DESeq2)

# Create a table for sample names. Define the factor (for HIOs it's what was injected) and set the sample names (in this case used from the sample_key.csv file)
sampleTable <- data.frame(condition = factor(sample_key$code_name))
rownames(sampleTable) <- colnames(txi$counts)

# Create DESeq dataset
dds <- DESeqDataSetFromTximport(txi, sampleTable, design = ~condition)
head(dds)
dds
# dds is now ready for DESeq() see DESeq2 vignette


#############################################
## The reads are now in a DESeq dataset. Proceed to exploratory analysis
#############################################

# There are lots of rows with zero reads
nrow(dds)
# Filter out rows with no counts
dds <- dds[rowSums(counts(dds)) > 1, ]
# 12000 rows have been removed
nrow(dds)

## Perform variance stabilizing transformations using rlog
# Should use argument 'blind = FALSE' but it's giving an error code. Leave it out for now.
# blind = FALSE only works when multiple replicates are there
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

# Build results table. Automatically compares fold change between STM over PBS (?)
res <- results(dds)
res
summary(res)

# Be more strict with results: lower false discovery rate threshold (padj) from 10% to 5%
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

# Raise log2 fold change threshold
resLFC1 <- results(dds, lfcThreshold = 1)
table(resLFC1$padj < 0.1)

# Results of padj < 0.05, log2 fold change > 1
res.sig <- results(dds, lfcThreshold = 1, alpha = 0.05)
# using results() doesn't seem to be working. Do it like this:
# Get genes increasing over 2 fold change, p value < 0.05
incr.sig <- res[res$log2FoldChange > 1 & res$pvalue < 0.05, ]
nrow(incr.sig)
incr.sig.df <- as.data.frame(incr.sig)

# Get genes decreasing > 2 fold change, p value < 0.05
# Giving an error: logical subscript contains NAs
decr.sig <- res.sig[res.sig$log2FoldChange < -1 & res.sig$pvalue < 0.05, ]
nrow(decr.sig)
decr.sig.df <- as.data.frame(decr.sig)


## Plotting results
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup = 'condition')

# MA plot
plotMA(res, ylim = c(-5,5), main = 'HIOs STM vs PBS\n MA plot')
# Label top significant gene on MA plot
topGene <- rownames(res)[which.min(res$padj)]
topGenep <- resOrdered[1,7]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2, pch = 1)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

# Add multiple labels to plot
# This adds labels to first 5 rows of ordered results (min p values)
with(resOrdered[1:5,], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2, pch = 1)
  text(baseMean, log2FoldChange, resOrdered[1:5, 7], pos=4, col="dodgerblue")
}) # pos = 1 (under), 2 (left), 3 (above), 4 (right)

head(res)
res.df <- as.data.frame(res) 

# Label a specific gene(s)
gene <- filter(res.df, symbol == 'MT1G') #'MT1G' 'MT2A'
# Filter one value: symbol == 'X', filter a list: symbol %in% c(a,b,c)
with(gene, {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2, pch = 1)
  text(baseMean, log2FoldChange, gene$symbol, pos=2, col="dodgerblue")
})

## Annotation
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


## Exporting results as csv file
resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file = 'STMoverPBS.csv')

## Gene clustering
library(genefilter)
# Create a list of genes with the top variance across samples (rows)
# The format is a list of numbers corresponding to row numbers in the dataset
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
# Make a matrix of the top genes variance (row ID is ENSG, columns are sample rlog expressions)
mat <- assay(rld)[topVarGenes, ]
# Make expression relative to the mean across all samples. Output is fold above mean for each sample
mat <- mat - rowMeans(mat)
# Get column data (sample IDs) for annotation
anno <- as.data.frame(colData(rld))
# Change labels from ENTREZ IDs to gene symbols
topVarGeneIDs <- res[rownames(mat), 7]
# Heatmap of the 20 genes with the highest variance between samples
pheatmap(mat, annotation_col = anno, labels_row = topVarGeneIDs)

# Make a cluster heatmap of a custom list of genes
gene_list <- c('MT1A', 'MT1E', 'MT1F', 'MT1G', 'MT1H', 'MT1X', 'MT2A')
x <- row.names(subset(res, symbol %in% gene_list))
mat1 <- assay(rld)[x, ]
mat1 <- mat1 - rowMeans(mat1)
IDs <- res[rownames(mat1), 7]
pheatmap(mat1, annotation_col = anno, labels_row = IDs)
head(assay(rld))

## Make a cluster heatmap of genes increased 2FC+, top 25 pvalue < 0.05
a <- incr.sig.df[order(incr.sig.df$pvalue), ]
a <- head(a, 25)
gene_list <- a$symbol
# Remove NAs from vector
gene_list <- gene_list[!is.na(gene_list)]
# Remove duplicates from vector
gene_list <- unique(gene_list)
x <- row.names(subset(res, symbol %in% gene_list))
mat1 <- assay(rld)[x, ]
mat1 <- mat1 - rowMeans(mat1)
IDs <- res[rownames(mat1), 7]
pheatmap(mat1, annotation_col = anno, labels_row = IDs)
head(assay(rld))


### GO analysis
# Need a dataframe with log2fold change and gene ID
geneList <- res.df[,c(2,8)]
head(geneList)
class(geneList) # It's a dataframe
# Select genes that are at least 2 fold changed
gene <- geneList[abs(geneList$log2FoldChange) > 1, ]
gene <- as.character(row.names(gene))
head(gene)
class(gene)

gene.df <- bitr(gene, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
gene <- as.character(gene.df$ENTREZID)
# Run GO analysis
# Input is a character vector of genes
ggo <- groupGO(gene = gene,
               OrgDb = org.Hs.eg.db,
               #keytype = 'ENSEMBL', # specify the format of gene terms in list
               ont = 'MF',  # MF (molecular function), BP (biol process), CC (cell compart)
               level = 3,
               readable = TRUE)
head(ggo)
?groupGO
# Plot results of GO
barplot(ggo, drop = TRUE, showCategory = 20, decreasing = TRUE)
# Problem: barplot is not in order of most enriched terms, always the same order
# decreasing = TRUE doesn't change it

## GO over-representation test
# input:
# gene = character vector of selected genes for analysis (use ENTREZID format)
# universe = character vector of all genes in dataset
ego <- enrichGO(gene = gene,
                universe = gene.a$ENTREZID,
               # keytype = 'ENSEMBL',
                OrgDb = org.Hs.eg.db,
                ont = 'BP',
                pAdjustMethod = 'BH',
                minGSSize = 10,
                pvalueCutoff = 0.01,
                qvalueCutoff = 0.05,
                readable = TRUE)
## Gives error every time: No gene set have size > 2 (whatever you set minGSSize to)
## Fixed error:
## Only input ENTREZID format, not working when you put in ENSEMBL (even if specified in keytype)
## 

# Simplify GO terms in over-representation
ego2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)

# Visualize over representation 
dotplot(ego2, showCategory = 20)

enrichMap(ego2, n = 20) # n = number of top nodes to look at, default is 50 

## categorySize can be scaled by 'pvalue' or 'geneNum'
cnetplot(ego2)
?cnetplot

plotGOgraph(ego)

res.p <- res.df[res.df$pvalue < 0.05,]
nrow(res.p)
