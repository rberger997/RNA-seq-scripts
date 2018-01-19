### Make a volcano plot from RNA seq data

# Load data
setwd("~/Desktop")
res <- read.table('results.txt', header = TRUE)
head(res)

# Data is in a dataframe with Gene name, log2foldchange, pvalue, padj

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main = 'Volcano Plot', xlim=c(-2.5,2.5)))

# Add colored points: red if padj<0.05, orange if log2FC>1, green if both
with(subset(res, padj<0.05), points(log2FoldChange, -log10(pvalue), pch=20, col='red'))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col='green'))
with(subset(res, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col='orange'))

# Label points with the textxy function from the calibrate plot
install.packages('calibrate')
library(calibrate)
with(subset(res, padj<0.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.5) )
