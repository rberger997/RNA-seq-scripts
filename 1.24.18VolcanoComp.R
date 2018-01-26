### Compare my RNA seq analysis with David's using volcano plot

setwd("~/Desktop/RNAseq stuff/old files")
res.dh <- read.csv('Styphimurium_over_PBS.csv')
head(res.dh)
# dataframe with genes in column X

### Make a volcano plot from RNA seq data

# Data is in a dataframe with Gene name, log2foldchange, pvalue, padj

# Make a basic volcano plot
with(res.dh, plot(log2FoldChange, -log10(pvalue), pch=20, col = 'gray', main = 'Volcano Plot\n HIOs STM vs PBS (David Analysis)'))

# Add colored points: red if padj<0.05, orange if log2FC>1, green if both
with(subset(res.dh, padj<0.05), points(log2FoldChange, -log10(pvalue), pch=20, col='red'))
with(subset(res.dh, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col='green'))
with(subset(res.dh, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col='orange'))

# Label points with the textxy function from the calibrate plot
install.packages('calibrate')
library(calibrate)
with(subset(res, padj<0.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.5) )

############################## Do the same with HIOs RNA seq data

# Have data in res.df
head(res.df)
res <- res.df
# Do with Davids analysis
#setwd("~/Desktop/RNAseq stuff/old files")
#res <- read.csv('STM over PBS1.csv')
head(res)

# Data is in a dataframe with Gene name, log2foldchange, pvalue, padj

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), 
               pch=20, 
               main = 'HIOs STM/PBS Gene expression\n Volcano Plot' 
               #, xlim=c(-10,10)
))

# Add colored points: red if padj<0.05, orange if log2FC>1, green if both
with(subset(res, padj<0.05), points(log2FoldChange, -log10(pvalue), pch=20, col='red'))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col='green'))
with(subset(res, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col='orange'))

# Label points with the textxy function from the calibrate plot
# install.packages('calibrate')
library(calibrate)
with(subset(res, padj<0.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=symbol, cex=.5) )


