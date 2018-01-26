#----------------------------------------------------------
#          RNA seq view pathway interactions
#              using ReactomePA package
#       
#                      Ryan Berger 
#                       1-26-18
#       
#  Purpose: Visualize gene fold changes in pathway of your
#           choice from RNA seq differential expression data. 
#
#  Input: RNA seq differential expression results dataframe.
#  Output: Graphic of gene pathway colored by fold change.
#
#  Source: https://bioconductor.org/packages/release/bioc/html/ReactomePA.html
#----------------------------------------------------------

# install package: 
# source("https://bioconductor.org/biocLite.R")
# biocLite("ReactomePA")

library(ReactomePA)

# Input is the pathway you want, organism, foldchange (a named vector of foldchange/entrezid)

#----------------------------------------------------------
#      Run pathway viewer with RNA seq data

# First need to set RNA seq dataframe into named vector:
# This function quickly converts dataframe to vector for pathway viewer
# Inputs are: rna seq dataframe, column # for log2fc, column # for entrezid
foldchgvec <- function(df, log2FC, entrez){
  y <- df[,c(log2FC, entrez)]
  y <- y[!duplicated(y[,2]), ]
  z <- structure(as.numeric(y[,1]), names = y[,2])
}

# Run function on dataframe
head(res.df) # log2FC is column 2, entrezid is column 8. So input is:
foldChange <- foldchgvec(res.df, 2, 8)

# Test for duplicates. False is good. True gives error in viewPathway()
any(duplicated(names(foldChange)))  

# Run pathway viewer
# Find pathways to look at by going to: https://reactome.org/PathwayBrowser/
# Assign pathway of interest to pathway object (spelling must match!)
pathway <- 'NOD1/2 Signaling Pathway'
library(ReactomePA)
viewPathway(pathway, 
            organism = 'human', 
            readable = TRUE, 
            foldChange = foldChange, 
            col.bin = 25,  # number of color bins
            vertex.label.color = 'black', # Text color
            vertex.label.cex = 1,  # Font size
            legend.x = -1.5,  # legend position (- is left, + is right)
            legend.y = -1,   # legend position (- is down, + is up)
            main = pathway
            )

?viewPathway  # Look at options of the pathway plots (not much there)
?netplot  # Use this to change the format of the plots

#----------------------------------------------------------
#  Example using IL17 pathway
# Set up vector of fold change data: these are made up
il17 <- c(5, 2, -2)
# names of vector are gene entrezIDs: these are IDs for IL17RA, IL17C, IL17A
names(il17) <- c(23765, 27189, 3605)

viewPathway('Interleukin-17 signaling', 
            organism = 'human', 
            readable = TRUE, 
            foldChange = il17)
# Common error: Error: !any(duplicated(names(foldChange))) is not TRUE
# Caused by having duplicate entrezIDs
# Check vector first by using
any(duplicated(names(il17)))  # If false, proceed to viewpathway
# Find pathways to look at by going to: https://reactome.org/PathwayBrowser/

# ----------------------------------------------------------
# #  Original way of making input vector. Was streamlined with foldchgvec() function.
# # command + shift + c to bring back code.
# 
# # Start with dataframe of RNA seq results
# # Need to convert to named vector: names(entrezID), vector(log2FC)
# head(res.df)
# # Take columns 2 (log2Foldchange) and 8 (entrezID) into new dataframe
# a <- res.df[,c(2,8)]
# head(a)
# nrow(a)  # 25079 rows
# # Remove duplicate entrezIDs
# a <- a[!duplicated(a$entrez), ]
# nrow(a)  # 17898 rows
# # Make a vector with values log2foldchange and names of entrezIDs
# head(a)
# b <- a$log2FoldChange
# head(b)
# names(b) <- a$entrez
# any(duplicated(names(b)))  # Test for duplicates. False is good.
