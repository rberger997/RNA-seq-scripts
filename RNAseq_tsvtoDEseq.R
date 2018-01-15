######################################################
##              RNA seq analysis pt I 
#       Kallisto aligned reads to DESeq dataset
#               Ryan Berger 1-15-18
# 
# This script will convert Kallisto aligned RNA seq files from .tsv format
# to a DEseq Dataset for analysis using the DESeq2 pipeline.
# (This script is designed to work with human sequences).
# 
# Before beginning:
#     Have all files/folders for analysis in ~/Desktop/RNAseq stuff/SamplesForAnalysis
#     Have a file named sample_key.csv in directory:
#         with filenames in column 'file_names'
#         sample IDs in column titled short_name
######################################################

# Need to import the .tsv output files using tximport

# First set up a gene reference (tx2gene) for annotation from ENSEMBL transcipt IDs to gene IDs

# If not installed, run next two lines
# source("https://bioconductor.org/biocLite.R")
# biocLite("EnsDb.Hsapiens.v75")
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


# Assign directory to RNA seq folder with alignment files
dir <- '~/Desktop/RNAseq stuff/SamplesForAnalysis'

# Load in sample key 
setwd('~/Desktop/RNAseq stuff/SamplesForAnalysis')
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
sampleTable <- data.frame(condition = factor(sample_key$short_name))
rownames(sampleTable) <- colnames(txi$counts)

# Create DESeq dataset
dds <- DESeqDataSetFromTximport(txi, sampleTable, design = ~condition)
head(dds)
dds
# dds is now ready for DESeq() see DESeq2 vignette


#############################################
## The reads are now in a DESeq dataset. Proceed to exploratory analysis using 
## RNAseq_DESeq2_analysis.R script
#############################################

