# DESeq2 R script for analyzing sex-biased expression

'''R
# First load in the packages required by this script
library(DESeq2)
#library(apeglm)
#library(writex1)
library(dplyr)

# IMPORTANT! set the working directory to the folder that contains the count matrix (which was generated using kallisto on the command line)
setwd("C:/Users/ypott/Desktop/Evans_Project/muelleri_differential_expression")

# (Optional) Make the working directory into a dataframe and then list all the files contained in that directory
dir <- "C:/Users/ypott/Desktop/Evans_Project/muelleri_differential_expression"
list.files(dir)

# Read in the counts matrix
counts <- read.table("muel_kallisto_countz_.isoform.counts.matrix")
# (Optional) Rename the columns. Make sure that the new names are in the same order as the original names
colnames(counts)[1:10] <- c("tad31","tad32","tad33","tad34","tad35","tad36","tad37","tad38","tad39","tad42")
View(counts)
# Read in a table that contains the sexes of all the samples. Make sure the sample names are in the same order as the counts matrix.
genotypes <- read.csv("sex_genotype_muel.csv", header=T)
# The sexes have to be a factor to be interpreted by DESeq
genotypes$Sex <- as.factor(genotypes$Sex)
View(genotypes)

# Construct the DESeq data set: Round the counts to whole numbers, indicate that the genotypes dataframe contains extra info about the samples, signify the condition you want to analyze (this variable MUST be a column in the genotypes dataframe)
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = genotypes,
                              design = ~ Sex)

# DESeq command does the differential expression analysis, including normalizing counts! This analysis is NOT pre-filtered
dds <- DESeq(dds)
# dim gives the dimensions of the dds object
dim(dds)
# relevel so that the log2FoldChange is relative to females (negative values are female biased)
dds$Sex <- relevel(dds$Sex, ref="F")
res <- results(dds)
# summary gives a quick overview of how many genes are female or male biased (and some other stuff)
summary(res)
View(res)

# Repeat the analysis, but this time pre-filter. Here I'm only keeping counts at least double the size of ? 
keep <- rowSums(counts(dds)) >= 2* length(colnames(dds))
View(keep)
# subset the dds dataframe, retaining only what was stored in 'keep'
dds <- dds[keep,]
# dim gets the dimensions of the dataframe
dim(dds)
# repeat the analysis with the subsetted data frame
dds$Sex <- relevel(dds$Sex, ref="F")
dds <- DESeq(dds)
res <- results(dds)
summary(res)
# Order the results by p-value (most to least significant)
resOrdered <- res[order(res$pvalue),]
summary(resOrdered)
head(resOrdered)
View(resOrdered)

#
plotMA(res, ylim=c(-2,2))

# PCA coloured by sex
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("Sex"))

# Volcano plot
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue')

# Get all male-biased contigs 
malebiased <- resOrdered[resOrdered$log2FoldChange > 0,]
head(malebiased)

# Get all female-biased contigs
femalebiased <- resOrdered[resOrdered$log2FoldChange < 0,]
head(femalebiased)

# Want to get all contigs that have high expression in females but ZERO expression in males
# basemean is the average normalized expression across all samples -> want means for F and M samples
resExpression <- res[order(res$baseMean),]
head(resExpression)

# Group the samples by sex
females <- (counts(dds, normalized=T)[,c('tad31','tad32','tad34','tad35','tad36','tad37','tad42')])
males <- (counts(dds, normalized=T)[,c('tad33','tad38','tad39')])

# Change the object type so that rowSums can be applied
datamales=as.data.frame(males)

# Find rows that ONLY contain zeroes (meaning no expresison in all males)
rowswithzeroes <- rowSums(datamales != 0) == 0
maleszeroexpression <- datamales[rowswithzeroes, ]
contigs <- rownames(maleszeroexpression)
# Keep only the female-biased contigs that have ZERO expression in males
femonlyexpression <- femalebiased[rownames(femalebiased) %in% contigs, ]

#Quickly view the contigs with female-biased expression and zero expression in all males in order of significance
head(femonlyexpression[order(femonlyexpression$pvalue),])

#Check how many of the contigs with female-biased expression have an expression level of ZERO in all males
dim(femonlyexpression)



# "Sanity check" to make sure 
norm <- as.data.frame(counts(dds, normalized=T))
rsquare <- data.frame(matrix(ncol = ncol(norm),
                             nrow = ncol(norm)))
# exact = FALSE was added in order to prevent an error
for(i in 1:(ncol(norm)-1)) {
  for(j in (i+1):ncol(norm)) {
    print(paste(i," ",j))
    x <- cor.test(norm[ , i],
                  norm[ , j],
                  method = 'spearman',
                  exact = FALSE)
    rsquare[i,j] <- x$estimate
  }
}

colnames(rsquare) <- colnames(norm)
rownames(rsquare) <- colnames(norm)
View(rsquare)

# View the dispersion plot
plotDispEsts(dds)
'''
