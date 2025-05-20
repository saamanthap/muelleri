# Analysis of masculinization of gene expression (Kallisto, DeSeq2)
### Copied from Ben's repo and annotated by me
```R
library(edgeR)
library(tximport)
library('edgeR')
library('rhdf5')
library('readxl')
library('ggplot2')
library(grid)
require('gridExtra')
library("org.Xl.eg.db")
library(PCAtools)
library("HTSFilter")
library(tidyverse)
library(purrr)
library(writexl)

setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/submitted/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2023_Kallisto_DeSeq2_done")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/submitted/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2023_Kallisto_DeSeq2_done"
list.files(dir)


f_files<- list.files(".", pattern = "Kallisto_DeSeq2_unfiltered.csv", full.names = T);f_files
# import into a list
myfiles = lapply(f_files, read.delim, sep = ",")

# rename the columns so they are sensible
colnames(myfiles[[1]]) <- c("gene","MF_ccdc_baseMean","MF_ccdc_logFC","MF_ccdc_lfcSE","MF_ccdc_stat","MF_ccdc_pvalue","MF_ccdc_padj")
colnames(myfiles[[2]]) <- c("gene","MF_dmrt1L_baseMean","MF_dmrt1L_logFC","MF_dmrt1L_lfcSE","MF_dmrt1L_stat","MF_dmrt1L_pvalue","MF_dmrt1L_padj")
colnames(myfiles[[3]]) <- c("gene","MF_dmrt1S_baseMean","MF_dmrt1S_logFC","MF_dmrt1S_lfcSE","MF_dmrt1S_stat","MF_dmrt1S_pvalue","MF_dmrt1S_padj")
colnames(myfiles[[4]]) <- c("gene","wtko_ccdc_baseMean","wtko_ccdc_logFC","wtko_ccdc_lfcSE","wtko_ccdc_stat","wtko_ccdc_pvalue","wtko_ccdc_padj")
colnames(myfiles[[5]]) <- c("gene","wtko_dmw_baseMean","wtko_dmw_logFC","wtko_dmw_lfcSE","wtko_dmw_stat","wtko_dmw_pvalue","wtko_dmw_padj")
colnames(myfiles[[6]]) <- c("gene","wtko_scan_baseMean","wtko_scan_logFC","wtko_scan_lfcSE","wtko_scan_stat","wtko_scan_pvalue","wtko_scan_padj")

library(plyr)
alldata<-join_all(myfiles, by = "gene", type = "full", match = "all")
library(ggplot2)
library(GGally)

# get rid of outliers
# https://www.r-bloggers.com/2020/01/how-to-remove-outliers-in-r/
boxplot(alldata$MF_ccdc_logFC, plot=FALSE)$out
# these outliers are the first quartile - 1.5 the interquartile range
# and the third quartile plus 1.5 the interquartile range
#MF_ccdc
outliers <- boxplot(alldata$MF_ccdc_logFC, plot=FALSE)$out
MF_ccdc_logFC_trim<-alldata[,c(1,3)]
if(any(outliers)) {
  MF_ccdc_logFC_trim<- MF_ccdc_logFC_trim[-which(MF_ccdc_logFC_trim$MF_ccdc_logFC %in% outliers),]
}
#MF_dmrt1L
outliers <- boxplot(alldata$MF_dmrt1L_logFC, plot=FALSE)$out
MF_dmrt1L_logFC_trim<-alldata[,c(1,9)]
if(any(outliers)) {
  MF_dmrt1L_logFC_trim<- MF_dmrt1L_logFC_trim[-which(MF_dmrt1L_logFC_trim$MF_dmrt1L_logFC %in% outliers),]
}
#MF_dmrt1S
outliers <- boxplot(alldata$MF_dmrt1S_logFC, plot=FALSE)$out
MF_dmrt1S_logFC_trim<-alldata[,c(1,15)]
if(any(outliers)) {
  MF_dmrt1S_logFC_trim<- MF_dmrt1S_logFC_trim[-which(MF_dmrt1S_logFC_trim$MF_dmrt1S_logFC %in% outliers),]
}

#wtko_dmw
outliers <- boxplot(alldata$wtko_dmw_logFC, plot=FALSE)$out
wtko_dmw_logFC_trim<- alldata[,c(1,27)]
if(any(outliers)) {
  wtko_dmw_logFC_trim<- wtko_dmw_logFC_trim[-which(wtko_dmw_logFC_trim$wtko_dmw_logFC %in% outliers),]
}
#wtko_scan
outliers <- boxplot(alldata$wtko_scan_logFC, plot=FALSE)$out
wtko_scan_logFC_trim<- alldata[,c(1,33)]
if(any(outliers)) {
  wtko_scan_logFC_trim<- wtko_scan_logFC_trim[-which(wtko_scan_logFC_trim$wtko_scan_logFC %in% outliers),]
}
#wtko_ccdc
outliers <- boxplot(alldata$wtko_ccdc_logFC, plot=FALSE)$out
wtko_ccdc_logFC_trim<- alldata[,c(1,21)]
if(any(outliers)) {
  wtko_ccdc_logFC_trim<- wtko_ccdc_logFC_trim[-which(wtko_ccdc_logFC_trim$wtko_ccdc_logFC %in% outliers),]
}

# combine no outlier files
# make a list of df
df_list <- list(MF_ccdc_logFC_trim, MF_dmrt1L_logFC_trim, MF_dmrt1S_logFC_trim, wtko_dmw_logFC_trim, wtko_scan_logFC_trim, wtko_ccdc_logFC_trim)
#merge all data frames in list
alldata_no_outliers <- df_list %>% reduce(full_join, by='gene')

colnames(alldata_no_outliers) <- c("gene","MF1","MF2",
                                   "MF3","dmw","scan",
                                   "ccdc")
library(ggplot2)

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
  #  geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}

my_fn2 <- function(data, mapping, method="p", use="pairwise", ...){
  
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue, 
  # zero to white, and one to red 
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  ggally_cor(data = data, mapping = mapping, ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

my_custom_smooth <- function(data, mapping, ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(color = I("blue")) + 
    geom_smooth(method = "lm", fill="blue", color="blue", ...)
  
  lmModel <- eval(substitute(lm(y ~ x, data = data), mapping))
  fs <- summary(lmModel)$fstatistic
  pValue <- pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
  
  if (pValue < 0.05) {
    p <- p + theme(
      panel.border = element_rect(
        color = "red", 
        size = 3,
        linetype = "solid",
        fill = "transparent"
      )
    )
  }
  
  p
}
p_ <- GGally::print_if_interactive
g<-ggpairs(alldata_no_outliers[,c(2:7)], 
        #upper = list(continuous = "density", combo = "box_no_facet"),
        #upper = list(continuous = wrap(ggally_cor, size = 2)), 
        # upper = list(continuous = my_fn2),
        #upper = "blank",
        lower = list(continuous = my_fn)) +
        #lower = list(continuous = my_custom_smooth)) +
  #theme_bw() +
  theme(strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid")) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank());g

box1_2 <- ggally_text("\nr = 0.249\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_3 <- ggally_text("\nr = -0.011\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_4 <- ggally_text("\nr = 0.287*\np = 0.096\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_5 <- ggally_text("\nr = -0.160\np = 0.822\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_6 <- ggally_text("\nr = 0.370*\np = 0.783\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_3 <- ggally_text("\nr = -0.021\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_4 <- ggally_text("\nr = 0.254\np = 0.208\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_5 <- ggally_text("\nr = -0.299*\np = 0.972\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_6 <- ggally_text("\nr = -0.044\np = 0.490\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_4 <- ggally_text("\nr = 0.460*\np = 0.003*\n",geom_text = ggplot2::aes(size = 6), color = I("red"))
box3_5 <- ggally_text("\nr = -0.053\np = 0.649\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_6 <- ggally_text("\nr = -0.320*\np = 0.859\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box4_5 <- ggally_text("\nr = -0.050\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box4_6 <- ggally_text("\nr = 0.158\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box5_6 <- ggally_text("\nr = -0.021\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
g[1, 2] <- box1_2
g[1, 3] <- box1_3
g[1, 4] <- box1_4
g[1, 5] <- box1_5
g[1, 6] <- box1_6
g[2, 3] <- box2_3
g[2, 4] <- box2_4
g[2, 5] <- box2_5
g[2, 6] <- box2_6
g[3, 4] <- box3_4
g[3, 5] <- box3_5
g[3, 6] <- box3_6
g[4, 5] <- box4_5
g[4, 6] <- box4_6
g[5, 6] <- box5_6
# small function to display plots only if it's interactive

p_(g)

ggsave(file="Kallisto_DeSeq2_sexrelated_pairwise_unfiltered_new.pdf", g, width=10, height=4)
```

# Analysis of masculinization of gene expression (Kallisto, EdgeR)
```R
library(edgeR)
library(tximport)
library('edgeR')
library('rhdf5')
library('readxl')
library('ggplot2')
library(grid)
require('gridExtra')
library("org.Xl.eg.db")
library(PCAtools)
library("HTSFilter")
library(tidyverse)
library(purrr)
library(writexl)

setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/submitted/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2023_Kallisto_EdgeR_done")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/submitted/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/2022_Kallisto_EdgeR_done"
list.files(dir)


f_files<- list.files(".", pattern = "Kallisto_edgeR_unfiltered.csv", full.names = T);f_files
# import into a list
myfiles = lapply(f_files, read.delim, sep = ",")

# rename the columns so they are sensible
cnames <- data.frame(MF_ccdc = c("gene","MF_ccdc_logFC","MF_ccdc_logCPM","MF_ccdc_PValue"),
                     MF_dmrt1L = c("gene","MF_dmrt1L_logFC","MF_dmrt1L_logCPM","MF_dmrt1L_PValue"),
                     MF_dmrt1S = c("gene","MF_dmrt1S_logFC","MF_dmrt1S_logCPM","MF_dmrt1S_PValue"),
                     wtko_ccdc = c("gene","wtko_ccdc_logFC","wtko_ccdc_logCPM","wtko_ccdc_PValue"),
                     wtko_dmw = c("gene","wtko_dmw_logFC","wtko_dmw_logCPM","wtko_dmw_PValue"),
                     wtko_scan = c("gene","wtko_scan_logFC","wtko_scan_logCPM","wtko_scan_PValue"))
colnames(myfiles[[1]]) <- c("gene","MF_ccdc_logFC","MF_ccdc_logCPM","MF_ccdc_PValue")
colnames(myfiles[[2]]) <- c("gene","MF_dmrt1L_logFC","MF_dmrt1L_logCPM","MF_dmrt1L_PValue")
colnames(myfiles[[3]]) <- c("gene","MF_dmrt1S_logFC","MF_dmrt1S_logCPM","MF_dmrt1S_PValue")
colnames(myfiles[[4]]) <- c("gene","wtko_ccdc_logFC","wtko_ccdc_logCPM","wtko_ccdc_PValue")
colnames(myfiles[[5]]) <- c("gene","wtko_dmw_logFC","wtko_dmw_logCPM","wtko_dmw_PValue")
colnames(myfiles[[6]]) <- c("gene","wtko_scan_logFC","wtko_scan_logCPM","wtko_scan_PValue")

library(plyr)
alldata<-join_all(myfiles, by = "gene", type = "full", match = "all")
library(ggplot2)
library(GGally)

# get rid of outliers
# https://www.r-bloggers.com/2020/01/how-to-remove-outliers-in-r/
boxplot(alldata$MF_ccdc_logFC, plot=FALSE)$out
# these outliers are the first quartile - 1.5 the interquartile range
# and the third quartile plus 1.5 the interquartile range
#MF_ccdc
outliers <- boxplot(alldata$MF_ccdc_logFC, plot=FALSE)$out
MF_ccdc_logFC_trim<-alldata[,c(1,2)]
if(any(outliers)) {
  MF_ccdc_logFC_trim<- MF_ccdc_logFC_trim[-which(MF_ccdc_logFC_trim$MF_ccdc_logFC %in% outliers),]
}
#MF_dmrt1L
outliers <- boxplot(alldata$MF_dmrt1L_logFC, plot=FALSE)$out
MF_dmrt1L_logFC_trim<-alldata[,c(1,5)]
if(any(outliers)) {
  MF_dmrt1L_logFC_trim<- MF_dmrt1L_logFC_trim[-which(MF_dmrt1L_logFC_trim$MF_dmrt1L_logFC %in% outliers),]
}
#MF_dmrt1S
outliers <- boxplot(alldata$MF_dmrt1S_logFC, plot=FALSE)$out
MF_dmrt1S_logFC_trim<-alldata[,c(1,8)]
if(any(outliers)) {
  MF_dmrt1S_logFC_trim<- MF_dmrt1S_logFC_trim[-which(MF_dmrt1S_logFC_trim$MF_dmrt1S_logFC %in% outliers),]
}
#wtko_dmw
outliers <- boxplot(alldata$wtko_dmw_logFC, plot=FALSE)$out
wtko_dmw_logFC_trim<- alldata[,c(1,14)]
if(any(outliers)) {
  wtko_dmw_logFC_trim<- wtko_dmw_logFC_trim[-which(wtko_dmw_logFC_trim$wtko_dmw_logFC %in% outliers),]
}
#wtko_scan
outliers <- boxplot(alldata$wtko_scan_logFC, plot=FALSE)$out
wtko_scan_logFC_trim<- alldata[,c(1,17)]
if(any(outliers)) {
  wtko_scan_logFC_trim<- wtko_scan_logFC_trim[-which(wtko_scan_logFC_trim$wtko_scan_logFC %in% outliers),]
}
#wtko_ccdc
outliers <- boxplot(alldata$wtko_ccdc_logFC, plot=FALSE)$out
wtko_ccdc_logFC_trim<- alldata[,c(1,11)]
if(any(outliers)) {
  wtko_ccdc_logFC_trim<- wtko_ccdc_logFC_trim[-which(wtko_ccdc_logFC_trim$wtko_ccdc_logFC %in% outliers),]
}

# combine no outlier files
# make a list of df
df_list <- list(MF_ccdc_logFC_trim, MF_dmrt1L_logFC_trim, MF_dmrt1S_logFC_trim, wtko_dmw_logFC_trim, wtko_scan_logFC_trim, wtko_ccdc_logFC_trim)
#merge all data frames in list
alldata_no_outliers <- df_list %>% reduce(full_join, by='gene')

colnames(alldata_no_outliers) <- c("gene","MF1","MF2",
                                   "MF3","dmw","scan",
                                   "ccdc")
library(ggplot2)

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
  #  geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}

my_fn2 <- function(data, mapping, method="p", use="pairwise", ...){
  
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue, 
  # zero to white, and one to red 
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  ggally_cor(data = data, mapping = mapping, ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

my_custom_smooth <- function(data, mapping, ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(color = I("blue")) + 
    geom_smooth(method = "lm", fill="blue", color="blue", ...)
  
  lmModel <- eval(substitute(lm(y ~ x, data = data), mapping))
  fs <- summary(lmModel)$fstatistic
  pValue <- pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
  
  if (pValue < 0.05) {
    p <- p + theme(
      panel.border = element_rect(
        color = "red", 
        size = 3,
        linetype = "solid",
        fill = "transparent"
      )
    )
  }
  
  p
}

p_ <- GGally::print_if_interactive
g<-ggpairs(alldata_no_outliers[,c(2:7)], 
        #upper = list(continuous = "density", combo = "box_no_facet"),
        #upper = list(continuous = wrap(ggally_cor, size = 2)), 
        upper = list(continuous = my_fn2),
        lower = list(continuous = my_fn)) +
        #lower = list(continuous = my_custom_smooth)) +
  #theme_bw() +
  theme(strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid")) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

box1_2 <- ggally_text("\nr = 0.265\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_3 <- ggally_text("\nr = -0.195\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_4 <- ggally_text("\nr = 0.307*\np = 0.086\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_5 <- ggally_text("\nr = 0.099\np = 0.268\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_6 <- ggally_text("\nr = 0.409*\np = 0.705\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_3 <- ggally_text("\nr = -0.156\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_4 <- ggally_text("\nr = 0.213\np = 0.272\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_5 <- ggally_text("\nr = 0.254\np = 0.048*\n",geom_text = ggplot2::aes(size = 6), color = I("red"))
box2_6 <- ggally_text("\nr = -0.126\np = 0.667\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_4 <- ggally_text("\nr = 0.196\np = 0.108\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_5 <- ggally_text("\nr = 0.028\np = 0.364\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_6 <- ggally_text("\nr = -0.231\np = 0.690\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box4_5 <- ggally_text("\nr = 0.054\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box4_6 <- ggally_text("\nr = 0.097\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box5_6 <- ggally_text("\nr = -0.012\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
g[1, 2] <- box1_2
g[1, 3] <- box1_3
g[1, 4] <- box1_4
g[1, 5] <- box1_5
g[1, 6] <- box1_6
g[2, 3] <- box2_3
g[2, 4] <- box2_4
g[2, 5] <- box2_5
g[2, 6] <- box2_6
g[3, 4] <- box3_4
g[3, 5] <- box3_5
g[3, 6] <- box3_6
g[4, 5] <- box4_5
g[4, 6] <- box4_6
g[5, 6] <- box5_6
# small function to display plots only if it's interactive

p_(g)

ggsave(file="Kallisto_edgeR_sexrelated_pairwise_unfiltered.pdf", g, width=10, height=4)

```
# Analysis of masculinization of gene expression (STAR, DeSeq2)
```R
library(edgeR)
library(tximport)
library('edgeR')
library('rhdf5')
library('readxl')
library('ggplot2')
library(grid)
require('gridExtra')
library("org.Xl.eg.db")
library(PCAtools)
library("HTSFilter")
library(tidyverse)
library(purrr)
library(writexl)

setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/submitted/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR_done")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/submitted/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR_done"
list.files(dir)

# import into a list
temp = list.files(pattern="STAR_DeSeq2_unfiltered.csv");temp
myfiles = lapply(temp, read.delim, sep = ",")

# rename the columns so they are sensible
#cnames <- data.frame(MF_ccdc = c("gene","MF_ccdc_logFC","MF_ccdc_logCPM","MF_ccdc_PValue"),
#                     MF_dmrt1L = c("gene","MF_dmrt1L_logFC","MF_dmrt1L_logCPM","MF_dmrt1L_PValue"),
#                     MF_dmrt1S = c("gene","MF_dmrt1S_logFC","MF_dmrt1S_logCPM","MF_dmrt1S_PValue"),
#                     wtko_ccdc = c("gene","wtko_ccdc_logFC","wtko_ccdc_logCPM","wtko_ccdc_PValue"),
#                     wtko_dmw = c("gene","wtko_dmw_logFC","wtko_dmw_logCPM","wtko_dmw_PValue"),
#                     wtko_scan = c("gene","wtko_scan_logFC","wtko_scan_logCPM","wtko_scan_PValue"))
colnames(myfiles[[1]]) <- c("gene","MF_ccdc_baseMean","MF_ccdc_logFC","MF_ccdc_lfcSE","MF_ccdc_stat","MF_ccdc_pvalue","MF_ccdc_padj")
colnames(myfiles[[2]]) <- c("gene","MF_dmrt1L_baseMean","MF_dmrt1L_logFC","MF_dmrt1L_lfcSE","MF_dmrt1L_stat","MF_dmrt1L_pvalue","MF_dmrt1L_padj")
colnames(myfiles[[3]]) <- c("gene","MF_dmrt1S_baseMean","MF_dmrt1S_logFC","MF_dmrt1S_lfcSE","MF_dmrt1S_stat","MF_dmrt1S_pvalue","MF_dmrt1S_padj")
colnames(myfiles[[4]]) <- c("gene","wtko_ccdc_baseMean","wtko_ccdc_logFC","wtko_ccdc_lfcSE","wtko_ccdc_stat","wtko_ccdc_pvalue","wtko_ccdc_padj")
colnames(myfiles[[5]]) <- c("gene","wtko_dmw_baseMean","wtko_dmw_logFC","wtko_dmw_lfcSE","wtko_dmw_stat","wtko_dmw_pvalue","wtko_dmw_padj")
colnames(myfiles[[6]]) <- c("gene","wtko_scan_baseMean","wtko_scan_logFC","wtko_scan_lfcSE","wtko_scan_stat","wtko_scan_pvalue","wtko_scan_padj")

library(plyr)
alldata<-join_all(myfiles, by = "gene", type = "full", match = "all")
library(ggplot2)
library(GGally)

# get rid of outliers
# https://www.r-bloggers.com/2020/01/how-to-remove-outliers-in-r/
boxplot(alldata$MF_ccdc_logFC, plot=FALSE)$out
# these outliers are the first quartile - 1.5 the interquartile range
# and the third quartile plus 1.5 the interquartile range
#MF_ccdc
outliers <- boxplot(alldata$MF_ccdc_logFC, plot=FALSE)$out
MF_ccdc_logFC_trim<-alldata[,c(1,3)]
MF_ccdc_logFC_trim<- MF_ccdc_logFC_trim[-which(MF_ccdc_logFC_trim$MF_ccdc_logFC %in% outliers),]
#MF_dmrt1L
outliers <- boxplot(alldata$MF_dmrt1L_logFC, plot=FALSE)$out
MF_dmrt1L_logFC_trim<-alldata[,c(1,9)]
MF_dmrt1L_logFC_trim<- MF_dmrt1L_logFC_trim[-which(MF_dmrt1L_logFC_trim$MF_dmrt1L_logFC %in% outliers),]
#MF_dmrt1S
outliers <- boxplot(alldata$MF_dmrt1S_logFC, plot=FALSE)$out
MF_dmrt1S_logFC_trim<-alldata[,c(1,15)]
MF_dmrt1S_logFC_trim<- MF_dmrt1S_logFC_trim[-which(MF_dmrt1S_logFC_trim$MF_dmrt1S_logFC %in% outliers),]
#wtko_dmw
outliers <- boxplot(alldata$wtko_dmw_logFC, plot=FALSE)$out
wtko_dmw_logFC_trim<- alldata[,c(1,27)]
wtko_dmw_logFC_trim<- wtko_dmw_logFC_trim[-which(wtko_dmw_logFC_trim$wtko_dmw_logFC %in% outliers),]
#wtko_scan
outliers <- boxplot(alldata$wtko_scan_logFC, plot=FALSE)$out
wtko_scan_logFC_trim<- alldata[,c(1,33)]
wtko_scan_logFC_trim<- wtko_scan_logFC_trim[-which(wtko_scan_logFC_trim$wtko_scan_logFC %in% outliers),]
#wtko_ccdc
outliers <- boxplot(alldata$wtko_ccdc_logFC, plot=FALSE)$out
wtko_ccdc_logFC_trim<- alldata[,c(1,21)]
wtko_ccdc_logFC_trim<- wtko_ccdc_logFC_trim[-which(wtko_ccdc_logFC_trim$wtko_ccdc_logFC %in% outliers),]


# combine no outlier files
# make a list of df
df_list <- list(MF_ccdc_logFC_trim, MF_dmrt1L_logFC_trim, MF_dmrt1S_logFC_trim, wtko_dmw_logFC_trim, wtko_scan_logFC_trim, wtko_ccdc_logFC_trim)
#merge all data frames in list
alldata_no_outliers <- df_list %>% reduce(full_join, by='gene')

colnames(alldata_no_outliers) <- c("gene","MF1","MF2",
                                   "MF3","dmw","scan",
                                   "ccdc")
library(ggplot2)

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
  #  geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}

my_fn2 <- function(data, mapping, method="p", use="pairwise", ...){
  
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue, 
  # zero to white, and one to red 
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  ggally_cor(data = data, mapping = mapping, ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

my_custom_smooth <- function(data, mapping, ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(color = I("blue")) + 
    geom_smooth(method = "lm", fill="blue", color="blue", ...)
  
  lmModel <- eval(substitute(lm(y ~ x, data = data), mapping))
  fs <- summary(lmModel)$fstatistic
  pValue <- pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
  
  if (pValue < 0.05) {
    p <- p + theme(
      panel.border = element_rect(
        color = "red", 
        size = 3,
        linetype = "solid",
        fill = "transparent"
      )
    )
  }
  
  p
}

p_ <- GGally::print_if_interactive
g<-ggpairs(alldata_no_outliers[,c(2:7)], 
        #upper = list(continuous = "density", combo = "box_no_facet"),
        #upper = list(continuous = wrap(ggally_cor, size = 2)), 
        upper = list(continuous = my_fn2),
        lower = list(continuous = my_fn)) +
        #lower = list(continuous = my_custom_smooth)) +
  #theme_bw() +
  theme(strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid")) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

box1_2 <- ggally_text("\nr = 0.001\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_3 <- ggally_text("\nr = -0.049\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_4 <- ggally_text("\nr = 0.118\np = 0.356\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_5 <- ggally_text("\nr = 0.132\np = 0.258\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_6 <- ggally_text("\nr = 0.176\np = 0.962\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_3 <- ggally_text("\nr = -0.033\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_4 <- ggally_text("\nr = 0.059\np = 0.587\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_5 <- ggally_text("\nr = 0.130\np = 0.218\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_6 <- ggally_text("\nr = -0.083\np = 0.608\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_4 <- ggally_text("\nr = 0.409*\np = 0.011*\n",geom_text = ggplot2::aes(size = 6), color = I("red"))
box3_5 <- ggally_text("\nr = -0.068\np = 0.558\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_6 <- ggally_text("\nr = -0.140\np = 0.513\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box4_5 <- ggally_text("\nr = 0.136\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box4_6 <- ggally_text("\nr = 0.118\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box5_6 <- ggally_text("\nr = -0.163\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
g[1, 2] <- box1_2
g[1, 3] <- box1_3
g[1, 4] <- box1_4
g[1, 5] <- box1_5
g[1, 6] <- box1_6
g[2, 3] <- box2_3
g[2, 4] <- box2_4
g[2, 5] <- box2_5
g[2, 6] <- box2_6
g[3, 4] <- box3_4
g[3, 5] <- box3_5
g[3, 6] <- box3_6
g[4, 5] <- box4_5
g[4, 6] <- box4_6
g[5, 6] <- box5_6
# small function to display plots only if it's interactive

p_(g)

ggsave(file="STAR_DeSeq2_sexrelated_pairwise_unfiltered.pdf", g, width=10, height=4)

```

# Analysis of masculinization of gene expression (STAR, EdgeR)
```R

library(edgeR)
library(tximport)
library('edgeR')
library('rhdf5')
library('readxl')
library('ggplot2')
library(grid)
require('gridExtra')
library("org.Xl.eg.db")
library(PCAtools)
library("HTSFilter")
library(tidyverse)
library(purrr)
library(writexl)

setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/submitted/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR_done")
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/submitted/2022_Supergene/2022_KO_tad_RNAseq/2022_EdgeR_and_DeSeq2/STAR_done"
list.files(dir)


# import into a list
temp = list.files(pattern="STAR_edgeR_unfiltered.csv");temp
myfiles = lapply(temp, read.delim, sep = ",")

# rename the columns so they are sensible
colnames(myfiles[[1]]) <- c("gene","MF_ccdc_logFC","MF_ccdc_logCPM","MF_ccdc_PValue")
colnames(myfiles[[2]]) <- c("gene","MF_dmrt1L_logFC","MF_dmrt1L_logCPM","MF_dmrt1L_PValue")
colnames(myfiles[[3]]) <- c("gene","MF_dmrt1S_logFC","MF_dmrt1S_logCPM","MF_dmrt1S_PValue")
colnames(myfiles[[4]]) <- c("gene","wtko_ccdc_logFC","wtko_ccdc_logCPM","wtko_ccdc_PValue")
colnames(myfiles[[5]]) <- c("gene","wtko_dmw_logFC","wtko_dmw_logCPM","wtko_dmw_PValue")
colnames(myfiles[[6]]) <- c("gene","wtko_scan_logFC","wtko_scan_logCPM","wtko_scan_PValue")

library(plyr)
alldata<-join_all(myfiles, by = "gene", type = "full", match = "all")
library(ggplot2)
library(GGally)

# get rid of outliers
# https://www.r-bloggers.com/2020/01/how-to-remove-outliers-in-r/
boxplot(alldata$MF_ccdc_logFC, plot=FALSE)$out
# these outliers are the first quartile - 1.5 the interquartile range
# and the third quartile plus 1.5 the interquartile range
#MF_ccdc
outliers <- boxplot(alldata$MF_ccdc_logFC, plot=FALSE)$out
MF_ccdc_logFC_trim<-alldata[,c(1,2)]
MF_ccdc_logFC_trim<- MF_ccdc_logFC_trim[-which(MF_ccdc_logFC_trim$MF_ccdc_logFC %in% outliers),]
#MF_dmrt1L
outliers <- boxplot(alldata$MF_dmrt1L_logFC, plot=FALSE)$out
MF_dmrt1L_logFC_trim<-alldata[,c(1,5)]
MF_dmrt1L_logFC_trim<- MF_dmrt1L_logFC_trim[-which(MF_dmrt1L_logFC_trim$MF_dmrt1L_logFC %in% outliers),]
#MF_dmrt1S
outliers <- boxplot(alldata$MF_dmrt1S_logFC, plot=FALSE)$out
MF_dmrt1S_logFC_trim<-alldata[,c(1,8)]
MF_dmrt1S_logFC_trim<- MF_dmrt1S_logFC_trim[-which(MF_dmrt1S_logFC_trim$MF_dmrt1S_logFC %in% outliers),]
#wtko_dmw
outliers <- boxplot(alldata$wtko_dmw_logFC, plot=FALSE)$out
wtko_dmw_logFC_trim<- alldata[,c(1,14)]
wtko_dmw_logFC_trim<- wtko_dmw_logFC_trim[-which(wtko_dmw_logFC_trim$wtko_dmw_logFC %in% outliers),]
#wtko_scan
outliers <- boxplot(alldata$wtko_scan_logFC, plot=FALSE)$out
wtko_scan_logFC_trim<- alldata[,c(1,17)]
wtko_scan_logFC_trim<- wtko_scan_logFC_trim[-which(wtko_scan_logFC_trim$wtko_scan_logFC %in% outliers),]
#wtko_ccdc
outliers <- boxplot(alldata$wtko_ccdc_logFC, plot=FALSE)$out
wtko_ccdc_logFC_trim<- alldata[,c(1,11)]
wtko_ccdc_logFC_trim<- wtko_ccdc_logFC_trim[-which(wtko_ccdc_logFC_trim$wtko_ccdc_logFC %in% outliers),]


# combine no outlier files
# make a list of df
df_list <- list(MF_ccdc_logFC_trim, MF_dmrt1L_logFC_trim, MF_dmrt1S_logFC_trim, wtko_dmw_logFC_trim, wtko_scan_logFC_trim, wtko_ccdc_logFC_trim)
#merge all data frames in list
alldata_no_outliers <- df_list %>% reduce(full_join, by='gene')

colnames(alldata_no_outliers) <- c("gene","MF1","MF2",
                                   "MF3","dmw","scan",
                                   "ccdc")
library(ggplot2)

my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
  #  geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}

my_fn2 <- function(data, mapping, method="p", use="pairwise", ...){
  
  # grab data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # calculate correlation
  corr <- cor(x, y, method=method, use=use)
  
  # calculate colour based on correlation value
  # Here I have set a correlation of minus one to blue, 
  # zero to white, and one to red 
  # Change this to suit: possibly extend to add as an argument of `my_fn`
  colFn <- colorRampPalette(c("blue", "white", "red"), interpolate ='spline')
  fill <- colFn(100)[findInterval(corr, seq(-1, 1, length=100))]
  
  ggally_cor(data = data, mapping = mapping, ...) + 
    theme_void() +
    theme(panel.background = element_rect(fill=fill))
}

my_custom_smooth <- function(data, mapping, ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(color = I("blue")) + 
    geom_smooth(method = "lm", fill="blue", color="blue", ...)
  
  lmModel <- eval(substitute(lm(y ~ x, data = data), mapping))
  fs <- summary(lmModel)$fstatistic
  pValue <- pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
  
  if (pValue < 0.05) {
    p <- p + theme(
      panel.border = element_rect(
        color = "red", 
        size = 3,
        linetype = "solid",
        fill = "transparent"
      )
    )
  }
  
  p
}

p_ <- GGally::print_if_interactive
g<-ggpairs(alldata_no_outliers[,c(2:7)], 
        #upper = list(continuous = "density", combo = "box_no_facet"),
        #upper = list(continuous = wrap(ggally_cor, size = 2)), 
        upper = list(continuous = my_fn2),
        lower = list(continuous = my_fn)) +
        #lower = list(continuous = my_custom_smooth)) +
  #theme_bw() +
  theme(strip.background = element_rect(
      color="white", fill="white", size=1.5, linetype="solid")) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
box1_2 <- ggally_text("\nr = 0.159\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_3 <- ggally_text("\nr = -0.232\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_4 <- ggally_text("\nr = 0.278*\np = 0.116\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_5 <- ggally_text("\nr = 0.123\np = 0.283\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box1_6 <- ggally_text("\nr = 0.212\np = 0.972\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_3 <- ggally_text("\nr = -0.033\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_4 <- ggally_text("\nr = 0.111\np = 0.162\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_5 <- ggally_text("\nr = 0.138\np = 0.215\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box2_6 <- ggally_text("\nr = -0.222\np = 0.785\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_4 <- ggally_text("\nr = 0.331*\np = 0.028*\n",geom_text = ggplot2::aes(size = 6), color = I("red"))
box3_5 <- ggally_text("\nr = -0.139\np = 0.675\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box3_6 <- ggally_text("\nr = -0.021\np = 0.286\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box4_5 <- ggally_text("\nr = 0.143\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box4_6 <- ggally_text("\nr = 0.025\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
box5_6 <- ggally_text("\nr = -0.167\n\n",geom_text = ggplot2::aes(size = 6), color = I("black"))
g[1, 2] <- box1_2
g[1, 3] <- box1_3
g[1, 4] <- box1_4
g[1, 5] <- box1_5
g[1, 6] <- box1_6
g[2, 3] <- box2_3
g[2, 4] <- box2_4
g[2, 5] <- box2_5
g[2, 6] <- box2_6
g[3, 4] <- box3_4
g[3, 5] <- box3_5
g[3, 6] <- box3_6
g[4, 5] <- box4_5
g[4, 6] <- box4_6
g[5, 6] <- box5_6
# small function to display plots only if it's interactive

p_(g)
ggsave(file="STAR_edgeR_sexrelated_pairwise_unfiltered.pdf", g, width=10, height=4)
