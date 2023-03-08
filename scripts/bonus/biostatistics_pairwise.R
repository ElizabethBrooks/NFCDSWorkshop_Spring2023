#!/usr/bin/env Rscript

##
# Working Directory
##

# set the working directory
#setwd("/YOUR/PATH/")
setwd("/Users/bamflappy/Desktop/NFCDSWorkshop_Fall2022-main/data")


##
# Packages
##

# install packages, if necessary
#install.packages("ggplot2")
#install.packages("ghibli")
#install.packages("ggVennDiagram")
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")

# import libraries
library(ggplot2)
library(ghibli)
library(ggVennDiagram)
library(edgeR)


##
# Data
##

# import gene count data
tribolium_counts <- read.csv("TriboliumCounts.csv", row.names="X")


##
# Pairwise Setup
##

# add grouping factor
group <- factor(c(rep("cntrl_4h",3), rep("treat_4h",3), rep("cntrl_24h",3), rep("treat_24h",3)))
# why can't we just do the following?
#group <- factor(colnames(tribolium_counts))

# begin to construct the DGE list object
list <- DGEList(counts=tribolium_counts,group=group)


##
# Pairwise Normalization
##

# plot the library sizes before normalization
barplot(list$samples$lib.size*1e-6, names=1:12, ylab="Library size (millions)")

# filter the list of gene counts based on expression levels
keep <- filterByExpr(list)

# view the number of filtered genes
table(keep)

# remove genes that are not expressed in either experimental condition
list <- list[keep, , keep.lib.sizes=FALSE]

# calculate scaling factors
list <- calcNormFactors(list)

# compute counts per million (CPM) using normalized library sizes
normList <- cpm(list, normalized.lib.sizes=TRUE)


##
# Plotting Palettes
##

# change the graphical parameters
par(mfrow=c(9,3))

# view all available ghibli palettes
for(i in names(ghibli_palettes)) print(ghibli_palette(i))

# close the plot and return the display to the default graphical parameters
dev.off()

# retrieve the vector of colors associated with PonyoMedium
ghibli_colors <- ghibli_palette("PonyoMedium", type = "discrete")

# vector with a subset of colors associated with PonyoMedium
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4])


##
# Pairwise Data Exploration
##

# vector of shape numbers for the MDS plot
points <- c(0,1,15,16)

# vector of colors for the MDS plot
colors <- rep(c(ghibli_colors[3], ghibli_colors[6]), 2)

# add extra space to right of plot area and change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

# MDS plot with distances approximating log2 fold changes
plotMDS(list, col=colors[group], pch=points[group])

# place the legend outside the right side of the plot
legend("topright", inset=c(-0.4,0), legend=levels(group), pch=points, col=colors)

# close the plot
dev.off()

# calculate the log CPM of the gene count data
logcpm <- cpm(list, log=TRUE)

# draw a heatmap of individual RNA-seq samples using moderated log CPM
heatmap(logcpm)


##
# Pairwise Fitting
##

# estimate common dispersion and tagwise dispersions to produce a matrix of pseudo-counts
list <- estimateDisp(list)

# plot dispersion estimates and biological coefficient of variation
plotBCV(list)


##
# Pairwise Contrasts
##

### 
## treat_4h vs cntrl_4h
###

# perform an exact test for treat_4h vs cntrl_4h
tested_4h <- exactTest(list, pair=c("cntrl_4h", "treat_4h"))

# view the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_4h))

# plot log-fold change against log-counts per million with DE genes highlighted
plotMD(tested_4h)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# create a results table of DE genes
resultsTbl_4h <- topTags(tested_4h, n=nrow(tested_4h$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
resultsTbl_4h$topDE <- "NA"

# identify significantly up DE genes
resultsTbl_4h$topDE[resultsTbl_4h$logFC > 1 & resultsTbl_4h$FDR < 0.05] <- "Up"

# identify significantly down DE genes
resultsTbl_4h$topDE[resultsTbl_4h$logFC < -1 & resultsTbl_4h$FDR < 0.05] <- "Down"

# create volcano plot
ggplot(data=resultsTbl_4h, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))


###
## treat_24h vs cntrl_24h
###

# perform an exact test for treat_24h vs cntrl_24h
tested_24h <- exactTest(list, pair=c("cntrl_24h", "treat_24h"))

# view the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_24h))

# plot log-fold change against log-counts per million with DE genes highlighted
plotMD(tested_24h)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# create a results table of DE genes
resultsTbl_24h <- topTags(tested_24h, n=nrow(tested_24h$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
resultsTbl_24h$topDE <- "NA"

# identify significantly up DE genes
resultsTbl_24h$topDE[resultsTbl_24h$logFC > 1 & resultsTbl_24h$FDR < 0.05] <- "Up"

# identify significantly down DE genes
resultsTbl_24h$topDE[resultsTbl_24h$logFC < -1 & resultsTbl_24h$FDR < 0.05] <- "Down"

# create volcano plot
ggplot(data=resultsTbl_24h, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
resultsTbl_24h.keep <- resultsTbl_24h$FDR < 0.05

# create filtered results table of DE genes
resultsTbl_24h_filtered <- resultsTbl_24h[resultsTbl_24h.keep,]


###
## treat_4h vs treat_24h
###

# perform an exact test for treat_4h vs treat_24h
tested_treat <- exactTest(list, pair=c("treat_24h", "treat_4h"))

# view the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_treat))

# plot log-fold change against log-counts per million with DE genes highlighted
plotMD(tested_treat)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# create a results table of DE genes
resultsTbl_treat <- topTags(tested_treat, n=nrow(tested_treat$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
resultsTbl_treat$topDE <- "NA"

# identify significantly up DE genes
resultsTbl_treat$topDE[resultsTbl_treat$logFC > 1 & resultsTbl_treat$FDR < 0.05] <- "Up"

# identify significantly down DE genes
resultsTbl_treat$topDE[resultsTbl_treat$logFC < -1 & resultsTbl_treat$FDR < 0.05] <- "Down"

# create volcano plot
ggplot(data=resultsTbl_treat, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
resultsTbl_treat.keep <- resultsTbl_treat$FDR < 0.05

# create filtered results table of DE genes
resultsTbl_treat_filtered <- resultsTbl_treat[resultsTbl_treat.keep,]


###
## cntrl_4h vs cntrl_24h
###

# perform an exact test for cntrl_4h vs cntrl_24h
tested_cntrl <- exactTest(list, pair=c("cntrl_24h", "cntrl_4h"))

# view the total number of differentially expressed genes at a p-value of 0.05
summary(decideTests(tested_cntrl))

# plot log-fold change against log-counts per million with DE genes highlighted
plotMD(tested_cntrl)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# create a results table of DE genes
resultsTbl_ncntrl <- topTags(tested_cntrl, n=nrow(tested_cntrl$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
resultsTbl_ncntrl$topDE <- "NA"

# identify significantly up DE genes
resultsTbl_ncntrl$topDE[resultsTbl_ncntrl$logFC > 1 & resultsTbl_ncntrl$FDR < 0.05] <- "Up"

# identify significantly down DE genes
resultsTbl_ncntrl$topDE[resultsTbl_ncntrl$logFC < -1 & resultsTbl_ncntrl$FDR < 0.05] <- "Down"

# create volcano plot
ggplot(data=resultsTbl_ncntrl, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
resultsTbl_cntrl.keep <- resultsTbl_ncntrl$FDR < 0.05

# create filtered results table of DE genes
resultsTbl_cntrl_filtered <- resultsTbl_ncntrl[resultsTbl_cntrl.keep,]


##
# Pairwise Results Exploration
##

# retrieve set of DE gene names for 24h contrast
geneSet_24h <- rownames(resultsTbl_24h_filtered)

# retrieve set of DE gene names for treat contrast
geneSet_treat <- rownames(resultsTbl_treat_filtered)

# retrieve set of DE gene names for cntrl contrast
geneSet_cntrl <- rownames(resultsTbl_cntrl_filtered)

# create combined list of DE gene names
list_venn <- list(h24 = geneSet_24h, 
                  treat = geneSet_treat, 
                  cntrl = geneSet_cntrl)

# create venn diagram
ggVennDiagram(list_venn, label_alpha=0.25, category.names = c("24h","treat","cntrl")) +
  scale_color_brewer(palette = "Paired")
