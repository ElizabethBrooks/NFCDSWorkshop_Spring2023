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
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")

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
# GLM Design
##

# import grouping factor
glm_targets <- read.csv(file="groupingFactors_tribolium.csv", row.names="sample")

# setup a design matrix
glm_group <- factor(paste(glm_targets$treatment, glm_targets$hours, sep="."))

# begin to construct the DGE list object
glm_list <- DGEList(counts=tribolium_counts, group=glm_group)

# add the sample names
colnames(glm_list) <- rownames(glm_targets)

# parametrize the experimental design with a one-way layout 
glm_design <- model.matrix(~ 0 + glm_group)

# add group names
colnames(glm_design) <- levels(glm_group)


##
# GLM Normalization
##

# filter the list of gene counts based on expression levels
glm_keep <- filterByExpr(glm_list)

# view the number of filtered genes
table(glm_keep)

# remove genes that are not expressed in either experimental condition
glm_list <- glm_list[glm_keep, , keep.lib.sizes=FALSE]

# calculate scaling factors
glm_list <- calcNormFactors(glm_list)

# compute counts per million (CPM) using normalized library sizes
norm_glm_list <- cpm(glm_list, normalized.lib.sizes=TRUE)


##
# GLM Fitting
##

# estimate common dispersion and tagwise dispersions to produce a matrix of pseudo-counts
glm_list <- estimateDisp(glm_list, glm_design, robust=TRUE)

# estimate the QL dispersions
glm_fit <- glmQLFit(glm_list, glm_design, robust=TRUE)

# plot the QL dispersions
plotQLDisp(glm_fit)


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

# view the selected color palette
ghibli_colors

# vector with a subset of colors associated with PonyoMedium
ghibli_subset <- c(ghibli_colors[3], ghibli_colors[6], ghibli_colors[4])


##
# GLM Contrasts
##

###
## treatment
###

# examine the overall effect of treatment
con.treatment <- makeContrasts(set.treatment = 
                               (treat.4h + treat.24h)/2
                               - (cntrl.4h + cntrl.24h)/2,
                               levels=glm_design)

# conduct gene wise statistical tests
anov.treatment <- glmTreat(glm_fit, contrast=con.treatment)

# view summary of DE genes
summary(decideTests(anov.treatment))

# create MD plot of DE genes
plotMD(anov.treatment)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# generate table of DE genes
tagsTbl_treatment <- topTags(anov.treatment, n=nrow(anov.treatment$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
tagsTbl_treatment$topDE <- "NA"

# identify significantly up DE genes
tagsTbl_treatment$topDE[tagsTbl_treatment$logFC > 1 & tagsTbl_treatment$FDR < 0.05] <- "UP"

# identify significantly down DE genes
tagsTbl_treatment$topDE[tagsTbl_treatment$logFC < -1 & tagsTbl_treatment$FDR < 0.05] <- "DOWN"

# create volcano plot
ggplot(data=tagsTbl_treatment, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
tagsTbl_treatment.glm_keep <- tagsTbl_treatment$FDR < 0.05

# create filtered results table of DE genes
tagsTbl_treatment_filtered <- tagsTbl_treatment[tagsTbl_treatment.glm_keep,]


###
## hours
###

# examine the overall effect of hours
con.hours <- makeContrasts(set.hours = 
                           (cntrl.24h + treat.24h)/2
                           - (cntrl.4h + treat.4h)/2,
                           levels=glm_design)

# conduct gene wise statistical tests
anov.hours <- glmTreat(glm_fit, contrast=con.hours)

# view summary of DE genes
summary(decideTests(anov.hours))

# create MD plot of DE genes
plotMD(anov.hours)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# generate table of DE genes
tagsTbl_hours <- topTags(anov.hours, n=nrow(anov.hours$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
tagsTbl_hours$topDE <- "NA"

# identify significantly up DE genes
tagsTbl_hours$topDE[tagsTbl_hours$logFC > 1 & tagsTbl_hours$FDR < 0.05] <- "UP"

# identify significantly down DE genes
tagsTbl_hours$topDE[tagsTbl_hours$logFC < -1 & tagsTbl_hours$FDR < 0.05] <- "DOWN"

# create volcano plot
ggplot(data=tagsTbl_hours, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
tagsTbl_hours.glm_keep <- tagsTbl_hours$FDR < 0.05

# create filtered results table of DE genes
tagsTbl_hours_filtered <- tagsTbl_hours[tagsTbl_hours.glm_keep,]


###
## interaction
###

# examine any interaction effect
con.interaction <- makeContrasts(set.interaction = 
                                 ((treat.4h + treat.24h)/2
                                 - (cntrl.4h + cntrl.24h)/2)
                                 - ((cntrl.24h + treat.24h)/2
                                 - (cntrl.4h + treat.4h)/2),
                                 levels=glm_design)

# conduct gene wise statistical tests
anov.interaction <- glmTreat(glm_fit, contrast=con.interaction)

# view summary of DE genes
summary(decideTests(anov.interaction))

# create MD plot of DE genes
plotMD(anov.interaction)

# add blue lines to indicate 2-fold changes
abline(h=c(-1, 1), col="blue")

# generate table of DE genes
tagsTbl_inter <- topTags(anov.interaction, n=nrow(anov.interaction$table), adjust.method="fdr")$table

# add column for identifying direction of DE gene expression
tagsTbl_inter$topDE <- "NA"

# identify significantly up DE genes
tagsTbl_inter$topDE[tagsTbl_inter$logFC > 1 & tagsTbl_inter$FDR < 0.05] <- "UP"

# identify significantly down DE genes
tagsTbl_inter$topDE[tagsTbl_inter$logFC < -1 & tagsTbl_inter$FDR < 0.05] <- "DOWN"

# create volcano plot
ggplot(data=tagsTbl_inter, aes(x=logFC, y=-log10(FDR), color = topDE)) + 
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(type = ghibli_subset, breaks = c("Up", "Down"))

# identify significantly DE genes by FDR
tagsTbl_inter.glm_keep <- tagsTbl_inter$FDR < 0.05

# create filtered results table of DE genes
tagsTbl_inter_filtered <- tagsTbl_inter[tagsTbl_inter.glm_keep,]


##
# GLM Results Exploration
##

# retrieve set of DE gene names for hours contrast
geneSet_hours <- rownames(tagsTbl_hours_filtered)

# retrieve set of DE gene names for interaction contrast
geneSet_interaction <- rownames(tagsTbl_inter_filtered)

# create combined glm_list of DE gene names
glm_list_venn <- list(hours = geneSet_hours, 
                          interaction = geneSet_interaction)

# create venn diagram
ggVennDiagram(glm_list_venn, label_alpha=0.25, category.names = c("hours","interaction")) +
  scale_color_brewer(palette = "Paired")
