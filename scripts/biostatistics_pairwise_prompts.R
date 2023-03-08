#!/usr/bin/env Rscript

##
# Working Directory
##

# set the working directory



##
# Packages
##

# install packages, if necessary


# import libraries



##
# Data
##

# import gene count data



##
# Pairwise Setup
##

# add grouping factor


# create DGE list object



##
# Pairwise Normalization
##

# plot the library sizes before normalization


# filter the list of gene counts based on expression levels


# view the number of filtered genes


# remove genes that are not expressed in either experimental condition


# calculate scaling factors


# compute counts per million (CPM) using normalized library sizes



##
# Plotting Palettes
##

# change the graphical parameters


# view all available ghibli palettes


# close the plot and return the display to the default graphical parameters


# retrieve the vector of colors associated with PonyoMedium


# view the selected color palette


# vector with a subset of colors associated with PonyoMedium



##
# Pairwise Data Exploration
##

# vector of shape numbers for the MDS plot


# vector of colors for the MDS plot


# add extra space to right of plot area and change clipping to figure


# MDS plot with distances approximating log2 fold changes


# place the legend outside the right side of the plot


# close the plot


# calculate the log CPM of the gene count data


# draw a heatmap of individual RNA-seq samples using moderated log CPM



##
# Pairwise Fitting
##

# estimate common dispersion and tagwise dispersions to produce a matrix of pseudo-counts


# plot dispersion estimates and biological coefficient of variation



##
# Pairwise Contrasts
##

### 
## treat_4h vs cntrl_4h
###

# perform an exact test for treat_4h vs cntrl_4h


# view the total number of differentially expressed genes at a p-value of 0.05


# plot log-fold change against log-counts per million with DE genes highlighted


# add blue lines to indicate 2-fold changes


# create a results table of DE genes


# add column for identifying direction of DE gene expression


# identify significantly up DE genes


# identify significantly down DE genes


# create volcano plot



###
## treat_24h vs cntrl_24h
###

# perform an exact test for treat_24h vs cntrl_24h


# view the total number of differentially expressed genes at a p-value of 0.05


# plot log-fold change against log-counts per million with DE genes highlighted


# add blue lines to indicate 2-fold changes


# create a results table of DE genes


# add column for identifying direction of DE gene expression


# identify significantly up DE genes


# identify significantly down DE genes


# create volcano plot


# identify significantly DE genes by FDR


# create filtered results table of DE genes



###
## treat_4h vs treat_24h
###

# perform an exact test for treat_4h vs treat_24h


# view the total number of differentially expressed genes at a p-value of 0.05


# plot log-fold change against log-counts per million with DE genes highlighted


# add blue lines to indicate 2-fold changes


# create a results table of DE genes


# add column for identifying direction of DE gene expression


# identify significantly up DE genes


# identify significantly down DE genes


# create volcano plot


# identify significantly DE genes by FDR


# create filtered results table of DE genes



###
## cntrl_4h vs cntrl_24h
###

# perform an exact test for cntrl_4h vs cntrl_24h


# view the total number of differentially expressed genes at a p-value of 0.05


# plot log-fold change against log-counts per million with DE genes highlighted


# add blue lines to indicate 2-fold changes


# create a results table of DE genes


# add column for identifying direction of DE gene expression


# identify significantly up DE genes


# identify significantly down DE genes


# create volcano plot


# identify significantly DE genes by FDR


# create filtered results table of DE genes



##
# Pairwise Results Exploration
##

# retrieve set of DE gene names for 24h contrast


# retrieve set of DE gene names for treat contrast


# retrieve set of DE gene names for cntrl contrast


# create combined list of DE gene names


# create venn diagram

