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
# GLM Design
##

# import grouping factor


# setup a design matrix


# create DGE list object


# add the sample names


# parametrize the experimental design with a one-way layout 


# add group names



##
# GLM Normalization
##

# filter the list of gene counts based on expression levels


# view the number of filtered genes


# remove genes that are not expressed in either experimental condition


# calculate scaling factors


# compute counts per million (CPM) using normalized library sizes



##
# GLM Fitting
##

# estimate common dispersion and tagwise dispersions to produce a matrix of pseudo-counts


# estimate the QL dispersions


# plot the QL dispersions



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
# GLM Contrasts
##

###
## treatment
###

# examine the overall effect of treatment


# conduct gene wise statistical tests


# view summary of DE genes


# create MD plot of DE genes


# add blue lines to indicate 2-fold changes


# generate table of DE genes


# add column for identifying direction of DE gene expression


# identify significantly up DE genes


# identify significantly down DE genes


# create volcano plot


# identify significantly DE genes by FDR


# create filtered results table of DE genes



###
## hours
###

# examine the overall effect of hours


# conduct gene wise statistical tests


# view summary of DE genes


# create MD plot of DE genes


# add blue lines to indicate 2-fold changes


# generate table of DE genes


# add column for identifying direction of DE gene expression


# identify significantly up DE genes


# identify significantly down DE genes


# create volcano plot


# identify significantly DE genes by FDR


# create filtered results table of DE genes



###
## interaction
###

# examine any interaction effect


# conduct gene wise statistical tests


# view summary of DE genes


# create MD plot of DE genes


# add blue lines to indicate 2-fold changes


# generate table of DE genes


# add column for identifying direction of DE gene expression


# identify significantly up DE genes


# identify significantly down DE genes


# create volcano plot


# identify significantly DE genes by FDR


# create filtered results table of DE genes



##
# GLM Results Exploration
##

# retrieve set of DE gene names for hours contrast


# retrieve set of DE gene names for interaction contrast


# create combined glm_list of DE gene names


# create venn diagram

