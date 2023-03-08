#!/bin/bash
#$ -M netid@nd.edu
#$ -m abe
#$ -r n
#$ -N bioinformatics_driver_jobOutput

# load the bio module for the servers
module load bio

# run the shell script for data prep and alignment
bash bioinformatics_prep.sh

# run the R script for quantifying transcript data
Rscript bioinformatics_quantification.R
