#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
genofile <- as.character(args[1])
imagfile <- as.character(args[2])
covfile <- as.character(args[3])

out_prefix <- as.character(args[4])

min_genetic <- as.integer(args[5])
max_order <- as.integer(args[6])
outlier_criteria <- as.double(args[7])

source("/ifs/loni/faculty/kclark/candl/kristi/PheWAS/Neuroimaging_PheWAS.R")

Neuroimaging_PheWAS(genofile, imagfile, covfile, out_prefix, min_genetic, max_order, outlier_criteria)