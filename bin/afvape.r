#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
    stop("Usage: baypass_median.r <env> <loci file> <max iter.> incorrect number of arguments")
}

vcf <- args[1]
popmap <- args[2]
total_permutations <- args[3]

if(!(file.exists(vcf) && (!file.info(vcf)$isdir))) stop("Second argument '<vcf>' file does not exist")
if(!(file.exists(popmap) && (!file.info(popmap)$isdir))) stop("Second argument '<popmap>' file does not exist")

message("Input meta (Arg 1:) ", vcf)
message("Input meta (Arg 2:) ", popmap)
message("Input meta (Arg 3:) ", )

library(vcfR)

vcf_in <- read.vcfR(vcf,verbose=F)
popmap <- read.table(popmap)

#build_vector_list <- 

# splitting analysis over chromosomes

