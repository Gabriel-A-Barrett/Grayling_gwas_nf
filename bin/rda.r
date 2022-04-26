#!/usr/bin/env Rscript
library(vcfR)
library(vegan)

rda <- function(vcf, pop, env)
{
    vcf <- read.vcfR(vcf)
    pop <- read.table(pop, header=FALSE, col.names = c("indv","populations"))
    genind <- vcfR2genind(vcf, pop=pop$populations)
    allele_frq = data.frame(rraf(genind, by_pop = TRUE, correction = FALSE), check.names = FALSE)
    
    #write.csv(allele_frq, file = "full_allele_frq.csv", row.names = TRUE)
    #allele_freqs_outlier = read.csv("full_allele_frq.csv", row.names = 1, check.names = FALSE)

    # Random Number Generator
    set.seed(123)

    # run rda on all environmental var. 
    rda <- rda(allele_frq ~ ., data = env, scale = FALSE)

    

}