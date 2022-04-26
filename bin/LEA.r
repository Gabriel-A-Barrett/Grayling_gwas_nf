#!/usr/bin/env Rscript

#library("LEA", lib=rlibs)
#library("dplyr", lib=rlibs)

lea <- function(env, ped, env_file)
{ if (is.character(env))
    if (env != "pit"){
        latentFactors = 5
    } else {
        latentFactors = 3
    }
    geno <- ped2geno(ped, output.file=paste0(env,".geno"), force=TRUE)
    lfmm <- ped2lfmm(ped, output.file=paste0(env,".lfmm"), force=TRUE)

    ################
    # Differentiation
    proj.snmf <- snmf(geno,K=latentFactors,entropy=T,ploidy=2,project="new",alpha=10,tolerance=0.0001,repetitions=15,iterations=1000000,CPU=25,percentage=.75)
    # fst values
    
    best <- which.min(cross.entropy(proj.snmf, K = latentFactors))
    fst.values <- fst(proj.snmf, K = latentFactors, run = best)
    # z-scores
    n <- dim(Q(proj.snmf, K = latentFactors, run = best))[1]
    fst.values[fst.values<0] <- 0.000001
    z.scores <- sqrt(fst.values*(n - 5)/(1 - fst.values))
    write.table(x=z.scores,file=paste(env,"_GD_zscores.txt"),quote=F,row.names=F,col.names=paste0(env,"_z"))
    ################


    ################
    # Association
    proj.lfmm <- lfmm(lfmm, env_file, K = latentFactors, repetitions = 15, project = "new", iterations = 100000, burnin = 50000, CPU = 25, missing.data = TRUE)
    # z-scores from all repititions
    zv <- data.frame(z.scores(proj.lfmm, K = latentFactors))
    zv_median <- zv %>% rowwise() %>% mutate(paste0(env,"_z"), median(c_across(everything()))) %>% select(paste0(env,"_z"))
    #write.table(x = zv, file = "abc_EA_zscores.txt",quote=F,row.names=F,col.names=c("lat_z1","lat_z2","lat_z3","lat_z4","lat_z5","lat_z6","lat_z7","lat_z8","lat_z9","lat_z10","lat_z11","lat_z12","lat_z13","lat_z14","lat_z15"))
    write.table(x=zv_median, file = paste(env,"_EA_zscores.txt"),quote=F,row.names=F,col.names=T)
    #################
}