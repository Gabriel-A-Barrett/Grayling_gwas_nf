#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
    stop("Usage: meta_2_env.r <params.meta> <baypass groups> <baypass groups order>")
}

meta <- args[1]
grp_character <- args[2]
grp_order <- args[3]

if(!(file.exists(meta) && (!file.info(meta)$isdir))) stop("First argument '<meta.csv>' file does not exist")
#if(!(file.exists(indv) && (!file.info(indv)$isdir))) stop("Second argument '<env.indv>' file does not exist")
if(!(file.exists(grp_character) && (!file.info(grp_character)$isdir))) stop("Third argument '<env.grp>' file does not exists")
if(!(file.exists(grp_order) && (!file.info(grp_order)$isdir))) stop("Fourth argument '<env.frq.pop_order>' file does not exist")

message("Input meta (Arg 1:) ", meta)
message("Input meta (Arg 2:) ", grp_character)
message("Input meta (Arg 3:) ", grp_order)

library(dplyr)
#library(stringr)
library(tidyr)

# Baypass groupings for correctly structuring baypass env. input based on vcf2baypass.pl
grps <- read.table(grp_character, header = FALSE, col.names = c("id","pop"))
grp_order <- read.table(grp_order, header = FALSE, sep="",col.names="order")

env <- sub('.grp', '',basename(grp_character))

# extract environmental data
meta <- read.table(meta, header = T , sep = ',') %>% select(id=1,env)

# recode discrete variables to numeric for standardization
if (env == "H1_coast_noncoast"){
    meta <- mutate(meta, env = as.numeric(str_replace_all(meta[,2], c("noncoast" = "0", "coast" = "1"))))
} else if (env == "H2_Aufeis_above_below_coast"){
    meta <- mutate(meta, env = as.numeric(str_replace_all(meta[,2], c("above" = "-1", "below" = "0", "coast" = "1"))))
} else if (env == "H3_WShd_up_down_lower"){
    meta <- mutate(meta, env = as.numeric(str_replace_all(meta[,2], c("up" = "-1", "down" = "0", "lower" = "1"))))
} else {
    meta <- mutate(meta, env = meta[,2])
}

# Standardize all records
meta <- mutate(meta, standard_env = (as.numeric(env) - mean(as.numeric(env),na.rm=T)) / sd(as.numeric(env),na.rm=T))[,c(1,4)]
    
# Remove problematic records to match indv in vcf
correct <- left_join(grps, meta, by = "id")

#~~~~~~~~~~~~~~~~~
# Write Environmental Files
#~~~~~~~~~~~~~~~~~

# lfmm format
write.table(file = paste0(env,".env"), x = correct[,c(3)], row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# baypass env. format
correct %>% 
    group_by(pop) %>% 
    summarise(avg = mean(as.numeric(standard_env))) %>% 
    arrange(factor(pop,levels=grp_order$order)) %>% 
    tidyr::pivot_wider(.,names_from = "pop", values_from = avg) %>%
    write.table(x=.,file=paste0(env,".bayenv"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
