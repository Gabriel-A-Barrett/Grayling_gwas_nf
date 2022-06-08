#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 5) {
    stop("Usage: baypass_median.r <env> <lea> <baypass> <snpeff vcf> <header key> provided incorrect number of arguments")
}

env <- args[1]
lea <- args[2]
baypass <- args[3]
vcf <- args[4]
header_key <- args[5]

if(!(file.exists(loci) && (!file.info(loci)$isdir))) stop("Second argument '<loci>' file does not exist")

message("Input meta (Arg 1:) ", env)
message("Input meta (Arg 2:) ", lea)
message("Input meta (Arg 3:) ", baypass)
message("Input meta (Arg 4:) ", header_key)

library(dplyr)
library(vcfR)

header_key <- read.table(header_key,sep="",header=FALSE,col.names=c("chrom","chrom_id")) 
    %>% mutate(chrom=gsub("[^0-9]+","",chrom))

lea <- read.table(lea,h=T,sep="")

baypass <- read.table(baypass,sep="\t",h=T) %>% select(!(n))

entap <- read.delim(entap, h=T, fill=T,sep="\t", na.strings=c("","NA"))  %>% 
    select(feature_id=Query.Sequence,EggNOG.Predicted.Gene,EggNOG.Description,EggNOG.GO.Biological)

vcfR_obj <- read.vcfR(vcf)
vcfR_df <- vcfR2tidy(vcfR_obj)
vcf_df_fix <- vcfR_df$fix %>% 
  select(chrom_id=CHROM,pos=POS,ANN) %>% 
  separate(col = ANN,into = c("allele","annotation","annotation_impact","gene_name",
                                   "gene_id","feature_type","feature_id","transcript_biotype",
                                   "rank","HGVS.c","HGVS.p","cDNA.pos/cDNA.length","cDNS.pos/cDNS.length",
                                   "AA.pos/AA.length","distance","err"), sep = "\\|",convert=T,remove=F,fill="right") %>%
  select(!(ANN)) %>% 
  cbind(.,BayPass,LEA) %>%
  # left join can add rows if there are multiple matches to the right table
  #left_join(.,entap, by = c("feature_id"="Query.Sequence")) %>%
  left_join(.,header_key,by="chrom_id") %>%
  mutate_all(na_if,"") #%>% select(-feature_id) #%>%  cbind(.,BayPass,LEA)

# Add row number within each group so that only the first record will be added
SNPs_df <- left_join(vcf_df_fix %>% group_by(feature_id) %>% dplyr::mutate(id1 = row_number()),
                        entap %>% group_by(feature_id) %>% dplyr::mutate(id2 = row_number()),
                        by = c("feature_id","id1"="id2")) %>% 
                        
                        # necessary for downstream processes?? forget where..... :(
                        ungroup()


png(paste0(env,"_univariateHist.png"))
par(mfrow=c(2,2))
for (stat in c("XtX","GD","BF","EA")) {

    column <- paste0(env,"_",stat)
    
    hist(SNPs_df[[column]], main = paste0(env,"_",stat))
}
dev.off()


