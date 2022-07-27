#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 6) {
    stop("Usage: baypass_median.r <env> <lea> <baypass> <snpeff vcf> <entap> <header key> provided incorrect number of arguments")
}

env <- args[1]
baypass <- args[2]
lea <- args[3]
vcf <- args[4]
entap <- args[5]
header_key <- args[6]

if(!(file.exists(lea) && (!file.info(lea)$isdir))) stop("Second argument '<lea>' file does not exist")
if(!(file.exists(baypass) && (!file.info(baypass)$isdir))) stop("Third argument '<baypass>' file does not exist")
if(!(file.exists(vcf) && (!file.info(vcf)$isdir))) stop("Fourth argument '<vcf>' file does not exist")
if(!(file.exists(entap) && (!file.info(entap)$isdir))) stop("Fifth argument '<entap>' file does not exist")
if(!(file.exists(header_key) && (!file.info(header_key)$isdir))) stop("Sixth argument '<header_key>' file does not exist")

message("Input meta (Arg 1:) ", env)
message("Input meta (Arg 2:) ", baypass)
message("Input meta (Arg 3:) ", lea)
message("Input meta (Arg 4:) ", vcf)
message("Input meta (Arg 5:) ", entap)
message("Input meta (Arg 6:) ", header_key)

library(dplyr)
library(vcfR)
library(tidyr)
library(ggplot2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# COMPILE INFORMATION INTO SNPs_df
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

header_key <- read.table(header_key,sep="",header=FALSE,col.names=c("chrom","chrom_id")) %>% 
    mutate(chrom=gsub("[^0-9]+","",chrom))

lea <- read.table(lea,h=T,sep="")

baypass <- read.table(baypass,sep="\t",h=T) #%>% select(!(ends_with("Pearson")))

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
  cbind(.,baypass,lea) %>% # breaks when argument contains different number of rows, could switch to left_join based on Chromosome_Position
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

# write.table(file=paste0(env,"_SNPs_df.txt"),x=SNPs_df,row.names=FALSE,col.names=TRUE,quote=FALSE)

png(paste0(env,"_univariateHist.png"))
par(mfrow=c(2,3))
for (stat in c("XtX","GD","BF","EA", "Pearson")) {

    column <- paste0(env,"_",stat)
    
    hist(SNPs_df[[column]], main = paste0(env,"_",stat))
}
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate Mahalanobis Distance 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SNPs_df$mahalanobis <- mahalanobis(select(SNPs_df,starts_with(paste0(env))), colMeans(select(SNPs_df,starts_with(paste0(env)))), cov(select(SNPs_df,starts_with(paste0(env)))))

SNPs_df$pvalue <- pchisq(SNPs_df$mahalanobis,df=4,lower.tail=FALSE)

# Mahalanobis Distance Histogram
png(paste0(env,"_MahaDist_pvalue_hist.png"))
hist(SNPs_df$pvalue, main = paste0(env," Maha pvalue Hist"))
dev.off()

# Outlier List
SNPs_df %>% 
    #group_by(gene_name = SNPs_df$gene_name) %>%
    #dplyr::summarise(gene_count=n()) %>% 
    #ungroup() %>%
    mutate(env = paste0(env), id = paste0(chrom_id,"_",pos)) %>%
    select(id, annotation, gene_name,
            annotation_impact,HGVS.p,
            EggNOG.Predicted.Gene,EggNOG.Description,
            EggNOG.GO.Biological,chrom,env,pvalue) %>% 
    filter(pvalue < .00001) %>% 
    arrange(pvalue) %>%
    write.table(file=paste0(env,"_candidates.txt"),x=.,quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Manhattan Plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#SNPs_df <- SNPs_df %>% left_join(.,chrom_key,by="chrom_id")

don <- SNPs_df %>% 
    select(chrom,pos,pvalue) %>% 
    mutate(chrom = factor(chrom, levels = c(1:25))) %>% 
    filter(!(chrom == "NA")) %>%
    group_by(chrom) %>%
    dplyr::summarise(chr_len=max(pos)) %>%
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    dplyr::left_join(SNPs_df, ., by=c("chrom"="chrom")) %>%
    mutate(BPcum=pos+tot) %>%
    arrange(chrom,pos)

axis_df <- don %>% 
    group_by(chrom) %>%    
    dplyr::summarize(center=(max(BPcum) + min(BPcum)) / 2)

png(paste0(env,"_manhattanplot.png"), width=1285, height=240)
manhattan <- ggplot(don, aes(x=BPcum, y=-log10(pvalue))) +
    # Show all points
    geom_point(aes(color=as.factor(chrom)), alpha=0.85, size=2.3) +
    # Rotate between grey and lightgrey colors
    scale_color_manual(values = rep(c("grey", "lightgrey"), 22 )) +
    #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2.7) + 
    
    # custom X axis
    scale_x_continuous(label = axis_df$chrom, breaks= axis_df$center, expand=c(0,0)) + xlab("") + ylab("") +
    # Custom Y axis
    scale_y_continuous(expand = c(0, 0), limits = c(min(-log10(don$pvalue)) -.25, max(-log10(don$pvalue)) + .5)) + ylab("") + ylab(paste0(env)) + #labs(subtitle=paste0("n = ",NumberOfIndividuals)) +  # remove space between plot area and x axis with expand ORIGINAL: max(don$harmony) + 1
  
    # Custom the theme:
    theme_bw() +
    theme( legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = NULL,
      axis.text = element_text(size = 20.5),
      axis.title = element_text(size = 20),
      axis.title.y=element_text(angle=90,hjust=1),
      axis.text.y=element_text(size=12)) + 
      geom_hline(yintercept = -log10(.000001), linetype="dashed", size = .35, colour = "red")

print(manhattan)

dev.off()