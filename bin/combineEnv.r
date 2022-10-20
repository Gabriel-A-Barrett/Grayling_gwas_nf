#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
    stop("Usage: combineEnv.r <meta> <genotypes> incorrect number of arguments")
}

meta <- args[1]
genotypes <- args[2]

name_candidates <- list.files(pattern = "\\_candidates.txt")

for (file in c(1:length(name_candidates))) {if(!(file.exists(name_candidates[file]) && (!file.info(name_candidates[file])$isdir))) stop(paste(name_candidates[file]," does not exist"))}
if(!(file.exists(meta) && (!file.info(meta)$isdir))) stop("First argument '<meta>'' file does not exist")
if(!(file.exists(genotypes) && (!file.info(genotypes)$isdir))) stop("First argument '<genotypes>'' file does not exist")

message("Input meta: ", paste(name_candidates, sep=" "))
message("Input meta: ", meta)
message("Input meta: ", genotypes)

library(plyr)
library(dplyr)
options(dplyr.summarise.inform = FALSE)
library(ggplot2)
library(tidyr)

read_candidates <- lapply(name_candidates, read.delim) # reads all dataframes in list and store in list read_candidates

pop <- read.table(meta,h=T,sep=",") %>% dplyr::select(ind = 1, pop = 2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write Summary DataFrame with number of candidates in each dataframe
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

summary <- data.frame(env=character(),candidates=numeric()) # initilize dataframe

for (file in c(1:length(read_candidates))) {
  env <- sub("_candidates.txt", '', basename(name_candidates[file]))
  
  summary[nrow(summary) + 1, 1] <- env
  summary[nrow(summary), 2] <- nrow(read_candidates[[file]])
}

write.table(x=summary %>% arrange(desc(candidates)),file="Summary_Candidates.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compile Candidates into one Dataframe
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

binded_candidates <- dplyr::bind_rows(read_candidates)

grouped_candidates <- binded_candidates %>% 
  group_by(id, .drop=FALSE) %>%
  dplyr::mutate(id = gsub(':',"_",id),
        env_pvalue = paste0(env,"=",pvalue, collapse = ","),
        n=n()) %>%
  arrange(desc(n)) %>%
  distinct(id, .keep_all=TRUE) %>%
  select("chrom_pos"=id, 
        n, annotation,
        impact=annotation_impact,
        gene=EggNOG.Predicted.Gene,
        description=EggNOG.Description,
        GO_terms=EggNOG.GO.Biological,
        env_pvalue, chrom)

write.table(x=grouped_candidates,file="FinalResults_Candidates.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VIS Pie Genotypes & impact level bar chart
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

filtered_candidates <- grouped_candidates %>% filter_at(vars(gene, description), any_vars(!is.na(.)))

# Match strings to color
colors <- c("coral", "chocolate", "burlywood", "bisque","black")
names <- c("Ref_Homo","Hetero_1","Hetero_2","Alt_Homo", "Missing")
names(colors) <- names

# VIS. Levels
adjusted_pop_order_list <- c("IMup.a", "IMup", "IMlow.a","IMlow", "IM3.a", "LL3.a", "Itk3", "Itk4", "Itk4.5", "Kup.a", "Kup2", "Kup5", "Kup6", "Kup7", "KLAS5", "Kup8", "Lkup.a", "LKup0", "LKup1", "LKup3", "OksHW.a", "Oks0", "OksZev", "Oks2", "Oks2.5", "Oks3","OksSag.a", "LSagHV", "LSag2", "LSagTC")
location_pop_order_list <- c("IMO", "IMR1", "IMR2", "BJANE", "LIMS", "IMR3", "IM3", "LL3", "Itk3", "Itk4", "Itk4.5", "GCL", "CGK", "Kup2", "KUS", "Toolik-S3", "LTER", "KupR2", "Toolik-TOAS", "Kup6", "KupR4", "Kup7", "KLAS5", "Kup", "Kup8", "LKup", "LKup0", "LKup1", "LKup3", "OksLCS", "UZ", "OKm1", "OksR1", "Oks0","OksLTER", "OksBH", "Oks2", "OksRN4", "Oks3", "OKS3", "OksSag", "CGO3", "LSag", "HV", "LSag2", "TC")

# Dataset with candidate SNP genotypes as characters for visualization
geno_long <- read.table(genotypes,h=T) %>% 
  mutate(chrom_pos = paste(CHROM,"_",POS, sep="")) %>% 
  right_join(.,grouped_candidates, by = "chrom_pos") %>% 
  #filter(id %in% filtered_candidates$chrom_pos) %>% 
  pivot_longer(cols = starts_with("Golden"), names_to = "ind", values_to = "geno") %>%
  left_join(.,pop,by="ind") %>%
  mutate(genotype = as.character(recode(geno, './.' = 'Missing', 
                                    '0/0' = 'Ref_Homo', 
                                    '1/0' = 'Hetero_1', 
                                    '0/1' = 'Hetero_2', 
                                    '1/1' = 'Alt_Homo',
                                    '.|.' = 'Missing', 
                                    '0|0' = 'Ref_Homo', 
                                    '1|0' = 'Hetero_1', 
                                    '0|1' = 'Hetero_2', 
                                    '1|1' = 'Alt_Homo',
                                    '2/1' = 'Hetero_1',
                                    '1/2' = 'Hetero_2',
                                    '2|1' = 'Hetero_1',
                                    '1|2' = 'Hetero_2',
                                    '2/2' = 'Alt_Homo',
                                    '2|2' = 'Alt_Homo',
                                    .default = "Missing"))) %>%
  mutate(level = factor(genotype, levels = rev(as.character(names)))) %>%
  arrange(pop,level) %>% 
  group_by(chrom_pos, pop, genotype) %>% 
  dplyr::mutate(n = dplyr::n()) %>%
  distinct(chrom_pos,.keep_all = TRUE)


# Bar Chart
grp_barchart <- function(df, x = n, y = env, fill, xlab, ylab, title="", sub="", caption="", legend_title="",genotypes = FALSE) {
  
  barchart <- ggplot2::ggplot(data=df, aes({{y}},{{x}})) + geom_bar(aes(fill={{fill}}),position="stack",stat="identity") + 
    labs(title=title, subtitle=sub, y=ylab, x=xlab, caption=caption, fill=legend_title) +  
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme_bw()
  
  # Condition = apply color pallete 
  if(genotypes) {
    barchart + scale_fill_manual(name = names, values = colors)} 
  else {
    barchart}
}

# Pie Chart 
grp_piechart <- function(df, y = n, fill = genotype) {

  dlply(df, .(pop), function(z)
    ggplot(z, aes(x="", y={{y}}, fill=as.factor({{fill}}))) +
      geom_bar(stat='identity', width=1, color = "black") +
      coord_polar(theta='y') +
      scale_fill_manual(name = names, values = colors) +
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank()))
}

# For every candidate SNP draw bar chart and pie chart of pop genotype proportions
for (candidate in unique(geno_long %>% distinct(chrom_pos) %>% pull(chrom_pos))) {
  
  if (file.exists(paste("./",candidate,"/")) == FALSE) {
  dir.create(paste0("./plots/",candidate,"/",collapse = ""),showWarnings = FALSE, recursive=TRUE)} # create candidate dir. for storing plots, collapse hidden space in front
  
  geno_long_candidate <- geno_long %>% filter(chrom_pos == candidate)
  
  png(file=paste0("./plots/",candidate,"/",candidate,"_geno_barchart.png",collapse = ""))
  barchart <- grp_barchart(df=geno_long_candidate, 
               y=n, x=factor(pop, level = location_pop_order_list), fill=genotype, # variables 
               ylab="population", xlab="count", legend_title="genotype", # Titles
               title = geno_long_candidate$Gene, sub = geno_long_candidate$description,
               genotypes=TRUE)
  
  print(barchart)
  dev.off()
  
  pies <- grp_piechart(df=geno_long %>% filter(chrom_pos == candidate))
  
  for (pop in c(1:length(pies))){
    png(file=paste0("./plots/",candidate,"/",names(pies[pop]),"_geno_piechart.png"))
    print(pies[[pop]])
    dev.off()
    }
}


# Annotation Stacked Bar Plot VIS
for (type in c("annotation", "EggNOG.GO.Biological")) {
  
  if (type == "EggNOG.GO.Biological") {
    eggnog_biol <- binded_candidates %>% 
      dplyr::mutate(GO = gsub("(*),.*", "\\1", .data[[type]])) %>%
      filter(!is.na(GO)) %>% # remove rows with no Gene Ontology term
      group_by(env, GO) %>% 
      dplyr::summarise(n=n()) %>% 
      grp_barchart(fill=GO, xlab="environment", ylab="counts") 
    
  } else {
    
    eggnog_annot <- binded_candidates %>% 
      group_by(env, .data[[type]]) %>% 
      dplyr::summarise(n=n()) %>% 
      grp_barchart(fill=.data[[type]], xlab="environment",ylab="counts") 
    
  }
}

png("EggNOG_StackedBarChart.png")
print(eggnog_biol)
dev.off()

png("Annotation_StackedBarChart.png")
print(eggnog_annot)
dev.off()
