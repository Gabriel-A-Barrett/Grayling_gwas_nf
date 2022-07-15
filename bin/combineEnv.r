#!/usr/bin/env Rscript

#args <- commandArgs(trailingOnly=TRUE)

name_candidates <- list.files(pattern = "\\_candidates.txt")

for (file in c(1:length(name_candidates))) {if(!(file.exists(name_candidates[file]) && (!file.info(name_candidates[file])$isdir))) stop(paste(name_candidates[file]," does not exist"))}

message("Input meta: ", paste(name_candidates, sep=" "))

read_candidates <- lapply(name_candidates, read.delim) # reads all dataframes in list and store in list read_candidates

library(dplyr)
options(dplyr.summarise.inform = FALSE)
library(ggplot2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write Summary DataFrame with number of candidates in each dataframe
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

summary <- data.frame(env=character(),candidates=numeric()) # initilize dataframe

for (file in c(1:length(read_candidates))) {
    env <- sub("_candidates.txt", '', basename(name_candidates[file]))

    summary[nrow(summary) + 1, 1] <- env
    summary[nrow(summary), 2] <- nrow(read_candidates[[file]])
}

write.table(x=summary,file="Summary_Candidates.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

binded_candidates <- dplyr::bind_rows(read_candidates)

binded_candidates %>% # bind rows to of list of dataframes
  group_by(id) %>%
  mutate(env_pvalue = paste0(env,"=",pvalue, collapse = ",")) %>%
  mutate(n=n()) %>% 
  arrange(desc(n)) %>%
  select("chrom:pos"=id,n,annotation,
         impact=annotation_impact,
         Gene=EggNOG.Predicted.Gene,
         description=EggNOG.Description,
         GO_terms=EggNOG.GO.Biological,
         env_pvalue) %>%
  write.table(x=.,file="FinalResults_Candidates.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

# Bar Chart Function

grp_barchart <- function(df, x = n, y = env, fill) {
  
  ggplot2::ggplot(data=df, aes(fill={{fill}},{{y}},{{x}})) + geom_bar(position="stack",stat="identity") + theme_bw() + ylab("counts") + xlab("environment")

  }

# Annotation Stacked Bar Plot VIS
png("Annotation_StackedBarChart.png")
par(mfrow=c(2,1))
for (type in c("annotation", "EggNOG.GO.Biological")) {

  if (type == "EggNOG.GO.Biological") {
    eggnog_biol <- binded_candidates %>% 
      mutate(GO = gsub("(*),.*", "\\1", .data[[type]])) %>% 
      group_by(env, GO) %>% 
      summarise(n=n()) %>% 
      grp_barchart(fill=GO) 
    
    } else {
  
    eggnog_annot <- binded_candidates %>% 
      group_by(env, .data[[type]]) %>% 
      summarise(n=n()) %>% 
      grp_barchart(fill=.data[[type]]) 
  
    }
  }

print(eggnog_biol)
print(eggnog_annot)

dev.off()