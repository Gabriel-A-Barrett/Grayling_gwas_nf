#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
    stop("Usage: write_grps.r <meta> <vcfIndvList> <first_env_column>")
}

meta <- args[1]
vcfIndvList <- args[2]
first_env_column <- as.integer(args[3])

if(!(file.exists(meta) && (!file.info(meta)$isdir))) stop("First argument '<meta.csv>' file does not exist")
if(!(file.exists(vcfIndvList) && (!file.info(vcfIndvList)$isdir))) stop("Second argument '<prefix.vcf.indv>' does not exist")

# goes to stdout (most likely)
message("Input meta          (Arg 1): ", meta)
message("VCF individual list (Arg 2): ", vcfIndvList)

library("dplyr")

meta <- read.table(file = meta,sep = ",", header = T) %>% select(id=1,pop=2,first_env_column:last_col()) # read in meta starting at the 9th column to the last in the dataset
vcfindv <- read.table(file = vcfIndvList, sep = '', header = F, col.names = 'id') # Correct for the number of individuals based on bcftools query

# write env proceeding 2nd col.  
for (i in 3: ncol(meta)) { 
  tmp <- meta[,c(1,2,i)] # subset to 3 col. # ID. Pop. Env.
  
  title <- toString(names(tmp)[3])
  filename <- paste("./",trimws(title),".grp") %>% gsub("[[:space:]]","",.) # Group file removing any white spaces
  
  # Remove problematic columns
  clean <- filter(tmp, !(tmp[,3] == "Na" | tmp[,3] == "" | tmp[,3] == "GCL"))

  l=unique(c(as.character(clean$pop)))
  
  grp <- data.frame(id=clean$id,grp=as.numeric(factor(clean$pop, levels=l)))
  correct <- merge(grp, vcfindv,by = "id",all.x = FALSE, all.y = FALSE)
  
  write.table(file = filename,x = correct, quote = F,sep = "\t",col.names = F,row.names = F)
  
  i <- i + 1 # move to the next column 
}

citation("dplyr")
sessionInfo()
