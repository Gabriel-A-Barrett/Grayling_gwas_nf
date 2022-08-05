#!/usr/bin/env Rscript
Sys.setlocale(category = "LC_ALL", locale = "pt_PT.CP1252")

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 2) {
    stop("Usage: populations.r <popSumStats> <fstSumStats> incorrect number of arguments")
}

popSumStats <- args[1]
fstSumStats <- args[2]

if(!(file.exists(popSumStats) && (!file.info(popSumStats)$isdir))) stop("First argument '<popSumStats>'' file does not exist")
if(!(file.exists(fstSumStats) && (!file.info(fstSumStats)$isdir))) stop("First argument '<fstSumStats>'' file does not exist")

message("Input meta (Arg 1:) ", popSumStats)
message("Input meta (Arg 2:) ", fstSumStats)

library(dplyr)
library(ggplot2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write Markdown Dataframe Within-Population
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

popSumStats <- read.table(popSumStats,sep="\t",header=F, 
                          skip=grep("# All positions (variant and fixed)",readLines(popSumStats),fixed=TRUE),
                          col.names=c("pop","private_sites","variant_sites","fixed_sites", "polymorphic_sites","percent_polymorphic_sites","indvs","indvs_var","indvs_stderr","p","p_var","p_stderr","obs_het","obs_het_var","obs_het_stderr","obs_hom","obs_hom_var","obs_hom_stderr","exp_het","exp_het_var","exp_het_stderr","exp_hom","exp_hom_var","exp_hom_stderr","pi","pi_var","pi_stderr","fis","fis_var","fis_stderr")) # grep returns number

phi_stats <- popSumStats %>%
  select(pop,obs_het,obs_het_var,exp_het,exp_het_var,pi,pi_var,fis,fis_var)

# Round all rows in columns 3 to 6 to 3 decimal places
phi_stats[,3:8] <- lapply(phi_stats[,3:8], round, 3)

# Combine Variation with respective column "\u00B1" = +/- 
phi_stats$Obs_Het <- paste(phi_stats$obs_het, "\u00B1", phi_stats$obs_het_var)
phi_stats$Exp_Het <- paste(phi_stats$exp_het, "\u00B1", phi_stats$exp_het_var)
phi_stats$Pi <- paste(phi_stats$pi, "\u00B1", phi_stats$pi_var)
phi_stats$Fis <- paste(phi_stats$fis, "\u00B1", phi_stats$fis_var)

# Rename columns, asterisks convert to subscripts in markdown language
phi_dataframe <- phi_stats %>%
  select("Population"=pop,"H~o~"=Obs_Het,"H~e~"=Exp_Het,"pi"=Pi,"F~is~"=Fis)

# Convert dataframe into markdown flavor
knitr::kable(phi_dataframe) %>% writeLines(.,"withinPopSumStats.md")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# point plot w/ variation Within-Population
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Turn into user input
factors <- c("IMup.a", "IMup", "IMlow.a","IMlow", "IM3.a", "LL3.a", "Itk3", "Itk4", "Itk4.5", "Kup.a", "Kup2", "Kup5", "Kup6", "Kup7", "KLAS5", "Kup8", "Lkup.a", "LKup0", "LKup1", "LKup3", "OksHW.a", "Oks0", "OksZev", "Oks2", "Oks2.5", "Oks3","OksSag.a", "LSagHV", "LSag2", "LSagTC")

withinPopGraph <- function(data=df,y=pop,x=.data[[stat]],var=.data[[paste0(stat,"_var")]],xlab="Population",limits=FALSE,ylab=stat){

  ggplot(data=df,aes(x={{x}},y={{y}})) +
    geom_point(stat="identity") + 
    geom_linerange(aes(xmin={{x}} - {{var}}, xmax={{x}} + {{var}})) + # wrap all variables with {{}} to prevent data masking issues. see https://rlang.r-lib.org/reference/topic-data-mask-ambiguity.html
    labs(x=xlab,y=ylab) + 
    theme_bw() 
  
  }

for (stat in c('fis', 'obs_het', 'pi')) {
    print(stat)
    df <- phi_stats %>% select(stat, paste0(stat,"_var"),pop)
    
    plot <- withinPopGraph()
    
    png(paste(stat,"_PointPlot.png"))    
    print(plot)
    dev.off()
  
}

# ~~~~~~~~~~~~~~~~~~~
# Fst HeatMap
# ~~~~~~~~~~~~~~~~~~~

fstSumStats <- read.table()