# Arctic Grayling Landscape Genomics

nextflow pipeline for reproducing analysis done under Arctic Grayling

## Background

In an effort to support reproducibility of bioinformatic workflows and simply make re-running code easier in the event a more samples come in, 
I strung together the Arctic grayling landscape genomics analysis into NextFlow. This NextFlow pipeline will improve 
analysis interpretation through presentation of raw code and order of processes. 

## Usage

`nextflow run main.nf -entry NF_GWAS -with-singulariy gwas-nf.sif`

## ToDo
	1. add harmonicDist() and write outlier gene list and manhattan plots
	1. add snpeff to container for process invoking
	2. add entap to container for process invoking
