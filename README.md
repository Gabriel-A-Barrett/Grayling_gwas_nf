# Arctic Grayling Landscape Genomics

nextflow pipeline for reproducing analysis done under Arctic Grayling

## Background

In an effort to support reproducibility of bioinformatic workflows and simply make re-running code easier in the event more samples come in, 
I strung together the Arctic grayling landscape genomics analysis into NextFlow. This NextFlow pipeline will improve 
analysis interpretation through presentation of raw code and order of processes. 

## Usage

`nextflow run main.nf -entry NF_GWAS -with-singulariy gwas-nf.sif`

## ToDo
* Add fastStructure and PCA Analysis and visualization
	* For nuetral, outlier, combined analysis
* add snpeff 
* add entap 
