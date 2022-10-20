# LandscapeGenomics_GWAS

Bioinformatics pipeline for population genomic analysis and visualization. the tool calculates population summary statistics and tests for genome-wide association and differentiation based on single nucleotide polymorphisms. 

## Pipeline Summary

1. Query individuals within VCF and calculate population summary statistics (process_radtags)
2. Write groups file (R)
2. Visualize population summary statistics (R)
3. Write individual corrected PED and VCF files (plink1.9,bcftools, vcf2baypass.pl)
4. Standardize and write environments (R)
5. Test Genome-Wide association and differentiation (BayPass, LEA)
6. Present candidates based on multivariate outlier detection (R) 
7. Combine across environment tests and evaluate for overlap (R)

## Usage

`nextflow run main.nf -c nextflow.config --meta samplesheet.csv --entap final_annotations_nocontam_lvl1.txt --vcf snpeff_annotated.vcf --headers_key ncbi_chromosome_ids.txt --first_env_column 3 --ggplot_indv_factors <indv_list.txt> -entry NF_GWAS -with-singulariy gwas-nf.sif`

| Column     | Description |
|------------| ----------- |
| id	     | Sample name. Must match VCF sample name program |
| population | A string describing the individual grouping |
| env        | A numerical value describing an environmental condition |

You will need to create a meta file that contains, in order, columns with individual ID's found within the VCF, Population, and an environmental variable or variables. You can also specify the column number in which environmental variables start and the tool will interpret that and every column afterwords as an environmental variable

Depends on final annotation output from EnTAP, snpEff VCF output, and chromosome key with chromosome id's in the first and numbers in the second column. 

Provide an optional list of individuals to order VIS off of. Manipulates population summary statistics VIS



