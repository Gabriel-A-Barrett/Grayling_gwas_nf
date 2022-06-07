#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Gabriel-A-Barrett/grayling-gwas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GitHub : https://github.com/Gabriel-A-Barrett/grayling-gwas    
-----------------------------------------------------------------------
*/

// Parameters
//params.vcf = 
//params.meta = 

include { GWAS } from './workflows/GWAS.nf'

workflow NF_GWAS {    
    GWAS ()
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir\n" : "Oops .. something went wrong" )
}



// use -entry to narrow focus to NF_GWAS workflow