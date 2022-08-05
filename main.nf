#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Gabriel-A-Barrett/grayling-gwas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GitHub : https://github.com/Gabriel-A-Barrett/grayling-gwas    
-----------------------------------------------------------------------
*/

<<<<<<< HEAD
log.info """\
         R N A T O Y   P I P E L I N E    
         =============================
         meta: ${params.meta}
         vcf : ${params.vcf}
         reads : ${params.reads}
         EnTAP: ${params.fullentap}
         """
         .stripIndent()

=======
>>>>>>> 5c0a2d26d78f250c7e7e3d4396a465f5996f244a

include { GWAS } from './workflows/GWAS.nf'

// when defining a workflow you must use -entry with nextflow run 
workflow NF_GWAS {    
    GWAS ()
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir\n" : "Oops .. something went wrong" )
}
