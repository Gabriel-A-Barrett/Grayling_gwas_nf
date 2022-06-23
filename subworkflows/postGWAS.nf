include {compositeStats; combineEnv} from '../modules/composite_stats'

workflow postGWAS {
    take:
    gwas_output
    
    main:
    compositeStats (gwas_output, params.ann_vcf, params.fullentap, params.headers_key)

    combineEnv (compositeStats.out.candidates.collect())

}