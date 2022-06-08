//vcf = "/mnt/d/nextflow_testing/Grayling/freebayes/consensus.vcf"
//meta = "/mnt/d/nextflow_testing/Grayling/meta.csv"

include { vcfListIndividuals } from "../modules/vcf_list_individuals"
include { meta2grps } from "../modules/meta_2_grps"
include { vcf2GWAS } from "../modules/vcf_2_GWAS"
include { meta2env } from "../modules/meta_2_env"
//include { functions } from '../modules/fun_library' get message cannot find a component with the name 'function'

workflow prepGWAS {

    //take: // used to bring variables in the workflow invocation statement
    //vcf // call from parameters
    //meta 

    main:
    vcfListIndividuals ( params.vcf ) // get list of indv. in consesus.vcf 
    meta2grps ( params.meta, vcfListIndividuals.out.vcfIndv ) // get groups for baypass 
    vcf2GWAS ( params.vcf, meta2grps.out.grp2bay ) // subset vcf file and write baypass format
    meta2env ( params.meta, vcf2GWAS.out.grps2env ) // write environmental files and match baypass geno pop. order
    
    emit: // allows for prepGWAS.out 
    lea = getInputforGWAS(meta2env.out.env2lea, vcf2GWAS.out.ped2env)
    bay = getInputforGWAS(meta2env.out.env2bay, vcf2GWAS.out.frq2bay)
    loci = vcf2GWAS.out.baypass_loci
}

def getInputforGWAS(fib, lit) { // write [env, [geno_file, env_file]]
    
    env_ch = fib
    .mix(lit)
    .map { env -> 
        def key = env.name.toString().tokenize('.').get(0)
       return tuple (key, env)}
    .groupTuple(size:2)

    return env_ch
}
