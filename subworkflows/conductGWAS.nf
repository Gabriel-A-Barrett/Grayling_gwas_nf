include {lea} from '../modules/lea'
include {baypass} from '../modules/baypass'
include {baypass_median} from '../modules/baypass_median'



workflow conductGWAS {
    
    Channel.from(1..3).set {bay_iter}
    maxiter = bay_iter.max()
    //Integer x = (maxiter.toInteger() * 2)
    
    take:
    lea_input
    bay_input
    loci

    main:
    lea (lea_input)
    baypass (bay_input, bay_iter)
        .transpose()
        .groupTuple(size:6) // turn into groovy variable based on maxiter
        .set {bay_output}
    baypass_median (bay_output, loci, maxiter) // [env, [betai_1, xtx_1, betai_2, xtx_2]]

    emit:
    gwas_output = baypass_median.out.med_bay
        .join(lea.out.lea_zscores)

}