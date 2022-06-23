include { prepGWAS } from '../subworkflows/preGWAS.nf'
include { conductGWAS } from '../subworkflows/conductGWAS.nf'
include { postGWAS } from '../subworkflows/postGWAS.nf'

bay_iter = channel.from(1..3)

workflow GWAS {

    main:
    prepGWAS ( )

    conductGWAS ( prepGWAS.out.lea, prepGWAS.out.bay, prepGWAS.out.loci )

    postGWAS ( conductGWAS.out.gwas_output )

}


