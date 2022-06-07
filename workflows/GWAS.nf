include { prepGWAS } from '../subworkflows/Prep_GWAS_Inputs.nf'
include { conductGWAS } from '../subworkflows/conductGWAS.nf'
//include { postGWAS } from '../subworkflows/'

bay_iter = channel.from(1..3)

workflow GWAS {

    main:
    prepGWAS ( )

    conductGWAS ( prepGWAS.out.lea, prepGWAS.out.bay)

   //postGWAS (conductGWAS.out)

}


