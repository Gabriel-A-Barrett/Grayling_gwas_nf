// Subworkflows
include { calcVisPop            } from '../subworkflows/populations'
include { prepGWAS              } from '../subworkflows/preGWAS.nf'
include { conductGWAS           } from '../subworkflows/conductGWAS.nf'
include { postGWAS              } from '../subworkflows/postGWAS.nf'
include { ped_faststructure_vis } from '../subworkflows/faststructure.nf'

bay_iter = channel.from(1..3)

workflow GWAS {

    main:
    
    if (params.calcPopSumStats) {
        calcVisPop ( )
    }
    
    prepGWAS ( )

    if (params.calcAdmixture) {
        ped_faststructure_vis ( prepGWAS.out.ped, params.modelcomplexityList )
    }
    
    conductGWAS ( prepGWAS.out.lea, prepGWAS.out.bay, prepGWAS.out.loci )

    postGWAS ( conductGWAS.out.gwas_output )


}