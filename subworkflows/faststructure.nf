include { fastStructure ; fastStructureVis } from '../modules/fastStructure'


workflow ped_faststructure_vis {

    take:
    ped
    listK

    main:
    fastStructure ( ped, listK )

}