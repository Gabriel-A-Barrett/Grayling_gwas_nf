include {lea} from '../modules/lea'
include {baypass} from '../modules/baypass'



workflow conductGWAS {
    
    Channel.from(1..3).set {bay_iter}
    
    take:
    lea_input
    bay_input

    main:
    lea (lea_input)
    baypass (bay_input, bay_iter)
}