include { populations ; populationsVis } from "../modules/populations"

workflow calcVisPop {

    main:
    populations (params.ann_vcf, params.meta)//channel.from("Kuparuk","Colville"))// provides segmentation fault error, trying himem, believe it had to with the variants included in the vcf
    populationsVis (populations.out.popSumStats, populations.out.fstSumStats)

}