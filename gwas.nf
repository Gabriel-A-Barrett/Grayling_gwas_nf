#!/usr/bin/env nextflow

Channel.fromPath( "data/vcf/test.vcf" )
    .into {fb_filtering_ch; prepEnv_vcf_ch}
lfmm_ch = Channel.fromPath ( "data/gwas/LEA/lfmm*" )
group_ch = Channel.fromPath ("data/gwas/BayPass/bp_fullgroup.txt")
meta_ch = Channel.fromPath ( "data/meta.csv")


process variant_filtering {
    input:
    path vcf from fb_filtering_ch
    path meta from params.meta
    path probSamples from params.probSamples

    output:
    //path 'consensus.recode.vcf' into clean_vcf_ch
    tuple path("consensus.vcf"), path("final_indv.txt") into clean_vcf_ch, vcf2baypass_ch
    path ("final_indv.txt") into vcfList_ch

    script:
    template 'fb_F_SITE-PIPELINE.sh'

    stub:
    """
    touch consensus.vcf.gz
    """
}

// split meta into env files for groups w/ data a.k.a genocide
process meta2grps {
    label 'little_demon'

    input:
    path meta from params.meta
    path writegrps_fun from params.writegrps_function
    path vcfList from vcfList_ch
    output:
    path "*.grp" into grp2baypass_ch, grp2meta2env_ch

    script:
    """
    #!/usr/bin/env Rscript
    source("${writegrps_fun}")
    writegrps("${meta}","${vcfList}")
    """
}

process vcf2baypass_ch {
    cache 'lenient'

    input:
    tuple path(vcf), path(indv) from vcf2baypass_ch
    each grp from grp2baypass_ch.flatten()
    //tuple path(vcf), path(indvs) from vcf2baypass_ch
    
    output:
    path ("*.vcf") into vcf2rda_ch
    path ("*.ped") into ped2lea_ch, ped2prune_ch
    path ("*.frq") into frq2bay_ch
    tuple val (grp), path ("*.indv"), path ("*.pop_order") into reorder_bay2env_ch
    path ("*.loci") into baypass_loci_ch
    
    script:
    """
    env="\$(basename ${grp} .grp)"
    awk -F"\t" '{print \$1}' ${grp} > "\${env}_indv.list"
    bgzip ${vcf} 
    bcftools view -I --force-samples -S "\${env}_indv.list" ${vcf}.gz > "\${env}_${vcf}"
    bcftools query -l "\${env}_${vcf}" > "\${env}.indv"
    plink --vcf "\${env}_${vcf}" --recode --allow-extra-chr --out "\${env}"
    cat "\${env}_${vcf}" | perl ${baseDir}/bin/vcf2baypass.pl ${grp} "\${env}.frq"
    """
}

process meta2env {
    cache 'lenient'
    input:
    path (meta2env) from params.meta2env_function
    path (meta) from params.meta
    tuple val (grps), path (indv) , path (grp_order) from reorder_bay2env_ch
    path (grp) from grp2meta2env_ch.collect()

    output:
    //tuple path("*.lfmm"), path("*.grp")
    path ("*.env") into env2lea_ch
    path ("*.bayenv") into env2bay_ch

    script:
    """
    #!/usr/bin/env Rscript
    source("${meta2env}")
    meta2env("${meta}", "${indv}", "${grps}", "${grp_order}")
    """
}


lea_ch = ped2lea_ch
    .mix ( env2lea_ch )
    .map { env -> 
        def key = env.name.toString().tokenize('.').get(0)
        return tuple(key, env)}
    .groupTuple(size:2)

process lea {
    label 'gwas'
    cache 'lenient'
    cpus "${task.cpus}"
    publishDir "${params.outdir}/gwas/output_files", mode: 'copy'

    input:
    // [env, [vcf,env],]
    tuple val(env), path(files) from lea_ch
    path fst_func from params.fst_function
    path lea_func from params.lea_function

    output:
    tuple val (env), path ("*_zscores.txt") into lea_zscores_ch
    stdout result

    script:
    """
    #!/usr/bin/env Rscript
    library("LEA", lib="${params.rlibrariesPath}")
    library("dplyr", lib="${params.rlibrariesPath}")
    source("${lea_func}")
    source('${fst_func}')
    lea("${env}","${files[0]}","${files[1]}")
    """
}

baypass_ch = frq2bay_ch
    .mix(env2bay_ch)
    .map { env -> 
        def key = env.name.toString().tokenize('.').get(0)
        return tuple(key, env)}
    .groupTuple(size:2)

channel.from(1..3).into {baypass_iter_ch; maxiter_ch}

maxiter_ch.max(). into {rep_ch; repgroup_ch}

process baypass {

    input:
    tuple val (env), path (files) from baypass_ch
    each rep from baypass_iter_ch

    output:
    tuple val (env), path ("baypass_core_${env}_${rep}_summary_{betai_reg,pi_xtx}.out") into baypass_rawoutput_ch, bay2spliter_ch
    
    script:
    """
    EPOP="\$(head -n 1 "${files[1]}" | awk '{print NF}')"
    POP="\$(head -n 1 "${files[0]}" | awk '{print NF}')"
    if ["\$EPOP" != "\$POP"]; then 
        echo "different number of populations in environment and genetic dataframes"
    else
        baypass -npop \${EPOP} -gfile "${files[0]}" -efile "${files[1]}" \
	    -nval 100 -thin 25 -burnin 50 -npilot 30 -pilotlength 10 \
	    -outprefix baypass_core_${env}_${rep} \
	    -nthreads ${task.cpus} -seed "\${RANDOM}"
    fi
    """
}

bay2merge_ch = bay2spliter_ch.transpose().groupTuple(size:6)

process merge_baypass {

    input:
    tuple val (env), path (files) from bay2merge_ch
    path (loci) from baypass_loci_ch
    path (merge_bay) from params.merge_baypass_function
    val (maxiter) from rep_ch

    output:
    tuple val (env), path ("${env}_med_baypass.txt") into baypass_med_ch

    script:
    """
    #!/usr/bin/env Rscript
    library("dplyr", lib="${params.rlibrariesPath}")
    library("matrixStats", lib="${params.rlibrariesPath}")
    source("${merge_bay}")
    merge_baypass("${env}","${maxiter}", "${loci}")
    """
}

uni2multi_ch = baypass_med_ch
    .join(lea_zscores_ch)
    .view()

process composite_stat {

    input:
    tuple val (env), path (bay), path (lea) from uni2multi_ch
    tuple path (ann_vcf), path (entap) from annotation_ch

    output:

    script:
    """
    """


}




/*
process rda {
    label 'gwas'
    cache 'lenient'

    input:
    tuple path(vcf), path(pop) from rda_ch

    output:
    stdout result
    //tuple( val(samplename), val(raw_filtered), file("${samplename}_${raw_filtered}_TSNEPlot.pdf")), emit: tsneplot_pdf
    script:
    """
    #!/usr/bin/env Rscript
    library("vcfR", lib="${params.rlibrariesPath}")
    library("vegan", lib="${params.rlibrariesPath}")
    source("${rda_func}")
    rda(${vcf},${pop},${env})
    """


}*/




