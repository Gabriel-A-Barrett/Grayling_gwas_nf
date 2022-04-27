#!/usr/bin/env nextflow

process index {
    cache 'lenient'
    label 'index'
    // denote a list of modules with : seperator 
    module '/isg/shared/modulefiles/bwa/0.7.17'

    input:
    path ref from params.reference

    output:
    // error if not a glob pattern
    path 'index.*' into index_ch

    script:
    """
    bwa index -p index ${ref}
    """

    stub:
    """
    touch index.amb
    """
}

// set channel on pooled forward and reverse files
Channel
    .fromFilePairs( params.fq_files, size: 2 ) // size: number of inputs expected in grouping
    .ifEmpty { error "Cannot find any reads matching: ${params.fq_files}" }
    .set { read_pairs_ch }

//println("barcodes_${read_pairs_ch.replaceAll(/\d+/, "")}.txt") can't apply .replaceAll on channel type


process fastp {
    cache 'lenient'
    tag "fastp on ${pair_id}" // associates each process execution with a custom label
    module '/isg/shared/modulefiles/fastp/0.23.2' 

    publishDir "${params.outdir}/01-fastp_trimmed", mode: 'copy', pattern: "${pair_id}.html" // publish only html reports

    input:
    // takes shared string from paired reads [pair_id_1, [/path/to/file/,/path/to/file],...]
    tuple val(pair_id), path(reads) from read_pairs_ch

    output:
    tuple val(pair_id), path("trim_${pair_id}_{1,2}.fq.gz") into fastp_ch, fastp_ch1 // place output into a channel for tunneling to next process
    //path("${pair_id}.{html,log}") // for publishDir

    // execute fastp as if on the command line. Assumes bash
    script:
    """
    fastp -w ${task.cpus} -i ${reads[0]} -I ${reads[1]} \
    -o trim_${reads[0]} -O trim_${reads[1]} \
    --correction -h ${pair_id}.html &> ${pair_id}.log
    """

    stub:
    """
    touch trim_Golden{001..012}_{1,2}.fq.gz
    touch Golden{001..012}_{1,2}.{html,log}
    """
}

/*
read_pairs_ch
    .groupTuple(by:0)
    .map {prefix, reads -> tuple(prefix, reads.sort{it.name} ) }
    .view {"read pairs $it"}
Channel
    .from(1..12)
    .map {fq -> tuple("Goldenchr", path("/some/path/foo.${chr}.indels.vcf"), path("/other/path/foo.snvs.${chr}.vcf")}
    .view{ "numbers $it"}
*/
Channel
    .fromPath(params.barcodes).set {barcodes}

//fastp_ch1.view {"fastp Channel $it"}
// takes items from channel and puts them into a sorted list
//barcodes.toSortedList().view {" barcodes $it"}

// split reads into respective individual files
process demultiplex {
    cache 'lenient'
    tag "demultiplexing on pool: ${pair_id}"
    module '/isg/shared/modulefiles/stacks/2.53'
    
    //publishDir "${params.outdir}/01-process", mode: 'copy', patter: '*.log'

    input:
    tuple val(pair_id), path(trimmed_reads) from fastp_ch
    path barcode from barcodes.collect() // place all symlinks barcodes in directory for grabbing, outputs all items in the channel

    output:
    // Missing value declared as output parameter: ind_id
    tuple val(pair_id), path('Golden????.{1,2}.fq.gz') into demult_output_ch
    //path "*.log" into demultLog_ch

    script:
    """
    pool="\$(echo ${pair_id} | sed 's/Golden//')"
    process_radtags -i gzfastq -1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]} \
    -b barcodes_"\${pool}".txt -o ./ \
    -c -q -r \
    --inline_null \
    --renz_1 sbfI --renz_2 mseI    
    """

    // replaces the actual process when the -stub-run or -stub command line option is invoked 
    stub:
    """
    touch Golden{1..3}A{01..12}.{1,2}.fq.gz
    """
}

// function to determine if paired end or not from https://github.com/nf-core/demultiplex/blob/dev/workflows/demultiplex.nf
def getFastqPairName(fqfile) {
    // variable sampleName with regex pattern to filter for with .find()
    def sampleName = (fqfile =~ /.*\/(.+).[12]\.fq\.gz/)
    // creates a matcher for the input String
    if (sampleName.find()) {
        //  prints the pattern matched following regex matcher
        return sampleName.group(1)
    }
    return fqfile
}
demult_output_ch
    //.collect ()
    .flatten ()
    .filter { it =~/.*fq.gz/ }
    // creates a new tuple with matching value
    .map { fastq -> [getFastqPairName(fastq), fastq] }
    .groupTuple ()
    //.view {"filter demultiplex channel $it"}
    .set { demult_pruned_ch }

process bwa {
    cache 'lenient'
    tag "aligning $pair_id"
    
    module '/isg/shared/modulefiles/bwa/0.7.17'
    module '/isg/shared/modulefiles/samtools/1.9'

    input:
    // add literal string that is assigned the forward or reverse read and conveys that in the string meaning
        // forward_trim_read=contains("F")
    tuple val(pair_id), path(t_d_reads) from demult_pruned_ch
    path index from index_ch
    output:
    
    path("${pair_id}.bam") into all_aligned_reads_ch
    path("${pair_id}.stat") into stat_write_report_ch

    """
    # change index to explicit literal string
    bwa mem -t ${task.cpus} -R "@RG\\tID:${pair_id}\\tSM:${pair_id}" index ${t_d_reads[0]} ${t_d_reads[1]} | \
    samtools view -@ ${task.cpus} -S -h -u - | \
    samtools sort -@ ${task.cpus} - > ${pair_id}.bam
    samtools index -@ ${task.cpus} ${pair_id}.bam
    samtools stat -@ ${task.cpus} ${pair_id}.bam > ${pair_id}.stat
    """
}

// writes individual reports, need to concatanate reports maybe through collect and bash for loop
// remove bam records with low alignment rates before calling variants
process bwa_stats {
    tag "bwa Stats: ${bam}"
    label 'stats'

    publishDir "${params.outdir}/bwa", mode: 'copy', pattern: 'bamStats.tsv'
    
    module '/isg/shared/modulefiles/samtools/1.9'

    input:
    path bam_stat from stat_write_report_ch

    output:
    path 'mapRate.tsv' into bam_report_ch
    'bamStats.tsv'
    
    script:
    """
    # Extract Individual stats into 1 tab seperated file
    write_BamReport.sh ${bam_stat}
    """
    
    stub:
    """
    touch bamstats.tsv
    """
}

/**********************
* Remove Records w/ low mapping rates in channel
************************/
all_aligned_reads_ch // convert channel to [matcher, bam]
    .map { bam -> 
    def key = bam.name.toString().tokenize('.').get(0)
    return tuple(key, bam)}
    .set{key_all_aligned_ch}

bam_report_ch // Channel correction based on stats file from samtools stats 
    .splitCsv(sep:'\t', header:false, skip: 1) // converts file into channel formating (each row). else prints file literal string
    .map { // assign variable to column subsets and convert 2nd col. to class float
        def key = it[0].toString().tokenize('.').get(0) // similar to cut -d 
        def mappingrate = it[1].toFloat()
        [ key, mappingrate ]}
    .filter ({ key, mappingrate -> mapping_rate >= .75}) // retain samples with a mapping greater than 75%
    .join(key_all_aligned_ch) // outputs [key, stat, bam]
    .map { it[2] } // retain only bam records (3rd column)
    .set {cleaned_aligned_reads_ch}

process variant_calling {
    cache 'lenient'
    publishDir "${params.outdir}/fb", mode:'copy'
    module '/isg/shared/modulefiles/freebayes/1.3.4'
    
    input:
    // from https://github.com/brwnj/freebayes-nf/blob/master/main.nf
    path aln from cleaned_aligned_reads_ch.collect()
    path ref from params.reference
    path idx from index_ch.collect()
    
    output:
    path 'fb.vcf.gz' into fb_filtering_ch, fb_popmap_ch

    script:
    """
    # Noah: write text file using bash
    # ls *.bam > bamlist.txt
    freebayes -f ${ref} \
    ${aln.collect { "--bam $it" }.join(" ")} \
    -m 30 -q 20 --min-coverage 100 --skip-coverage 50000 | \
    bgzip -c > fb.vcf.gz
    """

    stub:
    """
    touch fb.vcf.gz
    """
}

process write_popmap_for_vcfFiltering {
    label 'little_demon'

    module '/isg/shared/modulefiles/bcftools/1.9'

    input:
    path meta from params.meta
    path vcf from fb_popmap_ch

    output:
    path 'popmap.tsv' into popmap_ch, envPrep_popmap_ch

    script:
    """
    # bcftools query -l $vcf > vcfsamples.txt
    # awk -F',' 'NR==FNR{a[\$1];next} \$4 in a {print \$4 \$7}' vcfsamples.txt $meta > popmap.tsv
    awk -F',' 'NR>1 {print \$4 \$7}' $meta > popmap.tsv
    """
}

process variant_filtering {

    tag "filtering variants from $vcf"
    module '/isg/shared/modulefiles/vcftools/0.1.16'
    module '/isg/shared/modulefiles/bcftools/1.9'

    publishDir "${params.outdir}/fb", mode: 'copy'

    input:
    path vcf from fb_filtering_ch
    path popmap from popmap_ch
    // before filtering remove problematic samples, removed low mapping here instead of in bwa_stats

    path probSamples from params.probSamples

    output:
    path 'consensus.recode.vcf' into clean_vcf_ch

    // incorporate 3rd flag with list of filters to paste. where to make the decision: within process or params

    script:
    """
    fb_F_SITE-PIPELINE.sh ${vcf} ${popmap} ${probSamples} #|
    #perl ${baseDir}/bin/filter_hwe_by_pop.pl -v - -p ${popmap} -h 0.001 -c 0.5 -o 'consensus'
    #rm *Kup* *Sag* *Col*
    """

    stub:
    """
    touch consensus.vcf.gz
    """
}

/*process metaToEnv {
    tag 'little demon'
    cache: 'leniant'

    input:
    path meta from params.meta
    path metaFormatr from params.metaFormatr

    output:
    path 'BayPass*env' into bpBF_ch
    path 'BayPass_group.txt' into group_ch
    path "lfmm*env" into lfmm_ch

    // standardizing input environmental data and converting to BayPass (pop.) and LEA (ind.) env. formats
    script:
    """
    #!/bin/Rscript
    source("metaFormatr")
    metaFormatr("${meta}")
    """
}*/

/*Channel
*    .fromList ( ['lat','elev','pit','temp','abc','udl','cc'] )
*    .set {env_ch}
*/

// write a channel based on the env. lfmm's in the dir. 
lfmm_ch = Channel
    .fromPath ( params.lfmm )

// correct for individuals remaining in study
process vcf_EnvPrep {
    tag "Prepping $vcf for $env"
    label 'stats'

    publishDir "${params.outdir}/gwas/input_files", mode: 'copy'

    module '/isg/shared/modulefiles/bcftools/1.9'

    input:
    path vcf from clean_vcf_ch
    path lfmm from lfmm_ch
    path group from params.group
    path pop from envPrep_popmap_ch


    output:
    path "*.vcf" into vcf2RDA_ch 
    path "*.ped" into ped2LEA_ch, ped2Faststructure
    path "*.frq" into vcf2BayPass_ch
    path "*_m_lfmm.env" into env2LEA_ch
    path "*_popmap" into rda_popmap_ch

    script:
    """
    vcfToENV.sh ${vcf} ${group} ${lfmm}
    """
}

// create tuple using groovy
ped2LEA_ch
    .merge(env2LEA_ch)   
    .flatten()
    // writes [abc, vcf] and [abc, lfmm]; creates item in list that is shared ID
    .map { env -> 
        // collect String before '_' 
        def key = env.name.toString().tokenize('_').get(0)
        // print list with matcher
        return tuple(key, env) }
    // adds items to second list based on matcher
    .groupTuple(size:2)
    // specify channel for process reading 
    .set { lea_ch }

process lea {
    label 'gwas'
    cache 'lenient'

    publishDir "${params.outdir}/gwas/output_files", mode: 'copy'

    input:
    // [env, [vcf,env],]
    tuple val(env), path(files) from lea_ch
    path fst_func from params.fst_function
    path lea_func from params.lea_function

    output:
    path "*_GD_zscores.txt" into lea_GD
    path "*_EA_zscores.txt" into lea_EA

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

/*process rda {
    label: 'gwas'
    cache: 'lenient'

    input:
    path vcf from vcf2RDA_ch
    path env from env2RDA_ch
    path pop from rda_popmap_ch

    output:
    
    script:
    """
    #!/usr/bin/env Rscript
    library("vcfR", lib="${params.rlibrariesPath}")
    library("vegan", lib="${params.rlibrariesPath}")
    png("rda_plot.png")
    source("${rda_func}")
    rda(${vcf},${pop},${env})
    dev.off()
    """
}
*/

//vcfToBayPass_ch
//    .merge(envToBayPass_ch)

/*vcfToBayPass_ch
    .map {env -> 
    def key = env.name.toString().tokenize('.').get(0)
    return tuple(key, env) }
    .view { "$it" }
*/



