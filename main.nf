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
    
    //publishDir "${params.outdir}/01-process", mode = 'copy', pattern = '*.log'

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

// function to determine if paired end or not from demultiplex github
def getFastqPairName(fqfile) {
    def sampleName = (fqfile =~ /.*\/(.+).[12]\.fq\.gz/)
    if (sampleName.find()) {
        return sampleName.group(1)
    }
    return fqfile
}
demult_output_ch
    //.collect ()
    .flatten ()
    .filter { it =~/.*fq.gz/ }
    .map { fastq -> [getFastqPairName(fastq), fastq] }
    .groupTuple ()
    //.view {"filter demultiplex channel $it"}
    .set { demult_pruned_ch }

process bwa {
    cache 'lenient'
    tag "aligned reads from $pair_id"
    
    module '/isg/shared/modulefiles/bwa/0.7.17'
    module '/isg/shared/modulefiles/samtools/1.9'

    input:
    // add literal string that is assigned the forward or reverse read and conveys that in the string meaning
        // forward_trim_read=contains("F")
    tuple val(pair_id), path(t_d_reads) from demult_pruned_ch
    path index from index_ch
    output:
    
    path("${pair_id}.bam") into stats_aligned_reads_ch, vcf_aligned_reads_ch

    """
    # change index to explicit literal string
    bwa mem -t ${task.cpus} -R "@RG\\tID:${pair_id}\\tSM:${pair_id}" index ${t_d_reads[0]} ${t_d_reads[1]} | \
    samtools view -@ ${task.cpus} -S -h -u - | \
    samtools sort -@ ${task.cpus} - > ${pair_id}.bam
    samtools index -@ ${task.cpus} ${pair_id}.bam
    """
}

// remove bam records with low alignment rates before calling variants
process bwa_stats {
    tag "aligned reads from ${bam}"
    label 'stats'

    publishDir "${params.outdir}/bwa", mode = 'copy', pattern = "bamstats.tsv"
    
    module '/isg/shared/modulefiles/samtools/1.9'

    input:
    path bam from stats_aligned_reads_ch

    output:
    path 'bamstats.tsv' into bamStats_ch
    //tuple "${bam}", "${mapping_rate}"

    script:
    """
    samtools stats ${bam} > ${bam.baseName}.stats
    // compile individual bam stat reports into 1 file
    write_BamReport.sh ${bam.baseName}.stats
    """
    
    stub:
    """
    touch bamstats.tsv
    """
}

/* __________________
* R E A D  R E P O R T
*----------------------
*/
// Channel correction based on stats file from samtools stats 
bamStats_ch
    // converts file into channel formating. else prints file literal string
    .splitCsv(sep:'\t', header:false, skip: 1)
    // assign variable to column subsets and convert 2nd col. to class float
    .map ({
        def bam = path(it[0])
        def mapping_rate = it[1].toFloat()
        [ bam, mapping_rate ]
    })
    .filter ({ bam, mapping_rate -> mapping_rate >= .75}) // retain samples with a mapping greater than 75%
    // decompose channel into something that looks like one column 
    .flatten()
    // retain only bam records
    .filter { it =~ '/Golden*.bam/' }
    .set {cleaned_aligned_reads_ch}

//aligned_reads_ch
    //.filter { it =~ /^Golden1A06.bam/|/^Golden1A08.bam/|/^Golden1B06.bam/|/^Golden1B07.bam/|/^Golden1C07.bam/|/^Golden1C08.bam/|/^Golden1D08.bam/|/^Golden1E06.bam/|/^Golden1E07.bam/|/^Golden1F07.bam/|/^Golden1G05.bam/ }

process variant_calling {
    cache 'lenient'
    publishDir params.outdir, mode:'copy'
    module '/isg/shared/modulefiles/freebayes/1.3.4'
    
    input:
    path aligned_reads from vcf_aligned_reads_ch // provides points to written bam files in work dir
    path ref from params.reference
    // from https://github.com/brwnj/freebayes-nf/blob/master/main.nf
    path aln from cleaned_aligned_reads_ch.collect() // bam ID's to provide to freebayes
    path idx from index_ch.collect()
    
    output:
    path 'fb.vcf.gz' into fb_filtering_ch, fb_popmap_ch

    script:
    """
    # Noah: write text file using bash
    # ls *.bam > bamlist.txt
    freebayes -f $ref \
    ${aln.collect { "--bam $it" }.join(" ")} \
    -m 30 -q 20 --min-coverage 100 --skip-coverage 50000 | \
    bgzip -c > fb.vcf.gz
    """

    stub:
    """
    touch fb.vcf.gz
    """
}

process write_popmap {
    label 'demon'

    input:
    path meta from params.meta
    path vcf from fb_popmap_ch

    output:
    path 'popmap.tsv' into popmap_ch

    script:
    """
    bcftools query -l $vcf > vcfsamples.txt
    awk -F ',' NR==FNR{a[\$1];next} \$1 in a {print \$5 \$8} vcfsamples.txt $meta > popmap.tsv
    """

}

process variant_filtering {

    tag "filtering variants from $vcf"
    module '/isg/shared/modulefiles/vcftools/0.1.16'
    module '/isg/shared/modulefiles/bcftools/1.9'

    input:
    path vcf from fb_filtering_ch
    path popmap from popmap_ch
    path probSamples from params.probSamples

    output:
    path 'consensus.vcf.gz' into clean_vcf_ch

    // incorporate 3rd flag with list of filters to paste. where to make the decision: within process or params

    script:
    """
    variant_filtering.sh ${vcf} ${popmap} ${probsamples} 
    """

    stub:
    """
    touch consensus.vcf.gz
    """

}

Channel
    .from ('lat','elev','pit','temp','abc','udl','cc')
    .set {env_ch}

process vcf_EnvPrep {
    tag "Prepping $vcf for $env"
    label 'stats'

    module '/isg/shared/modulefiles/bcftools/1.9'

    input:
    val env from env_ch
    path vcf from clean_vcf_ch
    // baypass populations
    path group from params.group
    // LEA env. input file
    path lfmm from params.lfmm

    output:
    path "${env}_${vcf}" into vcfToRDA_ch, vcftoLEA_ch
    path "*.frq" into vcfToBayPass_ch

    script:
    """
    vcfToEnv.sh ${vcf} ${env} ${group}
    """
}
