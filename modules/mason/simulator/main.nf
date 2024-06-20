process MASON_SIMULATOR {

    // Job label
    // tag "${sample}"

    // Store results
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.{NGSWGS.R1.fastq,NGSWGS.R2.fastq,bam}"

    // Engine settings
    conda 'bioconda::mason=2.0.9'

    // Resources
    cpus 2

    // Process I/O
    input:
    tuple val(sample), path(vcf)
    val ref
    val ref_idx

    output:
    path "simulated_hap${sample}.NGSWGS.{R1,R2}.fastq",     emit: fastqs
    path "simulated_hap${sample}.NGSWGS.bam",               emit: bam

    // Job script
    script:
    unique_seed = (params.seed * sample ) % 2147483647       // that's (2^31)-1, the upper bound for mason
    """
    mason_simulator \
        -ir ${ref} \
        -iv ${vcf} \
        -o simulated_hap${sample}.NGSWGS.R1.fastq \
        -or simulated_hap${sample}.NGSWGS.R2.fastq \
        -oa simulated_hap${sample}.NGSWGS.bam \
        --seed ${unique_seed} \
        --num-threads ${task.cpus} \
        --num-fragments ${params.nb_frag} \
        --fragment-min-size ${params.fragment_min_size} \
        --fragment-max-size ${params.fragment_max_size} \
        --fragment-mean-size ${params.fragment_mean_size} \
        --fragment-size-std-dev ${params.fragment_size_std_dev} \
        --illumina-read-length ${params.illumina_read_length}
    """


}
