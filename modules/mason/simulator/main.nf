process MASON_SIMULATOR {

    // Job label
    // tag "${sample}"

    // Store results
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.NGSWGS.{R1,R2}.fastq"

    // Engine settings
    conda 'bioconda::mason=2.0.9'

    // Resources
    cpus 2

    // Process I/O
    input:
    val id
    path "simulated_hap${id}.fasta"         // haplotype sequence, e.g. simulated_1.fasta

    output:
    path "simulated_hap${id}.NGSWGS.{R1,R2}.fastq",    emit: fastq

    // Job script
    script:
    unique_seed = (params.seed * id) % 2147483647       // that's (2^31)-1, the upper bound for mason
    """
    mason_simulator \
        -ir simulated_hap${id}.fasta \
        -n ${params.nb_frag} \
        -o simulated_hap${id}.NGSWGS.R1.fastq \
        -or simulated_hap${id}.NGSWGS.R2.fastq \
        --seed ${unique_seed} \
        --num-threads ${task.cpus} \
        --fragment-min-size ${params.fragment_min_size} \
        --fragment-max-size ${params.fragment_max_size} \
        --fragment-mean-size ${params.fragment_mean_size} \
        --fragment-size-std-dev ${params.fragment_size_std_dev} \
        --illumina-read-length ${params.illumina_read_length}
    """


}
