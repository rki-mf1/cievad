process NANOSIM {

    // Job label
    // tag "${sample}"

    // Store results
    publishDir "${params.outdir}", mode: 'copy', pattern: "simulated_hap${id}.NANOWGS{_aligned_error_profile,_aligned_reads.fasta}"

    // Engine settings
    conda 'bioconda::nanosim=3.1.0'

    // Resources
    cpus 2

    // Process I/O
    input:
    val id
    path "simulated_hap${id}.fasta"         // haplotype sequence, e.g. simulated_1.fasta

    output:
    path "simulated_hap${id}.NANOWGS{_aligned_error_profile,_aligned_reads.fasta}",    emit: fastq

    // Job script
    script:
    unique_seed = (params.seed * id) % 2147483647       // that's (2^31)-1, the upper bound for mason
    """
    simulator.py genome \
        -dna_type ${params.dna_type} \
        -rg simulated_hap${id}.fasta \
        -c ${projectDir}/${params.model_prefix} \
        -b ${params.model_caller} \
        -med ${params.median_length} \
        -sd ${params.sd_length} \
        -n ${params.nb_reads} \
        -o simulated_hap${id}.NANOWGS \
        --seed ${unique_seed} \
        -t ${task.cpus}
    """
}
