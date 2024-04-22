process NANOSIM {

    // Job label
    tag "${sample}"

    // Store results
    publishDir "${params.outdir}", mode: 'copy', pattern: "simulated_hap${sample}.ONTWGS{_aligned_error_profile,_aligned_reads.fastq}"

    // Engine settings
    conda 'bioconda::nanosim=3.1.0'

    // Resources
    cpus 2

    // Process I/O
    input:
    tuple val(sample), path(fasta)

    output:
    path "simulated_hap${sample}.ONTWGS_aligned_error_profile",    emit: profile
    path "simulated_hap${sample}.ONTWGS_aligned_reads.fastq",      emit: fastq

    // Job script
    script:
    unique_seed = (params.seed * sample) % 2147483647       // that's (2^31)-1, the upper bound for mason
    """
    simulator.py genome \
        --ref_g        ${fasta} \
         -dna_type     ${params.dna_type} \
        --model_prefix ${projectDir}/${params.model_prefix} \
        --basecaller   ${params.model_caller} \
        --median_len   ${params.median_length} \
        --sd_len       ${params.sd_length} \
        --number       ${params.nb_reads} \
        --output       ${fasta.getSimpleName()}.ONTWGS \
        --seed         ${unique_seed} \
        --num_threads  ${task.cpus} \
        --fastq
    """
}
