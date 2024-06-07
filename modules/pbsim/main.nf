process PBSIM {

    // Job label
    tag "${sample}"

    // Store results
    publishDir "${params.outdir}", mode: 'copy', pattern: "simulated_hap${sample}.ONTWGS*.fastq"

    // Process I/O
    input:
    tuple val(sample), path(fasta)

    output:
    path "simulated_hap${sample}.ONTWGS_0001.fastq"

    // Job script
    script:
    unique_seed = (params.seed * sample) % 2147483647       // that's (2^31)-1, the upper bound for mason
    """
    ${projectDir}/${params.pbsim_bin} \
        --seed         ${unique_seed} \
        --strategy     ${params.pbsim_strategy} \
        --prefix       simulated_hap${sample}.ONTWGS \
        --genome       ${fasta} \
        --method       ${params.pbsim_method} \
        --qshmm        ${projectDir}/${params.pbsim_model} \
        --depth        ${params.pbsim_depth} \
        --length-max   ${params.pbsim_length_max} \
        --length-min   ${params.pbsim_length_min}
    """
}
