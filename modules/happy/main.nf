process HAPPY {
    // Job label
    tag "${sample}"

    // Store results
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.sompy.*"

    // Engine settings
    conda 'bioconda::hap.py=0.3.15'

    // Resources
    cpus 1

    // Process I/O
    input:
    tuple val(sample), path(truthset), path(callset)
    val ref
    val ref_idx

    output:
    path "simulated_hap${sample}.sompy.stats.csv",      emit: csv
    path "simulated_hap${sample}.sompy.metrics.json",   emit: json

    // Job script
    """
    som.py \
        --no-fixchr-truth \
        --no-fixchr-query \
        --normalize-all \
        -r ${ref} \
        -o simulated_hap${sample}.sompy \
        ${truthset} \
        ${callset}
    """


}