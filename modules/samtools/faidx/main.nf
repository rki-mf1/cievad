process SAMTOOLS_FAIDX {
    // Job label
    // tag "${sample}"

    // Engine settings
    conda 'bioconda::samtools=1.19.2'

    // Resources
    cpus 1

    // Process I/O
    input:
    path ref

    output:
    path "${ref}.fai",  emit: refidx

    // Job script
    """
    samtools faidx \
        ${ref} \
        -o ${ref}.fai
    """


}
