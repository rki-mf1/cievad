process BCFOOTLS_NORM {
    // Job label
    // tag "${sample}"

    // Engine settings
    conda 'bioconda::bcftools=1.19'

    // Resources
    cpus 1

    // Process I/O
    input:
    each vcffile
    path ref

    output:
    path "${vcffile.getSimpleName()}.normalized.vcf",    emit: norm_vcf

    // Job script
    """
    bcftools norm \
        --fasta-ref ${ref} \
        --check-ref s \
        --multiallelics -both \
        -o ${vcffile.getSimpleName()}.normalized.vcf \
        ${vcffile}
    """


}