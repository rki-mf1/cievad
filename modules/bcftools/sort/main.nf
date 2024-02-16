process BCFOOTLS_SORT {
    // Job label
    // tag "${sample}"

    // Store results
    publishDir "${params.outdir}", mode: 'copy', pattern: "*.sorted.vcf.gz"

    // Engine settings
    conda 'bioconda::bcftools=1.19'

    // Resources
    cpus 1

    // Process I/O
    input:
    path vcffile

    output:
    path "${vcffile.getSimpleName()}.sorted.vcf.gz",    emit: sort_vcf

    // Job script
    """
        bcftools sort \
            -o ${vcffile.getSimpleName()}.sorted.vcf.gz \
            -O z \
            ${vcffile}
    """


}