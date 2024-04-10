process BCFOOTLS_INDEX {
    // Job label
    // tag "${sample}"

    // Store results
    //publishDir "${params.outdir}", mode: 'copy', pattern: "*.tbi"

    // Engine settings
    conda 'bioconda::bcftools=1.19'

    // Resources
    cpus 1

    // Process I/O
    input:
    path vcffile

    output:
    path "${vcffile.getName()}.tbi",    emit: index_vcf

    // Job script
    """
        bcftools index -t ${vcffile}
    """


}