process MINIMAP2SAMTOOLS {

    // Job label
    tag "${fastq.getName()}"

    // Store results
    publishDir "${params.outdir}", mode: 'copy', pattern: "${fastq.getBaseName()}.sorted.*"

    // Engine settings
    conda 'bioconda::minimap2=2.28 bioconda::samtools=1.19.2'

    // Resources
    cpus 2

    // Process I/O
    input:
    path fastq
    val ref
    val ref_idx

    output:
    path "${fastq.getBaseName()}.sorted.bam"
    path "${fastq.getBaseName()}.sorted.bam.bai"

    // Job script
    script:
    """
    minimap2 \
        -ax map-ont \
        ${ref} \
        ${fastq} \
        > ${fastq.getBaseName()}.sam

    samtools view \
        -S -b -h \
        ${fastq.getBaseName()}.sam \
        > ${fastq.getBaseName()}.bam

    samtools sort \
        ${fastq.getBaseName()}.bam \
        > ${fastq.getBaseName()}.sorted.bam

    samtools index \
        ${fastq.getBaseName()}.sorted.bam
    """
}
