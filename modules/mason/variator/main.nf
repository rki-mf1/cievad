process MASON_VARIATOR {

    // Job label
    // tag "${sample}"

    // Engine settings
    conda 'bioconda::mason=2.0.9'

    // Resources
    cpus 1

    // Process I/O
    input:
    each sample
    file ref
    file ref_idx

    output:
    path "simulated_hap${sample}.fasta",    emit: fasta
    path "simulated_hap${sample}.vcf",      emit: vcf

    // Job script
    """
    mason_variator \
        --in-reference ${ref} \
        --out-fasta simulated_hap${sample}.fasta \
        --out-vcf simulated_hap${sample}.vcf \
        --seed ${params.seed} \
        --snp-rate 0.01 \
        --small-indel-rate 0.005 \
        --min-small-indel-size 1 \
        --max-small-indel-size 20 \
        --sv-indel-rate 0 \
        --sv-inversion-rate 0 \
        --sv-translocation-rate 0 \
        --sv-duplication-rate 0 \
        2> ${sample}.log
    """


}
