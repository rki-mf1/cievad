process MASON_VARIATOR {

    // Job label
    // tag "${sample}"

    // Store results
    publishDir "${params.outdir}", mode: 'copy', pattern: "simulated_hap*"

    // Engine settings
    conda 'bioconda::mason=2.0.9'

    // Resources
    cpus 1

    // Process I/O
    input:
    val sample
    val ref
    val ref_idx

    output:
    path "simulated_hap${sample}.fasta",    emit: fasta
    path "simulated_hap${sample}.vcf",      emit: vcf

    // Job script
    script:
    unique_seed = (params.seed * sample) % 2147483647       // that's (2^31)-1, the upper bound for mason_variator
    """
    mason_variator \
        --in-reference ${ref} \
        --out-fasta simulated_hap${sample}.fasta \
        --out-vcf simulated_hap${sample}.vcf \
        --seed ${unique_seed} \
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
