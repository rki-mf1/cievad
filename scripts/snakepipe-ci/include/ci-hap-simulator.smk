
rule ci_hap_simulator:
    input:
        ref   = config["REF"]
    output:
        fasta = config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}.fasta",
        vcf   = config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}.vcf",
    conda:
        config["HEAD_DIR"] + "/env/conda_ci.yaml"
    log:
        config["HEAD_DIR"] + "logs/ci/mason_variator.log"
    shell:
        """
            mason_variator \
                --in-reference {input.ref} \
                --out-fasta {output.fasta} \
                --out-vcf {output.vcf} \
                --seed {wildcards.sample} \
                --snp-rate 0.01 \
                --small-indel-rate 0.005 \
                --min-small-indel-size 1 \
                --max-small-indel-size 20 \
                --sv-indel-rate 0 \
                --sv-inversion-rate 0 \
                --sv-translocation-rate 0 \
                --sv-duplication-rate 0 \
                2>> {log}
        """