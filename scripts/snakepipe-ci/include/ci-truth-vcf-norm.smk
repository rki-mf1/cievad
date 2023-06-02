rule ci_truth_vcf_norm:
    input:
        rules.ci_hap_simulator.output.vcf
    output:
        temp(config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/simulated.normalized.vcf")
    params:
        ref = config["REF"]
    conda:
        config["HEAD_DIR"] + "/env/conda_ci.yaml"
    log:
        config["HEAD_DIR"] + "/logs/ci_truth_vcf_norm.log"
    shell:
        """
            bcftools norm \
                --check-ref s \
                --multiallelics -both \
                -f {params.ref} \
                -o {output} \
                {input} 
        """

rule ci_truth_vcf_sort:
    input:
        rules.ci_truth_vcf_norm.output
    output:
        config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/simulated.normalized.sorted.vcf.gz"
    conda:
        config["HEAD_DIR"] + "/env/conda_ci.yaml"
    shell:
        """
            bcftools sort \
                -o {output} \
                -O z \
                {input}
        """

rule ci_truth_vcf_index:
    input:
        rules.ci_truth_vcf_sort.output
    output:
        config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/simulated.normalized.sorted.vcf.gz.tbi"
    conda:
        config["HEAD_DIR"] + "/env/conda_ci.yaml"
    shell:
        """
            bcftools index -t {input}
        """