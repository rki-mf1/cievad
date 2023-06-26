
rule vcf_norm:
    input:
        rules.hap_simulator.output.vcf,
        ref=config["REF"]
    output:
        temp(config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.normalized.vcf")
    log:
        config["HEAD_DIR"] + "/logs/ngs/vcf_norm.hap{sample}.log"
    params:
        extra="--check-ref s --multiallelics -both"
    wrapper:
        "v2.0.0/bio/bcftools/norm"


rule vcf_sort:
    input:
        rules.vcf_norm.output
    output:
        config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.normalized.sorted.vcf.gz"
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"
    log:
        config["HEAD_DIR"] + "/logs/bcftools/vcf_sort.hap{sample}.log"
    shell:
        """
            bcftools sort \
                -o {output} \
                -O z \
                {input}
        """


rule vcf_index:
    input:
        rules.vcf_sort.output
    output:
        config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.normalized.sorted.vcf.gz.tbi"
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"
    shell:
        """
            bcftools index -t {input}
        """