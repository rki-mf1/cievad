
rule survivor_simsv_parameterfile:
    output:
        config["HEAD_DIR"] + "/data/configs/parameter_file"
    log:
        config["HEAD_DIR"] + "/logs/survivor_simsv_parameterfile.log"
    shell:
        """
            {HEAD_DIR}/bin/SURVIVOR simSV {output}
        """

rule survivor_simsv:
    input:
        rules.survivor_simsv_parameterfile.output
    params:
        snp_freq = 0.01,
        outmode  = 0
    output:
        fasta = config["HEAD_DIR"] + "/data/simulated_hap/simulated.fasta",
        vcf   = config["HEAD_DIR"] + "/data/simulated_hap/simulated.vcf",
        ins   = config["HEAD_DIR"] + "/data/simulated_hap/simulated.insertions.fa",
        bed   = config["HEAD_DIR"] + "/data/simulated_hap/simulated.bed"
    log:
        config["HEAD_DIR"] + "/logs/survivor_simsv.log"
    shell:
        """
            {HEAD_DIR}/bin/SURVIVOR simSV \
                {REF} \
                {input} \
                {params.snp_freq} \
                {params.outmode} \
                {HEAD_DIR}/data/simulated_hap/simulated
        """

rule truthset_norm_vcf:
    input:
        rules.survivor_simsv.output.vcf
    log:
        config["HEAD_DIR"] + "/logs/truthset_norm_vcf.log"
    conda:
        config["HEAD_DIR"] + "/env/conda_freebayes_and_bcftools.yaml"
    params:
        ref = config["REF"]
    output:
        temp(config["HEAD_DIR"] + "/data/simulated_hap/simulated.normalized.vcf")
    shell:
        """
            bcftools norm \
                --check-ref s \
                -f {params.ref} \
                -o {output} \
                {input}
        """

rule truthset_sort_vcf:
    input:
        rules.truthset_norm_vcf.output
    output:
        config["HEAD_DIR"] + "/data/simulated_hap/simulated.normalized.sorted.vcf.gz"
    conda:
        config["HEAD_DIR"] + "/env/conda_freebayes_and_bcftools.yaml"
    shell:
        """
            bcftools sort \
                -o {output.vcfgz} \
                -O z \
                {input}
        """

rule truthset_index_vcf:
    input:
        rules.truthset_sort_vcf.output
    output:
        config["HEAD_DIR"] + "/data/simulated_hap/simulated.normalized.sorted.vcf.gz.tbi"
    conda:
        config["HEAD_DIR"] + "/env/conda_freebayes_and_bcftools.yaml"
    shell:
        """
            bcftools index -t {output}
        """