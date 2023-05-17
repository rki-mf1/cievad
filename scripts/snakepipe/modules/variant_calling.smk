
rule freebayes:
    input:
        config["HEAD_DIR"] + "/data/alignment/simulated_reads_to_original_ref.sorted.bam"
    log:
        config["HEAD_DIR"] + "/logs/freebayes.log"
    conda:
        config["HEAD_DIR"] + "/env/conda_freebayes_and_bcftools.yaml"
    params:
        mac=10,
        maf=0.1,
        mic=20
    output:
        config["HEAD_DIR"] + "/data/variant_calling/callset.freebayes.vcf"
    shell:
        """
            freebayes \
                -f {REF} \
                --min-alternate-count {params.mac} \
                --min-alternate-fraction {params.maf} \
                --min-coverage {params.mic} \
                --pooled-continuous \
                --haplotype-length -1 \
                {input} > \
                {output}
        """

rule freebayes_norm_vcf:
    input:
        rules.freebayes.output
    log:
        config["HEAD_DIR"] + "/logs/bcftools.norm.log"
    conda:
        config["HEAD_DIR"] + "/env/conda_freebayes_and_bcftools.yaml"
    params:
        ref=config["REF"]
    output:
        temp(config["HEAD_DIR"] + "/data/variant_calling/callset.freebayes.normalized.vcf")
    shell:
        """
            bcftools norm \
                --check-ref s \
                -f {params.ref} \
                -o {output} \
                {input}
        """

rule freebayes_sort_vcf:
    input:
        rules.freebayes_norm_vcf.output
    output:
        config["HEAD_DIR"] + "/data/variant_calling/callset.freebayes.normalized.sorted.vcf.gz"
    conda:
        config["HEAD_DIR"] + "/env/conda_freebayes_and_bcftools.yaml"
    shell:
        """
            bcftools sort \
                -o {output} \
                -O z \
                {input}
        """

rule freebayes_index_vcf:
    input:
        rules.freebayes_sort_vcf.output
    output:
        config["HEAD_DIR"] + "/data/variant_calling/callset.freebayes.normalized.sorted.vcf.gz.tbi"
    conda:
        config["HEAD_DIR"] + "/env/conda_freebayes_and_bcftools.yaml"
    shell:
        """
            bcftools index -t {output}
        """