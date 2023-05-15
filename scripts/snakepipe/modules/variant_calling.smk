
rule freebayes:
    input:
        config["HEAD_DIR"] + "/data/alignment/simulated_reads_to_original_ref.sorted.bam"
    log:
        config["HEAD_DIR"] + "/logs/variantcalling.log"
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

rule bcftools_norm:
    input:
        rules.freebayes.output
    log:
        config["HEAD_DIR"] + "/logs/bcftools.norm.log"
    conda:
        config["HEAD_DIR"] + "/env/conda_freebayes_and_bcftools.yaml"
    params:
        ref=config["REF"]
    output:
        config["HEAD_DIR"] + "/data/variant_calling/callset.freebayes.normalized.vcf"
    shell:
        """
            bcftools norm \
                -f {params.ref} \
                -o {output} \
                {input}
        """
