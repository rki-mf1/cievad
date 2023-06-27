
rule vcf_index:
    input:
        config["HEAD_DIR"] + "/data/simulated_hap{sample}/callset.vcf.gz"
    output:
        config["HEAD_DIR"] + "/data/simulated_hap{sample}/callset.vcf.gz.tbi"
    wrapper:
        "v2.0.0/bio/bcftools/index"


rule vcf_evaluation:
    input:
        truthset = config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.normalized.sorted.vcf.gz",
        truthidx = config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.normalized.sorted.vcf.gz.tbi",
        callset  = config["HEAD_DIR"] + "/data/simulated_hap{sample}/callset.vcf.gz",
        callidx  = config["HEAD_DIR"] + "/data/simulated_hap{sample}/callset.vcf.gz.tbi"
    params:
        prefix   = config["HEAD_DIR"] + "/data/simulated_hap{sample}/eval.picard"
    output:
        detail   = config["HEAD_DIR"] + "/data/simulated_hap{sample}/eval.picard.variant_calling_detail_metrics",
        summary  = config["HEAD_DIR"] + "/data/simulated_hap{sample}/eval.picard.variant_calling_summary_metrics"
    conda:
        config["HEAD_DIR"] + "/env/conda_picard.yaml"
    log:
        config["HEAD_DIR"] + "/logs/picard/evaluation_hap{sample}.log"
    shell:
        """
            picard CollectVariantCallingMetrics \
                --INPUT {input.callset} \
                --DBSNP {input.truthset} \
                -O {params.prefix}
        """


rule report:
    input:
        expand(config["HEAD_DIR"] + "/data/simulated_hap{sample}/eval.picard.variant_calling_summary_metrics", sample=config["SAMPLES"])
    output:
        config["HEAD_DIR"] + "/results/variant_calling_summary_ngs"
    params:
        head_dir = config["HEAD_DIR"]
    shell:
        """
            sh {params.head_dir}/scripts/eval/picard_summary_of_summaries.sh {params.head_dir} > {output}
        """