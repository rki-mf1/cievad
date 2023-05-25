
rule evaluation:
    input:
        truthset = config["HEAD_DIR"] + "/data/simulated_hap/simulated.normalized.sorted.vcf.gz",
        truthidx = config["HEAD_DIR"] + "/data/simulated_hap/simulated.normalized.sorted.vcf.gz.tbi",
        callset  = config["HEAD_DIR"] + "/data/variant_calling/callset.{nb_frag}.freebayes.normalized.sorted.vcf.gz",
        callidx  = config["HEAD_DIR"] + "/data/variant_calling/callset.{nb_frag}.freebayes.normalized.sorted.vcf.gz.tbi"
    params:
        prefix = config["HEAD_DIR"] + "/data/evaluation/eval.{nb_frag}.freebayes"
    output:
        config["HEAD_DIR"] + "/data/evaluation/eval.{nb_frag}.freebayes.variant_calling_detail_metrics",
        config["HEAD_DIR"] + "/data/evaluation/eval.{nb_frag}.freebayes.variant_calling_summary_metrics"
    conda:
        config["HEAD_DIR"] + "/env/conda_picard.yaml"
    log:
        config["HEAD_DIR"] + "/logs/evaluation.log"
    shell:
        """
            picard CollectVariantCallingMetrics \
                --INPUT {input.callset} \
                --DBSNP {input.truthset} \
                -O {params.prefix}
        """