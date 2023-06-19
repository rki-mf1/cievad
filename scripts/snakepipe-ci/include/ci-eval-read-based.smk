
rule ci_evaluation_readbased:
    input:
        truthset = config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/simulated.normalized.sorted.vcf.gz",
        truthidx = config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/simulated.normalized.sorted.vcf.gz.tbi",
        callset  = config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/" + config["CALLSET"],
        callidx  = config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/" + config["CALLSET_IDX"]
    params:
        prefix   = config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/eval.picard.readbased"
    output:
        detail   = config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/eval.picard.readbased.variant_calling_detail_metrics",
        summary  = config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/eval.picard.readbased.variant_calling_summary_metrics"
    conda:
        config["HEAD_DIR"] + "/env/conda_picard.yaml"
    log:
        config["HEAD_DIR"] + "/logs/ci.evaluation.readbased.log"
    shell:
        """
            picard CollectVariantCallingMetrics \
                --INPUT {input.callset} \
                --DBSNP {input.truthset} \
                -O {params.prefix}
        """
