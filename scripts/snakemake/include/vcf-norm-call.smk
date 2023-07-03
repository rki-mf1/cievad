from snakemake.utils import min_version
min_version("6.0")

module VCFNorm_workflow:
    snakefile:
        config["HEAD_DIR"] + "/scripts/snakemake/modules/Snakefile.VCFNorm"


use rule norm_no_ref_check from VCFNorm_workflow as norm_callset with:
    input:
        config["HEAD_DIR"] + "/data/simulated_hap{sample}/callset.vcf.gz"
    output:
        temp(config["HEAD_DIR"] + "/data/simulated_hap{sample}/callset.normalized.vcf.gz")
    log:
        config["HEAD_DIR"] + "/logs/ngs/vcf_norm.call.hap{sample}.log"
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"


use rule sort from VCFNorm_workflow as sort_callset with:
    input:
        rules.norm_callset.output
    output:
        config["HEAD_DIR"] + "/data/simulated_hap{sample}/callset.normalized.sorted.vcf.gz"
    log:
        config["HEAD_DIR"] + "/logs/bcftools/vcf_sort.call.hap{sample}.log"
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"


use rule index from VCFNorm_workflow as index_callset with:
    input:
        rules.sort_callset.output
    output:
        config["HEAD_DIR"] + "/data/simulated_hap{sample}/callset.normalized.sorted.vcf.gz.tbi"
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"