from snakemake.utils import min_version
min_version("6.0")

module bcftools:
    snakefile:
        config["HEAD_DIR"] + "/snakemake/modules/bcftools/Snakefile"


use rule bcftools_norm_noref from bcftools as norm_callset with:
    input:
        vcf = config["HEAD_DIR"] + "/data/simulated_hap{sample}/callset.vcf.gz"
    output:
        vcf = temp(config["HEAD_DIR"] + "/data/simulated_hap{sample}/callset.normalized.vcf.gz")
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"
    log:
        config["HEAD_DIR"] + "/logs/bcftools_norm.hap{sample}.callset.log"


use rule bcftools_sort from bcftools as sort_callset with:
    input:
        vcf = rules.norm_callset.output
    output:
        vcf = config["HEAD_DIR"] + "/data/simulated_hap{sample}/callset.normalized.sorted.vcf.gz"
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"
    log:
        config["HEAD_DIR"] + "/logs/bcftools/bcftools_sort.hap{sample}.callset.log"


use rule bcftools_index from bcftools as index_callset with:
    input:
        vcf = rules.sort_callset.output
    output:
        tbi = config["HEAD_DIR"] + "/data/simulated_hap{sample}/callset.normalized.sorted.vcf.gz.tbi"
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"