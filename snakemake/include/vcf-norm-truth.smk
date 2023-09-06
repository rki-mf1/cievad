from snakemake.utils import min_version
min_version("6.0")

module bcftools:
    snakefile:
        config["HEAD_DIR"] + "/snakemake/modules/bcftools/Snakefile"


use rule bcftools_norm from bcftools as norm_truthset with:
    input:
        vcf = config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.vcf",
        ref = config["REF"]
    output:
        vcf = temp(config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.normalized.vcf")
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"
    log:
        config["HEAD_DIR"] + "/logs/bcftools_norm.hap{sample}.truthset.log"


use rule bcftools_sort from bcftools as sort_truthset with:
    input:
        vcf = rules.norm_truthset.output
    output:
        vcf = config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.normalized.sorted.vcf.gz"
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"
    log:
        config["HEAD_DIR"] + "/logs/bcftools_sort.hap{sample}.truthset.log"


use rule bcftools_index from bcftools as index_truthset with:
    input:
        vcf = rules.sort_truthset.output
    output:
        tbi = config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.normalized.sorted.vcf.gz.tbi"
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"