from snakemake.utils import min_version
min_version("6.0")

module VCFNorm_workflow:
    snakefile:
        config["HEAD_DIR"] + "/scripts/snakemake/modules/Snakefile.VCFNorm"


use rule norm from VCFNorm_workflow as norm_truthset with:
    input:
        vcf=rules.hap_simulator.output.vcf,
        ref=config["REF"]
    output:
        temp(config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.normalized.vcf")
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"
    log:
        config["HEAD_DIR"] + "/logs/bcftools/vcf_norm.truth.hap{sample}.log"


use rule sort from VCFNorm_workflow as sort_truthset with:
    input:
        rules.norm_truthset.output
    output:
        config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.normalized.sorted.vcf.gz"
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"
    log:
        config["HEAD_DIR"] + "/logs/bcftools/vcf_sort.truth.hap{sample}.log"


use rule index from VCFNorm_workflow as index_truthset with:
    input:
        rules.sort_truthset.output
    output:
        config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.normalized.sorted.vcf.gz.tbi"
    conda:
        config["HEAD_DIR"] + "/env/conda_bcftools.yaml"