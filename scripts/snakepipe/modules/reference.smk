
rule index_ref:
    input:
        config["REF"]
    output:
        expand(config["REF"] + ".{ext}", ext=["bwt", "pac", "ann", "amb", "sa"])
    conda:
        config["HEAD_DIR"] + "/env/conda_bwa_and_samtools.yaml"
    log:
        config["HEAD_DIR"] + "/logs/index_ref.log"
    shell:
        """
            bwa index {input}
        """