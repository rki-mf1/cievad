
rule short_read_alignment:
    input:
        r1  = config["HEAD_DIR"] + "/data/simulated_reads/R1.{nb_frag}.fastq",
        r2  = config["HEAD_DIR"] + "/data/simulated_reads/R2.{nb_frag}.fastq",
        idx = expand(config["REF"] + ".{ext}", ext=["bwt", "pac", "ann", "amb", "sa"])
    output:
        config["HEAD_DIR"] + "/data/alignment/simulated_reads_to_original_ref.{nb_frag}.sorted.bam"
    conda:
        config["HEAD_DIR"] + "/env/conda_bwa_and_samtools.yaml"
    log:
        config["HEAD_DIR"] + "/logs/short_read_alignment.log"
    shell:
        """
            bwa mem \
                -t 1 \
                {REF} \
                {input.r1} \
                {input.r2} | \
            samtools view \
                -@ 1 \
                -Sb \
                - | \
            samtools sort \
                -@ 1 \
                -o {output} \
                -
        """

rule short_read_alignment_bam_index:
    input:
        rules.short_read_alignment.output
    output:
        config["HEAD_DIR"] + "/data/alignment/simulated_reads_to_original_ref.{nb_frag}.sorted.bam.bai"
    conda:
        config["HEAD_DIR"] + "/env/conda_bwa_and_samtools.yaml"
    log:
        config["HEAD_DIR"] + "/logs/short_read_alignment_bam_index.log"
    shell:
        """
            samtools index \
                {input}
        """