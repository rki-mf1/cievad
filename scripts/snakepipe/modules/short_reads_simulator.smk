
rule short_read_simulation_mason:
    input:
        config["HEAD_DIR"] + "/data/simulated_hap/simulated.fasta"
    params:
        mason = config["HEAD_DIR"] + "/bin/mason_simulator"
    output:
        r1 =      config["HEAD_DIR"] + "/data/simulated_reads/R1.fastq",
        r2 =      config["HEAD_DIR"] + "/data/simulated_reads/R2.fastq",
        ra = temp(config["HEAD_DIR"] + "/data/simulated_reads/simulated.bam")
    log:
        config["HEAD_DIR"] + "/logs/read_simulation_mason.log"
    shell:
        """
            {params.mason} \
                -ir {input} \
                -n 3000 \
                -o {output.r1} \
                -or {output.r2} \
                --out-alignment {output.ra} \
                --num-threads 1 \
                --fragment-min-size 450 \
                --fragment-max-size 550 \
                --fragment-mean-size 500 \
                --fragment-size-std-dev 20 \
                --illumina-read-length 150
        """

rule short_read_simulation_bam_sort:
    input:
        rules.short_read_simulation_mason.output.ra
    output:
        config["HEAD_DIR"] + "/data/simulated_reads/simulated.sorted.bam"
    conda:
        config["HEAD_DIR"] + "/env/conda_bwa_and_samtools.yaml"
    shell:
        """
            samtools sort -o {output} {input}
        """

rule short_read_simulation_bam_index:
    input:
        rules.short_read_simulation_bam_sort.output
    output:
        config["HEAD_DIR"] + "/data/simulated_reads/simulated.sorted.bam.bai"
    conda:
        config["HEAD_DIR"] + "/env/conda_bwa_and_samtools.yaml"
    shell:
        """
            samtools index {input}
        """