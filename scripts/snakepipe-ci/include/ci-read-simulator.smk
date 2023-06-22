
rule ci_read_simulator:
    input:
        simu_hap = config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/simulated.fasta"
    output:
        r1       = config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/simulated.R1.fastq",
        r2       = config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/simulated.R2.fastq"
#        ra = temp(config["HEAD_DIR"] + "/data-ci/simulated_hap{sample}/simulated.bam")
    params:
        nb_frag = 3000
    threads:
        4
    conda:
        config["HEAD_DIR"] + "/env/conda_ci.yaml"
    log:
        config["HEAD_DIR"] + "/logs/ci_read_simulator.hap{sample}.log"
    shell:
        """
            mason_simulator \
                -ir {input} \
                -n {params.nb_frag} \
                -o {output.r1} \
                -or {output.r2} \
                --num-threads {threads} \
                --fragment-min-size 450 \
                --fragment-max-size 550 \
                --fragment-mean-size 500 \
                --fragment-size-std-dev 20 \
                --illumina-read-length 150
        """
#                --out-alignment {output.ra} \

#rule ci_simulated_read_sort:
#    input:
#        rules.ci_read_simulator.output.ra
#    output:
#        config["HEAD_DIR"] + "/data-ci/simulated_hap{samples}/simulated.sorted.bam"
#    conda:
#        config["HEAD_DIR"] + "/env/conda_ci.yaml"
#    shell:
#        """
#            samtools sort -o {output} {input}
#        """

#rule ci_simulated_read_index:
#    input:
#        rules.ci_simulated_read_sort.output
#    output:
#        config["HEAD_DIR"] + "/data-ci/simulated_hap{samples}/simulated.sorted.bam.bai"
#    conda:
#        config["HEAD_DIR"] + "/env/conda_ci.yaml"
#    shell:
#        """
#            samtools index {input}
#        """