
rule ngs_read_simulator:
    input:
        simu_hap = config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.fasta"
    output:
        r1       = config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.ngs.R1.fastq",
        r2       = config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.ngs.R2.fastq"
#        ra = temp(config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.bam")
    params:
        nb_frag  = config["NGS_NB_FRAGS"]
    threads:
        workflow.cores
    conda:
        config["HEAD_DIR"] + "/env/conda_mason.yaml"
    log:
        config["HEAD_DIR"] + "/logs/mason/mason_simulator.hap{sample}.log"
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

#rule ngs_read_sort:
#    input:
#        rules.ngs_simulator.output.ra
#    output:
#        config["HEAD_DIR"] + "/data/simulated_hap{samples}/simulated.sorted.bam"
#    conda:
#        config["HEAD_DIR"] + "/env/conda_bwa_and_samtools.yaml"
#    shell:
#        """
#            samtools sort -o {output} {input}
#        """

#rule ngs_read_index:
#    input:
#        rules.ngs_read_sort.output
#    output:
#        config["HEAD_DIR"] + "/data-ci/simulated_hap{samples}/simulated.sorted.bam.bai"
#    conda:
#        config["HEAD_DIR"] + "/env/conda_bwa_and_samtools.yaml"
#    shell:
#        """
#            samtools index {input}
#        """