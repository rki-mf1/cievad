
rule nanopore_read_simulator:
    input:
        simu_hap = config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.fasta",
        model    = expand(config["HEAD_DIR"] + "/" + config["MODEL_PREFIX"] + "{f}", f=["_aligned_reads.pkl", "_aligned_region.pkl", "_chimeric_info", "_error_markov_model", "_error_rate.tsv", "_first_match.hist", "_gap_length.pkl", "_ht_length.pkl", "_ht_ratio.pkl", "_match_markov_model", "_model_profile", "_reads_alignment_rate", "_strandness_rate", "_unaligned_length.pkl"])
    output:
        reads  = config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.nanopore_aligned_reads.fasta",
        errors = config["HEAD_DIR"] + "/data/simulated_hap{sample}/simulated.nanopore_aligned_error_profile"
    params:
        dna_type      = config["DNA_TYPE"],
        model_caller  = config["MODEL_CALLER"],
        model_prefix  = config["MODEL_PREFIX"],
        median_length = config["MEDIAN_LENGTH"],
        sd_length     = config["SD_LENGTH"],
        nb_reads      = config["NB_READS"],
        seed          = config["SEED"]
    threads:
        workflow.cores
    conda:
        config["HEAD_DIR"] + "/env/conda_nanosim.yaml"
    log:
        config["HEAD_DIR"] + "/logs/nanosim/nanosim.hap{sample}.log"
    shell:
        """
            simulator.py genome \
                -dna_type {params.dna_type} \
                -rg {input.simu_hap} \
                -c {params.model_prefix} \
                --b {params.model_caller} \
                -med {params.median_length} \
                -sd {params.sd_length} \
                -n {params.nb_reads} \
                -o data/simulated_hap{wildcards.sample}/simulated.nanopore \
                --seed {params.seed} \
                -t {threads} \
                > {log}
        """
