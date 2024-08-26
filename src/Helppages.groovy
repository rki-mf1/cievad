class Helper {
    def helpEval(version,params){
        String c_green = "\033[0;32m";
        String c_reset = "\033[0m";
        String c_yellow = "\033[0;33m";
        String c_blue = "\033[0;34m";
        String c_red = "\u001B[31m";
        String c_dim = "\033[2m";
        log.info """
        ____________________________________________________________________________________________
        
        ${c_blue}Robert Koch Institute, Genome Competence Center${c_reset}

        Workflow: cievad (${version}) - evaluation of callsets

        ${c_yellow}Minimal Usage Examples:${c_reset}

        nextflow run eval.nf -profile local,conda --callsets_dir <path/to/callsets>
        or
        nextflow run eval.nf -profile local,conda --sample_sheet <path/to/sample_sheet>

        ${c_yellow}Input parameter (required):${c_reset}

        ${c_green} --callsets_dir ${c_reset} Directory containing variant callsets for evaluation (files of format: callset_<X>.vcf[.gz]), where <X> is the index of the corresponding truthset.
        OR
        ${c_green} --sample_sheet ${c_reset} Sample sheet (.csv) with the header "index,truthset,callset". Every following line contains an index and matching truth- and callset.

        ${c_yellow}Other workflow parameter:${c_reset}

        ${c_green} --outdir ${c_reset} directory to save results in [default: ${params.outdir}]
        """
    }

    def helpHap(version,params){
        String c_green = "\033[0;32m";
        String c_reset = "\033[0m";
        String c_yellow = "\033[0;33m";
        String c_blue = "\033[0;34m";
        String c_red = "\u001B[31m";
        String c_dim = "\033[2m";
        log.info """
        ____________________________________________________________________________________________
        
        ${c_blue}Robert Koch Institute, Genome Competence Center${c_reset}

        Workflow: cievad (${version}) - haplotype generation

        ${c_yellow}Minimal Usage Example:${c_reset}
        
        nextflow run hap.nf -profile local,conda --reference <cievad/path/to/ref>

        ${c_yellow}Input parameter (required):${c_reset}

        ${c_green} --reference ${c_reset} reference genome (.fasta) used for the generation of synthetic sequencing data

        ${c_yellow}Other workflow parameter:${c_reset}

        ${c_green} --n ${c_reset} number of synthetic samples to be generated [default: ${params.n}]
        ${c_green} --read_type ${c_reset} type of synthetic reads to be generated (options: ngs, ont) [default: ${params.read_type}]
        ${c_green} --outdir ${c_reset} directory to save results in [default: ${params.outdir}]
        ${c_green} --report_nstretches ${c_reset} flag to activate reporting of N stretches in the reference genome [default: ${params.report_nstretches}]

        ${c_yellow}Next Generation Sequencing parameter, optional if [--read_type ngs] is supplied ${c_reset}

        ${c_green} --nb_frag ${c_reset} number of fragments per sample [default: ${params.nb_frag}]
        ${c_green} --fragment_min_size ${c_reset} minimum size of fragments [default: ${params.fragment_min_size}]
        ${c_green} --fragment_max_size ${c_reset} maximum size of fragments [default: ${params.fragment_max_size}]
        ${c_green} --fragment_mean_size ${c_reset} mean size of fragments [default: ${params.fragment_mean_size}]
        ${c_green} --fragment_size_std_dev ${c_reset} standard deviation for fragment size [default: ${params.fragment_size_std_dev}]
        ${c_green} --illumina_read_length ${c_reset} read length of synthetic illumina reads [default: ${params.illumina_read_length}]
        ${c_green} --mason_additional_args ${c_reset} additional arguments for mason_simulator [default: ${params.mason_additional_args}]

        ${c_yellow}Nanopore Sequencing parameter, optional if [--read_type ont] is supplied ${c_reset}

        ${c_green} --dna_type ${c_reset} used DNA type (options: linear, circular) [default: ${params.dna_type}]
        ${c_green} --model_prefix ${c_reset} path and prefix of a NanoSim model [default: ${params.model_prefix}]
        ${c_green} --model_caller ${c_reset} algorithm to conduct the basecalling [default: ${params.model_caller}]
        ${c_green} --median_length ${c_reset} median length of the synthetic reads [default: ${params.median_length}]
        ${c_green} --sd_length ${c_reset} standard deviation of the synthetic read lengths [default: ${params.sd_length}]
        ${c_green} --nb_reads ${c_reset} number of synthetic reads per sample [default: ${params.nb_reads}]
        """
    }
}
