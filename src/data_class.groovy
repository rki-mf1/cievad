class Helper {
    def helpEval(version){
        String c_green = "\033[0;32m";
        String c_reset = "\033[0m";
        String c_yellow = "\033[0;33m";
        String c_blue = "\033[0;34m";
        String c_red = "\u001B[31m";
        String c_dim = "\033[2m";
        log.info """
        ____________________________________________________________________________________________
        
        ${c_blue}Robert Koch Institute, MF1 Bioinformatics${c_reset}

        Workflow: cievad (${version})

        ${c_yellow}Minimal Usage Examples:${c_reset}

        nextflow run eval.nf -profile local,conda --callsets_dir <path/to/callsets>
        or
        nextflow run eval.nf -profile local,conda --sample_sheet <path/to/sample_sheet>

        ${c_yellow}Data Input, required:${c_reset}

        ${c_green} --callsets_dir ${c_reset} Directory containing variant callsets for evaluation (naming format: callset_<X>.vcf[.gz]).
        OR
        ${c_green} --sample_sheet ${c_reset} Sample sheet (.csv) with the header ("index","truthset","callset"), mapping corresponding truth- and callsets.
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
        
        ${c_blue}Robert Koch Institute, MF1 Bioinformatics${c_reset}

        Workflow: cievad (${version})

        ${c_yellow}Minimal Usage Example:${c_reset}
        
        nextflow run hap.nf -profile local,conda

        ${c_yellow}Individual Parameter, required:${c_reset}

        ${c_green} --n ${c_reset} number of synthetic samples to be generated
        ${c_green} --reference ${c_reset} reference used for the generation of synthetic sequencing data
        ${c_green} --read_type ${c_reset} type of resulting WGS synthetic reads (options: ngs, ont)

        ${c_yellow}Output Directory, required:${c_reset}

        ${c_green} --outdir ${c_reset} directory to save results in

        ${c_yellow}Next Generation Sequencing (WGS) Parameter, required if [--read_type ngs] supplied ${c_reset}

        ${c_green} --nb_frag ${c_reset} number of fragments per sample [default: ${params.nb_frag}]
        ${c_green} --fragment_min_size ${c_reset} minimum size of fragments [default: ${params.fragment_min_size}]
        ${c_green} --fragment_max_size ${c_reset} maximum size of fragments [default: ${params.fragment_max_size}]
        ${c_green} --fragment_mean_size ${c_reset} mean size of fragments [default: ${params.fragment_mean_size}]
        ${c_green} --fragment_size_std_dev ${c_reset} standard deviation for fragment size [default: ${params.fragment_size_std_dev}]
        ${c_green} --illumina_read_length ${c_reset} read length of synthetic illumina reads [default: ${params.illumina_read_length}]

        ${c_yellow}Nanopore Sequencing (WGS) Parameter, required if [--read_type ont] supplied ${c_reset}

        ${c_green} --dna_type ${c_reset} used DNA type [default: ${params.dna_type}]
        ${c_green} --model_prefix ${c_reset} path and prefix of the used model (e.g.: ${params.model_prefix})
        ${c_green} --model_caller ${c_reset} algorithm to conduct the basecalling [default: ${params.model_caller}]
        ${c_green} --median_length ${c_reset} median length of the resulting synthetic reads [default: ${params.median_length}]
        ${c_green} --sd_length ${c_reset} standard deviation length of the resulting synthetic reads [default: ${params.sd_length}]
        ${c_green} --nb_reads ${c_reset} number of synthetic reads [default: ${params.nb_reads}]
        """
    }
}