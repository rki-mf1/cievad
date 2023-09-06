import os
import python.configGenerators


def run_hap(args):
    print("Running haplotype simulation...\n")

    if args.command == 'hap':
    
        if args.config is not None:

            os.system('snakemake -p --use-conda --cores 1 --configfile ' + args.config + ' -s ' + args.snakefile)

        else:

            python.configGenerators.generate_hap_config(args)

            os.system('snakemake -p --use-conda --cores 1 --configfile configs/snake_config_haplotype.yaml -s ' + args.snakefile)


def run_ngs(args):
    print("Running NGS read simulation...\n")

    if args.command == 'ngs':

        if args.config is not None:
            
            os.system('snakemake -p --use-conda --cores ' + str(args.threads) + '--configfile ' + args.config + ' -s ' + args.snakefile)

        else:

            python.configGenerators.generate_ngs_config(args)

            os.system('snakemake -p --use-conda --cores ' + str(args.threads) + ' --configfile configs/snake_config_ngs.yaml -s ' + args.snakefile)


def run_ampli(args):
    print("Running amplicon and NGS read simulation...\n")

    if args.command == 'ampli':
    
        if args.config is not None:

            os.system('snakemake -p --use-conda --cores 1 --configfile ' + args.config + ' -s ' + args.snakefile)

        else:

            python.configGenerators.generate_ampli_config(args)

            os.system('snakemake -p --use-conda --cores 1 --configfile configs/snake_config_amplicon.yaml -s ' + args.snakefile)


def run_nanopore(args):
    print("Running Nanopore read simulation...\n")

    if args.command == 'nano':

        if args.config is not None:

            os.system('snakemake -p --use-conda --cores ' + str(args.threads) + ' --configfile ' + args.config + ' -s ' + args.snakefile)

        else:

            python.configGenerators.generate_nanopore_config(args)

            os.system('snakemake -p --use-conda --cores ' + str(args.threads) + ' --configfile configs/snake_config_nanopore.yaml -s ' + args.snakefile)


def run_eval(args):
    print("Running VCF file-based evaluation of variants...\n")

    if args.command == 'eval':

        if args.config is not None:

            os.system('snakemake -p --use-conda --cores 1 --configfile ' + args.config + ' -s ' + args.snakefile)

        else:

            python.configGenerators.generate_eval_config(args)

            os.system('snakemake -p --use-conda --cores 1 --configfile configs/snake_config_eval.yaml -s ' + args.snakefile)