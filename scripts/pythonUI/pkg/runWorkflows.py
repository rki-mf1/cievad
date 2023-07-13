import os


def run_hap_simu(args):
    print("Running haplotype simulation...\n")

    os.system('snakemake -p --use-conda --cores 1 -s ' + args.snakefile)


def run_ngs_simu(args):
    print("Running NGS read simulation...\n")

    os.system('snakemake -p --use-conda --cores ' + str(args.threads) + ' -s ' + args.snakefile)


def run_nanopore_simu(args):
    print("Running Nanopore read simulation...\n")

    os.system('snakemake -p --use-conda --cores ' + str(args.threads) + ' -s ' + args.snakefile)


def run_vcf_eval(args):
    print("Running VCF file-based evaluation of variants...\n")

    os.system('snakemake -p --use-conda --cores 1 -s ' + args.snakefile)