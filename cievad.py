# ----------------------------------------------------------------------------------------
# SETUP
# ----------------------------------------------------------------------------------------
import os
import sys
import argparse
import time
from python.runWorkflows import run_hap, run_ngs, run_ampli, run_nanopore, run_eval

if sys.version_info.major != 3:
    print("Error: Abort: This UI requires python3.")
    exit(1)


# ----------------------------------------------------------------------------------------
# PARSER
# ----------------------------------------------------------------------------------------
if __name__ == "__main__":
    __version_info__ = ('0','1','0')
    __version__ = '.'.join(__version_info__)

    parser = argparse.ArgumentParser(
            prog='cievad',
            description='CIEVaD - A tool suite to facilitate continuous integration and evaluation of variant detection.',
            epilog='For more help and bug reports please refer to the GitHub repository.')
    parser.add_argument('--version', action='version', version="%(prog)s ("+__version__+")")

    subparsers = parser.add_subparsers(help='sub-command help', dest='command')

    # ----------------
    # SUB PARSERS    |
    # ----------------

    # parser for haplotype simulation
    parser_hap = subparsers.add_parser('hap',
            help='Module to generate haplotypes from a given reference.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_hap.add_argument(
            '-s', '--snakefile',
            help='Path to the Snakefile.',
            default='snakemake/hap/Snakefile')
    parser_hap_group1 = parser_hap.add_argument_group('Run with config', 'Use a config file (yaml) to generate haplotypes.')
    parser_hap_group1.add_argument(
            '-c', '--config',
            metavar='FILE',
            default = None,
            help='Path to a config file for the snakemake pipeline.')
    parser_hap_group2 = parser_hap.add_argument_group('Run with parameter', 'Specify parameters to generate haplotypes.')
    parser_hap_group2.add_argument(
            '-d', '--head-dir',
            metavar='DIR',
            default= os.path.realpath(os.path.dirname(__file__)),
            help='Root directory path.')
    parser_hap_group2.add_argument(
            '-n', '--nb-samples',
            type = int,
            metavar='INT',
            default = 10,
            help='Specify the number of samples to be simulated.')
    parser_hap_group2.add_argument(
            '--seed',
            type = int,
            metavar='INT',
            default = int(round(time.time())),
            help='Specify a random seed. Default is current system time in seconds.')
    parser_hap_group2.add_argument(
            '-r', '--reference',
            metavar='FASTA',
            help='Path to reference genome.')
    parser_hap.set_defaults(func=run_hap)

    # parser for NGS read simulation
    parser_ngs = subparsers.add_parser('ngs',
            help='Module to generate NGS reads from a given reference.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_ngs.add_argument(
        '-s', '--snakefile',
        help='Path to the Snakefile.',
        default='snakemake/ngs/Snakefile')
    parser_ngs.add_argument(
        '-t', '--threads',
        help='Number of CPU threads for the task.',
        metavar='INT',
        default = 1)
    parser_ngs_group1 = parser_ngs.add_argument_group('Run with config', 'Use a config file (yaml) to generate NGS reads.')
    parser_ngs_group1.add_argument(
            '-c', '--config',
            metavar='FILE',
            default = None,
            help='Path to a config file for the snakemake pipeline.')
    parser_ngs_group2 = parser_ngs.add_argument_group('Run with parameter', 'Specify parameters to generate NGS reads.')
    parser_ngs_group2.add_argument(
            '-d', '--head-dir',
            metavar='DIR',
            default= os.path.realpath(os.path.dirname(__file__)),
            help='Root directory path.')
    parser_ngs_group2.add_argument(
            '-n', '--nb-samples',
            type = int,
            metavar='INT',
            default = 10,
            help='Specify the number of samples to be simulated.')
    parser_ngs_group2.add_argument(
            '--seed',
            type = int,
            metavar='INT',
            default = int(round(time.time())),
            help='Specify a random seed. Default is current system time in seconds.')
    parser_ngs_group2.add_argument(
            '-f', '--nb-frags',
            type = int,
            metavar='INT',
            default = 3000,
            help='Specify the number of genomic fragments used for the reads simulation. This INT*2 will result in the total number of NGS reads.')
    parser_ngs.set_defaults(func=run_ngs)

    # parser for generating amplicons and NGS reads
    parser_ampli = subparsers.add_parser('ampli',
            help='Module to generate amplicons and NGS reads from a given reference.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_ampli.add_argument(
            '-s', '--snakefile',
            help='Path to the Snakefile.',
            default='snakemake/amplicon/Snakefile')
    parser_ampli_group1 = parser_ampli.add_argument_group('Run with config', 'Use a config file (yaml) to generate amplicons and NGS reads.')
    parser_ampli_group1.add_argument(
            '-c', '--config',
            metavar='FILE',
            default = None,
            help='Path to a config file for the snakemake pipeline.')
    parser_ampli_group2 = parser_ampli.add_argument_group('Run with parameter', 'Specify parameters to generate amplicons and NGS reads.')
    parser_ampli_group2.add_argument(
            '-d', '--head-dir',
            metavar='DIR',
            default= os.path.realpath(os.path.dirname(__file__)),
            help='Root directory path.')
    parser_ampli_group2.add_argument(
            '-n', '--nb-samples',
            type = int,
            metavar='INT',
            default = 10,
            help='Specify the number of samples to be simulated.')
    parser_ampli_group2.add_argument(
            '--seed',
            type = int,
            metavar='INT',
            default = int(round(time.time())),
            help='Specify a random seed. Default is current system time in seconds.')
    parser_ampli_group2.add_argument(
            '-r', '--reference',
            metavar='FASTA',
            help='Path to reference genome.')
    parser_ampli_group2.add_argument(
            '-p', '--primers',
            metavar = 'BED',
            help='Path to primer file.')
    parser_ampli.set_defaults(func=run_ampli)

    # parser for generating nanopore reads
    parser_nanopore = subparsers.add_parser('nano',
            help='Module to generate Oxford-Nanopore-style long reads from a given reference.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_nanopore.add_argument(
            '-s', '--snakefile',
            help='Path to the Snakefile.',
            default='snakemake/nanopore/Snakefile')
    parser_nanopore.add_argument(
        '-t', '--threads',
        help='Number of CPU threads for the task.',
        metavar='INT',
        default = 1)
    parser_nanopore_group1 = parser_nanopore.add_argument_group('Run with config', 'Use a config file (yaml) to generate ONT-style long reads.')
    parser_nanopore_group1.add_argument(
            '-c', '--config',
            metavar='FILE',
            default = None,
            help='Path to a config file for the snakemake pipeline.')
    parser_nanopore_group2 = parser_nanopore.add_argument_group('Run with parameter', 'Specify parameters to generate ONT-style long reads.')
    parser_nanopore_group2.add_argument(
            '-d', '--head-dir',
            metavar='DIR',
            default= os.path.realpath(os.path.dirname(__file__)),
            help='Root directory path.')
    parser_nanopore_group2.add_argument(
            '-n', '--nb-samples',
            type = int,
            metavar='INT',
            default = 10,
            help='Specify the number of samples to be simulated.')
    parser_nanopore_group2.add_argument(
            '--seed',
            type = int,
            metavar='INT',
            default = int(round(time.time())),
            help='Specify a random seed. Default is current system time in seconds.')
    parser_nanopore_group2.add_argument(
            '-m', '--model-prefix',
            metavar='STR',
            default = 'aux/nanosim_model/human_NA12878_DNA_FAB49712_guppy/training',
            help='Specify a path (relative to the HEAD_DIR) to the prefix of a nanosim model.')
    parser_nanopore_group2.add_argument(
            '-g', '--model-caller',
            metavar='STR',
            default = 'guppy',
            help='Specify a caller of the nanosim model specified with -m.')
    parser_nanopore_group2.add_argument(
            '-y', '--dna-type',
            metavar='STR',
            default = 'linear',
            help='Specify a dna type for the nanosim simulator.')
    parser_nanopore_group2.add_argument(
            '-l', '--median-length',
            type = int,
            metavar='INT',
            default = 5000,
            help='Specify a median read length for the nanosim simulator.')
    parser_nanopore_group2.add_argument(
            '-a', '--sd-length',
            type = float,
            metavar='FLOAT',
            default = 1.05,
            help='Specify a standard deviation of the read length for the nanosim simulator.')
    parser_nanopore_group2.add_argument(
            '-r', '--nb-reads',
            type = int,
            metavar='INT',
            default = 180,
            help='Specify the number of long reads to be simulated per sample.')
    parser_nanopore.set_defaults(func=run_nanopore)

    # parser for variant evaluation
    parser_eval = subparsers.add_parser('eval',
            help='Module for variant set evaluation.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_eval.add_argument(
            '-s', '--snakefile',
            help='Path to the Snakefile.',
            default='snakemake/eval/Snakefile')
    parser_eval_group1 = parser_eval.add_argument_group('Run with config', 'Use a config file (yaml) to evaluate a variant callset.')
    parser_eval_group1.add_argument(
            '-c', '--config',
            metavar='FILE',
            default = None,
            help='Path to a config file for the snakemake pipeline.')
    parser_eval_group2 = parser_eval.add_argument_group('Run with parameter', 'Specify parameters to evaluate a variant callset.')
    parser_eval_group2.add_argument(
            '-d', '--head-dir',
            metavar='DIR',
            default= os.path.realpath(os.path.dirname(__file__)),
            help='Root directory path.')
    parser_eval_group2.add_argument(
            '-n', '--nb-samples',
            type = int,
            metavar='INT',
            default = 10,
            help='Specify the number of samples to be simulated.')
    parser_eval.set_defaults(func=run_eval)

    # ---------------
    # PARSE ARGS    |
    # ---------------

    args = parser.parse_args()
    args.func(args) if len(sys.argv)>1 else print("Error: Abort: Too few arguments. See help page: python vc_benchmark.py --help")