# ----------------------------------------------------------------------------------------
# SETUP
# ----------------------------------------------------------------------------------------
import sys
import argparse
from pkg.configGenerators import generate_config_hap_simu, generate_config_ngs_simu, generate_config_eval, generate_config_nanopore_simu

if sys.version_info.major != 3:
    print("Error: Abort: This UI requires python3.")
    exit(1)


# ----------------------------------------------------------------------------------------
# PARSER
# ----------------------------------------------------------------------------------------
if __name__ == "__main__":
    __version_info__ = ('0','0','1')
    __version__ = '.'.join(__version_info__)

    parser = argparse.ArgumentParser(
            prog='vc_benchmark',
            description='This is a tool suite to facilitate variant benchmarking using the snakemake workflow language.',
            epilog='For more help and bug reports please refer to the GitHub repository.')
    parser.add_argument('--version', action='version', version="%(prog)s ("+__version__+")")

    subparsers = parser.add_subparsers(help='sub-command help')

    # parser for general arguments to be inherited by other parsers
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
            'HEAD_DIR',
            help='Root directory path.')
    parent_parser.add_argument(
            '-n', '--nb-samples',
            type = int,
            metavar='INT',
            default = 10,
            help='Specify the number of samples to be simulated.')

    # parser for generating a config for haplotype simulation
    parser_config_hap_simu = subparsers.add_parser('config-hap-simu',
            help='Module to create the configuration file for haplotype simulation.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            parents=[parent_parser], 
            aliases=['confhap'])
    parser_config_hap_simu.add_argument(
            'REF',
            help='Path to reference genome.')
    parser_config_hap_simu.set_defaults(func=generate_config_hap_simu)

    # parser for generating a config for NGS read simulation
    parser_config_ngs_simu = subparsers.add_parser('config-ngs-simu',
            help='Module to create the configuration file for NGS read simulation.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            parents=[parent_parser], 
            aliases=['confngs'])
    parser_config_ngs_simu.add_argument(
            '-f', '--nb-frags',
            type = int,
            metavar='INT',
            default = 3000,
            help='Specify the number of genomic fragments used for the reads simulation. This INT*2 will result in the total number of NGS reads.')
    parser_config_ngs_simu.set_defaults(func=generate_config_ngs_simu)

    # parser for generating a config for nanopore read simulation
    parser_config_nanopore_simu = subparsers.add_parser('config-nanopore-simu',
            help='Module to create the configuration file for Oxford-Nanopore-style long read simulation.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            parents=[parent_parser], 
            aliases=['confnano'])
    parser_config_nanopore_simu.add_argument(
            '-m', '--model-prefix',
            metavar='STR',
            default = 'aux/nanosim_model/human_NA12878_DNA_FAB49712_guppy/training',
            help='Specify a path (relative to the HEAD_DIR) to the prefix of a nanosim model.')
    parser_config_nanopore_simu.add_argument(
            '-c', '--model-caller',
            metavar='STR',
            default = 'guppy',
            help='Specify a caller of the nanosim model specified with -m.')
    parser_config_nanopore_simu.add_argument(
            '-d', '--dna-type',
            metavar='STR',
            default = 'linear',
            help='Specify a dna type for the nanosim simulator.')
    parser_config_nanopore_simu.add_argument(
            '-l', '--median-length',
            type = int,
            metavar='INT',
            default = 5000,
            help='Specify a median read length for the nanosim simulator.')
    parser_config_nanopore_simu.add_argument(
            '-s', '--sd-length',
            type = float,
            metavar='FLOAT',
            default = 1.05,
            help='Specify a standard deviation of the read length for the nanosim simulator.')
    parser_config_nanopore_simu.add_argument(
            '-r', '--nb-reads',
            type = int,
            metavar='INT',
            default = 180,
            help='Specify the number of long reads to be simulated per sample.')
    parser_config_nanopore_simu.add_argument(
            '--seed',
            type = int,
            metavar='INT',
            default = 479,
            help='Specify a random seed for the nanosim simulator.')
    parser_config_nanopore_simu.set_defaults(func=generate_config_nanopore_simu)

    # parser for generating a config for variant evaluation
    parser_config_eval = subparsers.add_parser('config-eval',
            help='Module to create the configuration file for variant evaluation.',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            parents=[parent_parser], 
            aliases=['confeval'])
    parser_config_eval.set_defaults(func=generate_config_eval)

    args = parser.parse_args()
    args.func(args) if len(sys.argv)>1 else print("Error: Abort: Too few arguments. See help page: python vc_benchmark.py --help")
