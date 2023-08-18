import os
import math
from pkg.myUtil import trim_trailing_slash


def generate_config_hap_simu(args):
    
    # check if the appropriate config directory already exists, create if not
    cfg_dir  = "configs/snakemake/simu/hap"
    dir_exists = os.path.exists(cfg_dir)
    if not dir_exists:
        os.makedirs(cfg_dir)
        print("New config directory is created at " + cfg_dir + "!")

    # generate config
    with open(cfg_dir + "/snake_config.yaml", 'w') as config :

        config.write('HEAD_DIR:\n')
        config.write('   ' + trim_trailing_slash(args.HEAD_DIR) + '\n\n')
        config.write('REF:\n')
        config.write('   ' + args.REF + '\n\n')
        config.write('SAMPLES:\n')

        padding = int(math.log2(args.nb_samples))
        for i in range(1, (args.nb_samples)+1):
            config.write('   - \"' + str(i).zfill(padding) + '\"\n')

        print("New config file is created at " + cfg_dir + "/snake_config.yaml\n")


def generate_config_ngs_simu(args):
    
    # check if the appropriate config directory already exists, create if not
    cfg_dir  = "configs/snakemake/simu/ngs"
    dir_exists = os.path.exists(cfg_dir)
    if not dir_exists:
        os.makedirs(cfg_dir)
        print("New config directory is created at " + cfg_dir + "!")

    # generate config
    with open(cfg_dir + "/snake_config.yaml", 'w') as config :

        config.write('HEAD_DIR:\n')
        config.write('   ' + trim_trailing_slash(args.HEAD_DIR) + '\n\n')
        config.write('NGS_NB_FRAGS:\n')
        config.write('   ' + str(args.nb_frags) + '\n\n')
        config.write('SAMPLES:\n')

        padding = int(math.log2(args.nb_samples))
        for i in range(1, (args.nb_samples)+1):
            config.write('   - \"' + str(i).zfill(padding) + '\"\n')

        print("New config file is created at " + cfg_dir + "/snake_config.yaml\n")


def generate_config_ampli_simu(args):
    
    # check if the appropriate config directory already exists, create if not
    cfg_dir  = "configs/snakemake/simu/amplicons"
    dir_exists = os.path.exists(cfg_dir)
    if not dir_exists:
        os.makedirs(cfg_dir)
        print("New config directory is created at " + cfg_dir + "!")

    # generate config
    with open(cfg_dir + "/snake_config.yaml", 'w') as config :

        config.write('HEAD_DIR:\n')
        config.write('   ' + trim_trailing_slash(args.HEAD_DIR) + '\n\n')
        config.write('REF:\n')
        config.write('   ' + args.REF + '\n\n')
        config.write('PRIMER:\n')
        config.write('   ' + args.PRIMER + '\n\n')
        config.write('SAMPLES:\n')

        padding = int(math.log2(args.nb_samples))
        for i in range(1, (args.nb_samples)+1):
            config.write('   - \"' + str(i).zfill(padding) + '\"\n')

        print("New config file is created at " + cfg_dir + "/snake_config.yaml\n")


def generate_config_eval(args):
    
    # check if the appropriate config directory already exists, create if not
    cfg_dir  = "configs/snakemake/eval/ngs/vcf"
    dir_exists = os.path.exists(cfg_dir)
    if not dir_exists:
        os.makedirs(cfg_dir)
        print("New config directory is created at " + cfg_dir + "!")

    # generate config
    with open(cfg_dir + "/snake_config.yaml", 'w') as config :

        config.write('HEAD_DIR:\n')
        config.write('   ' + trim_trailing_slash(args.HEAD_DIR) + '\n\n')
        config.write('SAMPLES:\n')

        padding = int(math.log2(args.nb_samples))
        for i in range(1, (args.nb_samples)+1):
            config.write('   - \"' + str(i).zfill(padding) + '\"\n')

        print("New config file is created at " + cfg_dir + "/snake_config.yaml\n")


def generate_config_nanopore_simu(args):
    
    # check if the appropriate config directory already exists, create if not
    cfg_dir  = "configs/snakemake/simu/nanopore"
    dir_exists = os.path.exists(cfg_dir)
    if not dir_exists:
        os.makedirs(cfg_dir)
        print("New config directory is created at " + cfg_dir + "!")

    # generate config
    with open(cfg_dir + "/snake_config.yaml", 'w') as config :

        config.write('HEAD_DIR:\n')
        config.write('   ' + trim_trailing_slash(args.HEAD_DIR) + '\n\n')
        config.write('MODEL_PREFIX:\n')
        config.write('   ' + str(args.model_prefix) + '\n\n')
        config.write('DNA_TYPE:\n')
        config.write('   ' + str(args.dna_type) + '\n\n')
        config.write('MODEL_CALLER:\n')
        config.write('   ' + str(args.model_caller) + '\n\n')
        config.write('MEDIAN_LENGTH:\n')
        config.write('   ' + str(args.median_length) + '\n\n')
        config.write('SD_LENGTH:\n')
        config.write('   ' + str(args.sd_length) + '\n\n')
        config.write('NB_READS:\n')
        config.write('   ' + str(args.nb_reads) + '\n\n')
        config.write('SEED:\n')
        config.write('   ' + str(args.seed) + '\n\n')
        config.write('SAMPLES:\n')

        padding = int(math.log2(args.nb_samples))
        for i in range(1, (args.nb_samples)+1):
            config.write('   - \"' + str(i).zfill(padding) + '\"\n')

        print("New config file is created at " + cfg_dir + "/snake_config.yaml\n")