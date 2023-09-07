import math
from python.myUtil import trim_trailing_slash, mkdir_if_not_present


def generate_hap_config(args):
    
    # check if configs directory exists, create if not
    mkdir_if_not_present("configs")

    # generate config
    with open("configs/snake_config_haplotype.yaml", 'w') as config :

        config.write('HEAD_DIR:\n')
        config.write('   ' + trim_trailing_slash(args.head_dir) + '\n\n')
        config.write('REF:\n')
        config.write('   ' + args.reference + '\n\n')
        config.write('SAMPLES:\n')

        padding = int(math.log2(args.nb_samples))
        for i in range(1, (args.nb_samples)+1):
            config.write('   - \"' + str(i).zfill(padding) + '\"\n')

        print("New config file is created at configs/snake_config_haplotype.yaml\n")


def generate_ngs_config(args):
    
    # check if configs directory exists, create if not
    mkdir_if_not_present("configs")

    # generate config
    with open("configs/snake_config_ngs.yaml", 'w') as config :

        config.write('HEAD_DIR:\n')
        config.write('   ' + trim_trailing_slash(args.head_dir) + '\n\n')
        config.write('NGS_NB_FRAGS:\n')
        config.write('   ' + str(args.nb_frags) + '\n\n')
        config.write('SAMPLES:\n')

        padding = int(math.log2(args.nb_samples))
        for i in range(1, (args.nb_samples)+1):
            config.write('   - \"' + str(i).zfill(padding) + '\"\n')

        print("New config file is created at configs/snake_config_ngs.yaml\n")


def generate_ampli_config(args):
    
    # check if configs directory exists, create if not
    mkdir_if_not_present("configs")

    # generate config
    with open("configs/snake_config_amplicon.yaml", 'w') as config :

        config.write('HEAD_DIR:\n')
        config.write('   ' + trim_trailing_slash(args.head_dir) + '\n\n')
        config.write('REF:\n')
        config.write('   ' + args.reference + '\n\n')
        config.write('PRIMER:\n')
        config.write('   ' + args.primers + '\n\n')
        config.write('SAMPLES:\n')

        padding = int(math.log2(args.nb_samples))
        for i in range(1, (args.nb_samples)+1):
            config.write('   - \"' + str(i).zfill(padding) + '\"\n')

        print("New config file is created at configs/snake_config_amplicon.yaml\n")


def generate_nanopore_config(args):
    
    # check if configs directory exists, create if not
    mkdir_if_not_present("configs")

    # generate config
    with open("configs/snake_config_nanopore.yaml", 'w') as config :

        config.write('HEAD_DIR:\n')
        config.write('   ' + trim_trailing_slash(args.head_dir) + '\n\n')
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

        print("New config file is created at configs/snake_config_nanopore.yaml\n")


def generate_eval_config(args):
    
    # check if configs directory exists, create if not
    mkdir_if_not_present("configs")

    # generate config
    with open("configs/snake_config_eval.yaml", 'w') as config :

        config.write('HEAD_DIR:\n')
        config.write('   ' + trim_trailing_slash(args.head_dir) + '\n\n')
        config.write('SAMPLES:\n')

        padding = int(math.log2(args.nb_samples))
        for i in range(1, (args.nb_samples)+1):
            config.write('   - \"' + str(i).zfill(padding) + '\"\n')

        print("New config file is created at configs/snake_config_eval.yaml\n")
