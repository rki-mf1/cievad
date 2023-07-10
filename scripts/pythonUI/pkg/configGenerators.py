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
