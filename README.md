# imsmp-variant-calling-benchmark
Workflows for benchmarking and continuous integration of variant calls, primarily tailored towards the RKI Sars-Cov-2 genome reconstruction and IMS-MP project

## Contents:
1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Help](#help)


## Requirements:

| Requirement | Tested with |
| --- | --- |
| 64 bits POSIX-compliant operating system | Ubuntu 20.04.5 LTS |
| [conda](https://docs.conda.io/en/latest/) | vers. 23.3.1 |
| [snakemake](https://snakemake.readthedocs.io/en/stable/) | vers. 3.13.3 |

Prior to the installation make sure your system meets all the requirements. Conda is required to launch individual steps of the pipeline. Snakemake is the current choice of workflow language for this pipeline and can easily be installed via conda as well. <br>

⚠️ On the system this pipeline was developend on we faced issues with snakemake versions 7 and later due to their internal shell execution settings. If similar issues occur on your system please refer to the tested versions listed above. <br>
⚠️ An older version of bcftools did not split multi-allelic sites correctly when using _bcftools norm_. This issue was not observed in the more recent version 1.17. The provided conda environment is set up accordingly.

## Installation:

The pipeline only needs to be downloaded, no installation of any internal components required.
```
git clone https://github.com/rki-mf1/imsmp-variant-calling-benchmark.git
```

One way to install snakemake is to use the conda package manager:
```
conda create -y -n snakemake
conda activate snakemake
conda install -y -c bioconda snakemake
```

## Usage:

Running any of the benchmarking workflows consists of two main steps:
   1. creating a configuration _yaml_ file
   2. running a snakemake workflow

### Config files
The current setup to launch the workflows of this repository is to create or edit a configuration file in _yaml_ format for each of the Snakemake workflows.
All parameters required for a workflow that have to be set within the config file are listed is the corresponding Snakemake file's `PARAMS` section and can be found in this repository's [wiki](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/wiki) together with some more explanation. <br>

An example configuration file `scripts/snakepipe/snake_config.yaml` could look like the following:

```
HEAD_DIR:
   <path to this cloned repository>
REF:
   <path to my reference genome>
NB_FRAGMENTS:
   - 3000
   - 6000
```

⚠️ Without further ado, please create the configuration file within the same folder as the Snakemake worfklow file (explained below).

### Run snakemake
With a configuration file in place you can simply run the snakemake pipeline via:
```
snakemake -p --use-conda --cores 1 -s scripts/snakepipe/Snakemake
```

The `Snakemake` files, the files that define a workflow, are configured in a way s.t. they find the configuration file automatically if present in the same folder.

⚠️ Without further ado, please run the pipelines from the command line at the top folder of this project.

## Help:

Please visit the project [wiki](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/wiki) for more information, help and FAQs.
