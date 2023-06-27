# imsmp-variant-calling-benchmark
<img alt="https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg?style=flat-square" src="https://img.shields.io/badge/snakemake-≥5.7.0-brightgreen.svg?style=flat-square">
Workflows for benchmarking and continuous integration of variant calls based on the Snakemake workflow language.

## Contents:
1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Help](#help)


## Requirements:

| Requirement | Tested with |
| --- | --- |
| 64 bits POSIX-compliant operating system | Ubuntu 20.04.5 LTS |
| [conda](https://docs.conda.io/en/latest/) | vers. 23.5.0 |
| [snakemake](https://snakemake.readthedocs.io/en/stable/) | vers. 7.25.3 |

Prior to the installation make sure your system meets all the requirements. Conda is required to launch individual steps of the pipeline. Snakemake is the current choice of workflow language for this pipeline and can easily be installed via conda as well. <br>


## Installation:

The pipeline only needs to be downloaded, no installation of any internal components required.
```
git clone https://github.com/rki-mf1/imsmp-variant-calling-benchmark.git
```

One way to install snakemake is to use the conda package manager:
```
conda create -y -n snakemake
conda activate snakemake
conda install -y -c bioconda snakemake">=5.7"
```

## Usage:

Running any of the benchmarking workflows consists of two main steps:
   1. creating a configuration _yaml_ file
   2. running a snakemake workflow

### Config files
The current setup to launch the workflows of this repository is to create or edit a configuration file in _yaml_ format for each of the Snakemake workflows.
All parameters required for a workflow that have to be set within the config file are listed is the corresponding Snakemake file's `PARAMS` section and can be found in this repository's [wiki](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/wiki/Confiluration-files) together with some more explanation. <br>

An example configuration file `snake_config.yaml` could look like the following:

```
HEAD_DIR:
   <path to this cloned repository>
REF:
   <path to my reference genome>
```

:construction: (TODO) Explain where configs have to go <br>
:construction: (TODO) There will likely be a simplification, e.g. a python module that generates this config in the right location 

### Run snakemake
With a configuration file in place you can simply run the snakemake pipeline via:
```
snakemake -p --use-conda --cores 4 -s scripts/snakemake/simu/ngs/Snakemake
```

The `Snakemake` files, the files that define a workflow, are configured in a way s.t. they find the configuration file automatically.

⚠️ Without further ado, please run the pipelines from the command line at the top folder of this project.

<details><summary>⚠️ Caution with outdated software packages </summary>
It is highly recommended to let the workflows utilize their designated conda environments (`--use-conda`) even if the required software is already available on the system. Outdated software packages might break the functionality of certain workflows (e.g. older versions of bcftools do not split multi-allelic sites correctly).
</details>

## Help:

Please visit the project [wiki](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/wiki) for more information, help and FAQs.
