# imsmp-variant-calling-benchmark
![Static Badge](https://img.shields.io/badge/conda-required-blue)
![Static Badge](https://img.shields.io/badge/python-≥3.2-blue)
![Static Badge](https://img.shields.io/badge/snakemake-≥6.0.0-blue)

This repository provides a tool suite for a simple, streamlined and rapid evaluation of variant callsets (SNPs/indels). Without loss of generality, the modules were implemented for the purpose of benchmarking variant callers, evaluating genomic data workflows or continuous integration of the aforementioned. The underlying workflows utilize the conda package management system and the Snakemake workflow language.

## Contents:
1. [System requirements](#system-requirements)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Help](#help)


## System requirements:

This tools suit was developed under Linux/UNIX but might work on other operating systems.
However, only 64 bits POSIX-compliant operating systems will be officially supported here.
Having any derivative of the `conda` package management system installed is the only strict on system requirement.
Having a recent `Snakemake` version installed will be required but the next section also provides instructions on how Snakemake can easily be installed via conda.

<details><summary> See tested setups: </summary>
   
| Requirement | Tested with |
| --- | --- |
| 64 bits POSIX-compliant operating system | Ubuntu 20.04.5 LTS |
| [conda](https://docs.conda.io/en/latest/) | vers. 23.5.0 |
| [snakemake](https://snakemake.readthedocs.io/en/stable/) | vers. 7.25.3 |

</details>

## Installation:

1. Download the repository:
```
git clone https://github.com/rki-mf1/imsmp-variant-calling-benchmark.git
```

2. [Optional] Install Snakemake if not yet on your system. You can use the conda environment description file provided in this repository:
```
conda env create -f env/conda_snakemake7.yaml
conda activate snakemake7
```

## Usage:

Running any of the workflows consists of two main steps:
   1. creating a configuration _yaml_ file
   2. running a snakemake workflow

### Config files
Each workflow is defined by a _Snakefile_ which again requires a configuration file in _yaml_ format.
The designated location within this repository and all parameters of a workflow that have to be set within the configuration file are listed in the Snakefile's `CONFIG` and `PARAMS` section, respectively.
This repository's [wiki page](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/wiki/Confiluration-files) contains more explanation, details and examples of configuration files. <br>

An example configuration file `snake_config.yaml` could look like the following:

```
HEAD_DIR:
   <path to this cloned repository>
REF:
   <path to my reference genome>
SAMPLES:
   - "001"
   - "002"
   - "003"
```

:construction: (TODO) There will likely be a simplification, e.g. a python module that generates the configs in the right location 

### Run snakemake
With a configuration file in place running a Snakemake workflow could look like the following:
```
snakemake -p --use-conda --cores 4 -s scripts/snakemake/simu/ngs/Snakefile
```

The `Snakefile`, the files that define a workflow, are configured in a way s.t. they find the configuration files automatically.

<details><summary>⚠️ Please run Snakemake from the root directory </summary>
Without further ado, please run the Snakemake workflows from a terminal at the top folder (root directory) of this project.
Otherwise relative paths within the workflows might be invalid.
</details>

<details><summary>⚠️ Caution with outdated software packages </summary>
It is highly recommended to let the workflows utilize their designated conda environments (--use-conda) even if the required software is already available on the system.
Outdated software packages might break the functionality of certain workflows (e.g. older versions of bcftools do not split multi-allelic sites correctly).
</details>

## Help:

:mag: Please visit the project [wiki](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/wiki) for more information, help and FAQs. <br>
:hammer: Please file issues, bug reports and questions to the [issues](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/issues) section.
