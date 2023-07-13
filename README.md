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

This tool suite provides multiple features to simulate training data and evaluate callsets. 
Running any of the features consists of two principal steps:
   1. Generating a configuration file
   2. Running a Snakemake workflow

A full list of features, their respective modules in the python command line interface (CLI) and a description of output files can be found in the [wiki](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/wiki) of this repository.
The following two subsections are showcasing one particular feature of the tool suite.

### Configuration files
Each functional feature of the tool suite requires the user to generate a configuration file in _yaml_ format.
For instance, the _haplotype simulation_ feature requires a configuration file with designated parameters for the haplotype simulation workflow in the second step.
The python CLI facilitates generating such a configuration file.
For instance, in order to generate a configuration file for the haplotype simulation (with default parameters) run:
```
python scripts/pythonUI/vc_benchmark.py config-hap-simu </path/to/this/repository/> </path/to/a/reference/genome.fasta>
```
All modules that generate configuration files (starting with "config-") can be found at the help page of the python CLI:
   
```
python scripts/pythonUI/vc_benchmark.py --help
```

### Snakemake workflows
The second step of each functional feature is running a Snakemake workflow.
For instance, the _haplotype simulation_ feature uses the configuration file generated in the first step to run its corresponding Snakemake workflow.
The python CLI facilitates running such a Snakemake workflow.
For instance, in order to run the Snakemake workflow for the haplotype simulation (with default parameters) run:
```
python scripts/pythonUI/vc_benchmark.py run-hap-simu
```
All modules that run a Snakemake workflow (starting with "run-") can be found at the help page of the python CLI:
```
python scripts/pythonUI/vc_benchmark.py --help
```

<details><summary>⚠️ Run commands from the root directory </summary>
Without further ado, please run the commands from a terminal at the top folder (root directory) of this repository.
Otherwise relative paths within the workflows might be invalid.
</details>


## Help:

Visit the project [wiki](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/wiki) for more information, help and FAQs. <br>
Please file issues, bug reports and questions to the [issues](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/issues) section.
