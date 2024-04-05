![Static Badge](https://img.shields.io/badge/requires-conda-blue)
![Static Badge](https://img.shields.io/badge/requires-nextflow-blue)

# CIEVaD
<ins>C</ins>ontinuous <ins>I</ins>ntegration and <ins>E</ins>valuation for <ins>Va</ins>riant <ins>D</ins>etection. This repository provides a tool suite for simple, streamlined and rapid creation and evaluation of genomic variant callsets. It is primarily designed for continuous integration of variant detection software and a plain containment check between sets of variants. The tools suite utilizes the _conda_ package management system and _nextflow_ workflow language.

## Contents:
1. [System requirements](#system-requirements)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Help](#help)


## System requirements:

This tool suite was developed under Linux/UNIX and is the only officially supported operating system here.
Having any derivative of the `conda` package management system installed is the only strict system requirement.
Having a recent `nextflow` (≥20.04.0) version installed is required to execute the workflows.
For an installation instruction of `nextflow` via conda see [Installation](#installation).

<details><summary>🛠️ See tested setups: </summary>
   
| Requirement | Tested with |
| --- | --- |
| 64 bits operating system | Ubuntu 20.04.5 LTS |
| [Conda](https://docs.conda.io/en/latest/) | vers. 23.5.0, 24.1.2|
| [Nextflow](https://nextflow.io/) | vers. 20.04.0 |

</details>


## Installation:

1. Download the repository:
```
git clone https://github.com/rki-mf1/imsmp-variant-calling-benchmark.git
```

2. [Optional] Install nextflow if not yet on your system. For good practise you should use a new conda environment:
```
conda deactivate
conda create -n cievad -c bioconda nextflow
conda activate cievad
```


## Usage:

This tool suite provides multiple workflows to generate synthetic sequencing data and evaluate sets of predicted variants (callsets).
A full list of workflows, their respective modules in the python command line interface (CLI) and a detailed description of input and output files can be found in this [wiki](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/wiki) page of the repository.
The current list of principal functionality is:
* Generating synthetic haplotypes from a given reference genome
* Generating synthetic NGS reads from a given haplotype
* Generating synthetic amplicon sequences from a given reference genome and generating synthetic NGS reads from the amplicons
* Generating synthetic long-reads from a given haplotype
* Evaluate compliance between sets of variants

The repository provides a simple CLI for a convenient application-like user experience with the underlying Snakemake workflows.
The CLI is started from the root directory via
```
python cievad.py --help
```
and each individual module provides another help page via its sub-command
```
python cievad.py <module> --help
```

```
nextflow run hap.nf -profile local,conda
```

<details><summary>⚠️ Run commands from the root directory </summary>
Without further ado, please run the commands from a terminal at the top folder (root directory) of this repository.
Otherwise relative paths within the workflows might be invalid.
</details>


## Help:

Visit the project [wiki](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/wiki) for more information, help and FAQs. <br>
Please file issues, bug reports and questions to the [issues](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/issues) section.
