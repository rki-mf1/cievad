[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A520.04.0-23aa62.svg)](https://www.nextflow.io/)

# CIEVaD
<ins>C</ins>ontinuous <ins>I</ins>ntegration and <ins>E</ins>valuation for <ins>Va</ins>riant <ins>D</ins>etection. This repository provides a tool suite for simple, streamlined and rapid creation and evaluation of genomic variant callsets. It is primarily designed for continuous integration of variant detection software and a plain containment check between sets of variants. The tools suite utilizes the _conda_ package management system and _nextflow_ workflow language.

## Contents:
1. [System requirements](#system-requirements)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Help](#help)


## System requirements:

This tool suite was developed for Linux and is the only officially supported operating system here.
Having any derivative of the conda package management system installed is the only strict system requirement.
A recent version (≥20.04.0) of nextflow is required to execute the workflows, but can easily be installed via conda.
For an installation instruction of nextflow via conda see [Installation](#installation).

<details><summary>🖥️ See list of tested setups: </summary>
   
| Requirement | Tested with |
| --- | --- |
| 64 bits Linux operating system | Ubuntu 20.04.5 LTS |
| [Conda](https://docs.conda.io/en/latest/) | vers. 23.5.0, 24.1.2|
| [Nextflow](https://nextflow.io/) | vers. 20.04.0, 23.10.1 |

</details>


## Installation:

1. Download the repository:
```
git clone https://github.com/rki-mf1/cievad.git
```

2. [Optional] Install nextflow if not yet on your system. For good practise you should use a new conda environment:
```
conda deactivate
conda create -n cievad -c bioconda nextflow
conda activate cievad
```


## Usage:

This tool suite provides multiple functional features to generate synthetic sequencing data, generate sets of ground truth variants (truthsets) and evaluate sets of predicted variants (callsets).
There are two main workflows, `hap.nf` and `eval.nf`. 
Both workflows are executed via the nextflow command line interface (CLI).
<details><summary>⚠️ Run commands from the root directory: </summary>
Without further ado, please run the commands from a terminal at the top folder (root directory) of this repository.
Otherwise relative paths within the workflows might be invalid.
</details>

### Generating haplotype data
The minimal command to generate haplotype data is
```
nextflow run hap.nf -profile local,conda
```

This generates the following data within the `<project_root>/results/` directory:
- a haplotype (FASTA), which is a copy of the provided reference sequence but deviates by a set of synthetic genomic variants
- the variant set (VCF) of synthetic genomic variants in the haplotype
- a set of reads (FASTQ) representing a sequencing experiment from the haplotype

### Evaluating variant calls
The minimal command to evaluate the accordance between a truthset (generated data) and a callset is
```
nextflow run eval.nf -profile local,conda --callsets_dir <path/to/callsets>
```
where `--callsets_dir` is the parameter to specify a folder containing the callset VCF files.
Currently, a callset within this folder has to follow the naming convention `callset_<X>.vcf[.gz]` where _\<X\>_ is the integer of the corresponding truthset.
Alternatively, one can provide a sample sheet of comma separated values (CSV file) with the columns "index", "truthset" and callset", where "index" is an integer from 1 to n (number of samples) and "callset"/"truthset" are paths to the pairwise matching VCF files.
Callsets can optionally be _gzip_ compressed.
The command for the sample sheet input is
```
nextflow run eval.nf -profile local,conda --sample_sheet <path/to/sample_sheet>
```

This generates the following data within the `<project_root>/results/` directory:
- a report (CSV, JSON) about accordance between the synthetic variant set and a given corresponding callset
- a report (CSV) with statistis across all tested individuals

### Tuning the workflow parameters
Many internal settings can be adjusted at the nextflow level.
The parameters to adjust the workflows are listed on their respective help pages.
To inspect the help pages type `--help` after the script name.
Parameters can be adjusted via the CLI or within the _nextflow.config_ file.
Mind that parameters provided by the CLI will overwrite parameters set in config.

## Help:

Visit the project [wiki](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/wiki) for more information, help and FAQs. <br>
Please file issues, bug reports and questions to the [issues](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/issues) section.
