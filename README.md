![Static Badge](https://img.shields.io/badge/requires-conda-blue)
![Static Badge](https://img.shields.io/badge/requires-nextflow-blue)

# CIEVaD
<ins>C</ins>ontinuous <ins>I</ins>ntegration and <ins>E</ins>valuation for <ins>Va</ins>riant <ins>D</ins>etection. This repository provides a tool suite for simple, streamlined and rapid creation and evaluation of genomic variant callsets. It is primarily designed for continuous integration of variant detection software and a plain containment check between sets of variants. The tools suite utilizes the _conda_ package management system and _nextflow_ workflow language.

## Contents:
1. [System requirements](#system-requirements)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Output](#output)
5. [Help](#help)


## System requirements:

This tool suite was developed for Linux and is the only officially supported operating system here.
Having any derivative of the conda package management system installed is the only strict system requirement.
A recent version (≥20.04.0) of nextflow is required to execute the workflows, but can easily be installed via conda.
For an installation instruction of nextflow via conda see [Installation](#installation).

<details><summary>🛠️ See list of tested setups: </summary>
   
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
The current list and roadmap of principal functionality is:
* [x] Generating synthetic haplotypes from a given reference genome. This returns a haplotype sequence (FASTA) and its set of variants (VCF) with respect to the reference.
* [x] Generating synthetic NGS reads from a given haplotype
* [ ] Generating synthetic amplicon sequences from a given reference genome and generating synthetic reads from those amplicons
* [x] Generating synthetic long-reads from a given haplotype
* [x] Evaluate compliance between sets of variants

### Generating haplotype data
The minimal command to generate haplotype data is
```
nextflow run hap.py -profile local,conda
```

### Evaluating variant calls
The minimal command to evaluate the accordance between a truthset (generated data) and a callset is
```
nextflow run eval.nf -profile local,conda --callsets_dir <path/to/callsets>
```
where `--callsets_dir` is the parameter to specify a folder containing the callset VCF files.
Currently, a callset within this folder has to follow the naming convention `callset_<X>.vcf[.gz]` where _\<X\>_ is the integer of the corresponding truthset.
Callsets can optionally be _gzip_ compressed.

🚧 For convenience, the `eval.nf` will get an option to provide a sample sheet as an alternative input format in the future.

<details><summary>⚠️ Run commands from the root directory </summary>
Without further ado, please run the commands from a terminal at the top folder (root directory) of this repository.
Otherwise relative paths within the workflows might be invalid.
</details>

### Tuning the workflows via CLI parameters
Many internal command line arguments can be adjusted at the nextflow level.
For instance, if long-reads in the fashion of Oxford Nanopore Technologies should be generated instead of NGS reads (default), use the `--read_tpye` parameter from the workflow's nextflow CLI:
```
nextflow run hap.nf -profile local,conda --read_type ont
```

### Tuning the workflows via the config file
\<TODO\>

## Output
\<TODO\>

## Help:

Visit the project [wiki](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/wiki) for more information, help and FAQs. <br>
Please file issues, bug reports and questions to the [issues](https://github.com/rki-mf1/imsmp-variant-calling-benchmark/issues) section.
