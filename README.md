# Snakemake workflow: `damid-seq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.12.0-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/damid-seq/actions/workflows/main.yml/badge.svg)](https://github.com/niekwit/damid-seq/actions/workflows/main.yml)
[![Documentation Status](https://readthedocs.org/projects/damid-seq/badge/?version=latest)](https://damid-seq.readthedocs.io/en/latest/?badge=latest)
[![CodeFactor](https://www.codefactor.io/repository/github/niekwit/damid-seq/badge)](https://www.codefactor.io/repository/github/niekwit/damid-seq)
[![DOI](https://zenodo.org/badge/708194033.svg)](https://zenodo.org/doi/10.5281/zenodo.10737672)

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## Aim

`damid-seq` is a containerized Snakemake pipeline for reproducible and scalable analysis of single/paired-end DamID-seq short read Illumina data.


## Documentation

Documentation of how to use `damid-seq` can be found at https://damid-seq.readthedocs.io/en/latest/.

### Overview of documentation

* DamID
    - DamID principle
    - Experimental considerations
* Software requirements
    - Conda/Mamba
    - Snakemake
    - Apptainer
    - Snakefetch
* Usage
    - Preparing raw sequencing data
    - Sample meta data and analysis settings
    - Configuration of Snakemake
    - Dry-run of the analysis
    - Visualization of the workflow
    - Running the analysis
    - Report of the results
* Output
    - Quality control
    - Visualization of damid-seq data
    - Peaks
* Citation
    - Workflow
    - Software used in workflow