# DamMapper

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.12.0-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/damid-seq/actions/workflows/main.yml/badge.svg)](https://github.com/niekwit/damid-seq/actions/workflows/main.yml)
[![Documentation Status](https://readthedocs.org/projects/damid-seq/badge/?version=latest)](https://damid-seq.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/708194033.svg)](https://zenodo.org/doi/10.5281/zenodo.10737672)

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## Aim

`DamMapper` is a containerized Snakemake pipeline for reproducible and scalable analysis of single/paired-end DamID short read Illumina data.


## Documentation

Documentation of how to use `DamMapper` can be found at https://damid-seq.readthedocs.io/en/latest/.

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
    - Archive of the analysis
* Output files
    - Quality control
    - Visualization of damid data
    - Peak-related plots
* Citation
    - Workflow
    - Software used in workflow
* literature
    - DamID
* Source code
    - Reporting issues