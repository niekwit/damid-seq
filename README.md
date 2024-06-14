# Snakemake workflow: `damid-seq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.12.0-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/damid-seq/actions/workflows/main.yml/badge.svg)](https://github.com/niekwit/damid-seq/actions/workflows/main.yml)
[![CodeFactor](https://www.codefactor.io/repository/github/niekwit/damid-seq/badge)](https://www.codefactor.io/repository/github/niekwit/damid-seq)
[![DOI](https://zenodo.org/badge/708194033.svg)](https://zenodo.org/doi/10.5281/zenodo.10737672)

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## Aim

`damid-seq` is a containerized Snakemake pipeline for reproducible analysis of single/paired-end DamID-seq short read Illumina data.

The core of the pipeline is the Perl script [damidseq_pipeline](https://github.com/owenjm/damidseq_pipeline), which is a great tool for the first steps of analysing DamID-seq data. However, it does not process biological replicate data, and is not written with deployment to server, cluster, grid and cloud environments in mind.

`damid-seq` implements the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system, which overcomes the above issues. In addition, we have added many features to the DamID-seq analysis workflow.

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