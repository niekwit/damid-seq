# Snakemake workflow: `damid-seq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.25.0-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/damid-seq/actions/workflows/main.yml/badge.svg)](https://github.com/niekwit/damid-seq/actions/workflows/main.yml)
[![CodeFactor](https://www.codefactor.io/repository/github/niekwit/damid-seq/badge)](https://www.codefactor.io/repository/github/niekwit/damid-seq)
[![DOI](https://zenodo.org/badge/708194033.svg)](https://zenodo.org/doi/10.5281/zenodo.10737672)

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

# CONTENTS

* [Aim](https://github.com/niekwit/damid-seq?tab=readme-ov-file#aim)
* [DamID](https://github.com/niekwit/damid-seq?tab=readme-ov-file#damid)
* [Experimental considerations](https://github.com/niekwit/damid-seq?tab=readme-ov-file#experimental-considerations)
* [Requirements](https://github.com/niekwit/damid-seq?tab=readme-ov-file#requirements)
* [Dependency graph of Snakemake rules](https://github.com/niekwit/damid-seq?tab=readme-ov-file#dependency-graph-of-snakemake-rules)
* [Installation of Conda/Mamba](https://github.com/niekwit/damid-seq?tab=readme-ov-file#installation-of-condamamba)
* [Installation of Snakemake](https://github.com/niekwit/damid-seq?tab=readme-ov-file#installation-of-snakemake)
* [Cloning `damid-seq` GitHub repository](https://github.com/niekwit/damid-seq?tab=readme-ov-file#cloning-github-repository)
* [Preparing your data](https://github.com/niekwit/damid-seq?tab=readme-ov-file#preparing-your-data)
* [Sample meta data and analysis settings](https://github.com/niekwit/damid-seq?tab=readme-ov-file#sample-meta-data-and-analysis-settings)
* [Configuration of Snakemake](https://github.com/niekwit/damid-seq?tab=readme-ov-file#configuration-of-snakemake)
* [Running the analysis with test data](https://github.com/niekwit/damid-seq?tab=readme-ov-file#running-the-analysis-with-test-data)
* [Dry-run of the analysis](https://github.com/niekwit/damid-seq?tab=readme-ov-file#dry-run-of-the-analysis)
* [Running the analysis](https://github.com/niekwit/damid-seq?tab=readme-ov-file#running-the-analysis)
* [Report of the results](https://github.com/niekwit/damid-seq?tab=readme-ov-file#report-of-the-results)
* [References](https://github.com/niekwit/damid-seq?tab=readme-ov-file#references)

## Aim

`damid-seq` is a Snakemake pipeline for reproducible analysis of single/paired-end DamID-seq short read Illumina data.

The core of the pipeline is the Perl script [damidseq_pipeline](https://github.com/owenjm/damidseq_pipeline), which is a great tool for the first steps of analysing DamID-seq data. However, it does process biological replicate data, and is not written with deployment to server, cluster, grid and cloud environments in mind.

`damid-seq` implements the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system, which overcomes the above issues. In addition, we have added many features to the DamID-seq analysis workflow.

The output of `damid-seq` is as follows:

1. Quality control of the adapter trimmed sequencing data using FastQC/MultiQC.

2. Bigwig files for visualisation of binding in genome browsers, such IGV.

4. PCA and correlation plots for checking consistency of biological replicates

5. Identified and annotated peaks using MACS2 and/or find_peaks.pl

6. Profile plot/heatmap to visualise binding around genomic features, such as transcription start sites, usingh deeptools

## DamID

TO DO

![DamID principle (Adapted from Van den Ameele et al. 2019 Current Opinion in Neurobiology)](/images/damid.png)
Figure adapted from Van den Ameele et al. 2019 Current Opinion in Neurobiology

## Experimental considerations

TO DO

## Requirements

`damid-seq` has been extensively tested on GNU/Linux-based operating systems, so we advice to run your analysis on for example Ubuntu or Fedora.

Hardware requirements differ for the kind of data that needs to be analysed: for the analysis of mammalian data sets, > 32GB of RAM is recommended. Much less RAM is needed for analysis for data from organisms with much smaller genomes, such as _Drosophila_.

## Dependency graph of Snakemake rules

![Dependency graph of rules](/images/rule_graph.png)

## Installation of Conda/Mamba

For reproducible analysis, `damid-seq` uses Conda environments in the Snakemake workflow.

Please follow the instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for a detailed guide to install Conda/Mamba.

## Installation of Snakemake

To install Snakemake create the following environment with `mamba`:

```shell
$ mamba create -n snakemake snakemake pandas
```

Activate the environment as follows:

```shell
$ mamba activate snakemake
```

## Cloning `damid-seq` GitHub repository

To obtain the workflow code, run the following command:

```shell
$ git clone https://github.com/niekwit/damid-seq.git -b v0.3.0 path/to/analysis/directory
```

This will clone the latest version of the workflow to a specified directory:

```
path/to/analysis/directory
├── config
│   ├── config.yaml
│   ├── README.md
│   └── samples.csv
├── images
│   └── rule_graph.png
├── LICENSE
├── README.md
└── workflow
    ├── envs
    │   ├── damid.yaml
    │   ├── deeptools.yaml
    │   ├── peak_calling.yaml
    │   ├── R.yaml
    │   └── trim.yaml
    ├── report
    │   ├── annotated_peaks.rst
    │   ├── correlation.rst
    │   ├── distance_to_tss.rst
    │   ├── feature_distributions.rst
    │   ├── heatmap.rst
    │   ├── mapping_rates.rst
    │   ├── pca.rst
    │   ├── profile_plot.rst
    │   ├── scree.rst
    │   └── workflow.rst
    ├── rules
    │   ├── bedgraph_processing.smk
    │   ├── bed.smk
    │   ├── damid.smk
    │   ├── deeptools.smk
    │   ├── fastqc.smk
    │   ├── motifs.smk
    │   ├── peak_calling.smk
    │   ├── plotting.smk
    │   ├── resources.smk
    │   └── trimming.smk
    ├── scripts
    │   ├── annotate_peaks.R
    │   ├── average_bigwig.py
    │   ├── average_wig.py
    │   ├── convert_bed2fasta.py
    │   ├── create_annotation_file.R
    │   ├── create_background_fasta.py
    │   ├── create_blacklist.py
    │   ├── damidseq_pipeline.py
    │   ├── filter_overlapping_peaks.py
    │   ├── general_functions.smk
    │   ├── get_resource.sh
    │   ├── install_homer_genome.py
    │   ├── peak_annotation_plots.R
    │   ├── plot_mapping_rates.R
    │   ├── plot_PCA.R
    │   ├── resources.py
    │   ├── run_find_peaks.py
    │   └── trim_galore.py
    └── Snakefile
```

## Preparing your data

Place your raw sequencing in path/to/analysis/directory/reads:

Each biological replicate should be placed in a seperate directory as follows:

```
path/to/analysis/directory/reads
├── reads
│   ├── exp1
│   │   ├── Dam.fastq.gz
│   │   ├── HIF1A.fastq.gz
│   │   └── HIF2A.fastq.gz
│   ├── exp2
│   │   ├── Dam.fastq.gz
│   │   ├── HIF1A.fastq.gz
│   │   └── HIF2A.fastq.gz
│   └── exp3
│       ├── Dam.fastq.gz
│       ├── HIF1A.fastq.gz
│       └── HIF2A.fastq.gz
```

The above example is for single-end (SE) data, where the file extension should be .fastq.gz.

Paired-end (PE) data should end with \_R1\_001.fastq.gz/\_R2\_001.fastq.gz for read 1 and read 2, respectively.

The Dam only control should always be Dam.fastq.gz (SE) or Dam.\_R1\_001.fastq.gz/Dam.\_R2\_001.fastq.gz

## Sample meta data and analysis settings

The `config` directory contains `samples.csv` with sample meta data as follows:

| sample    | genotype | treatment |
|-----------|----------|-----------|
|HIF1A      | WT       | Hypoxia   |
|HIF2A      | WT       | Hypoxia   |
|Dam        | WT       | Hypoxia   | 

`config.yaml` in the same directory contains the settings for the analysis:

```yaml
genome: dm6
ensembl_genome_build: 110
plasmid_fasta: none
fusion_genes: FBgn0038542,FBgn0085506 # Genes from these proteins will be removed from the analysis
damidseq_pipeline:
  binsize: 300
  extra: "" # extra argument for damidseq_pipeline
deeptools:
  bamCoverage: # bam to bigwig conversion for QC
    binSize: 10
    normalizeUsing: RPKM
    extra: ""
  matrix: # Settings for computeMatrix
    mode: scale-regions # scale-regions or reference-point
    referencePoint: TSS # TSS, TES, center (only for reference-point mode)
    regionBodyLength: 6000
    upstream: 3000
    downstream: 3000
    binSize: 10
    averageTypeBins: mean
    regionsFileName: "" # BED or GTF file(s) with regions of interest (optional, whole genome if not specified)
    no_whole_genome: False # If True, will omit whole genome as region and only use regionsFileName(s)
    extra: "" # Any additional parameters for computeMatrix
  plotHeatmap:
    interpolationMethod: auto
    plotType: lines # lines, fill, se, std
    colorMap: viridis # https://matplotlib.org/2.0.2/users/colormaps.html
    alpha: 1.0
    extra: "" 
peak_calling_perl:
  run: True
  iterations: 5 # N argument
  fdr: 0.01
  fraction: 0 # Fraction of random fragments to consider per iteration
  min_count: 2 # Minimum number of reads to consider a peak
  min_quantile: 0.95 # Minimum quantile for considering peaks
  step: 0.01 # Stepping for quantiles
  unified_peaks: max # Method for calling peak overlaps. 'min': call minimum overlapping peak area. 'max': call maximum overlap as peak
  extra: ""
  overlapping_peaks:
    max_size: 10 # Maximum size of peaks to be extended
    extend_by: 40 # Number of bp to extend peaks on either side
    keep: 2 # Minimum number peaks that must overlap to keep
peak_calling_macs2:
  run: False
  mode: narrow
  qvalue: 0.05
  broad_cutoff: 0.1
  extra: ""
resources: # computing resources
  trim:
    cpu: 8
    time: 60
  fastqc:
    cpu: 4
    time: 60
  damid:
    cpu: 24
    time: 720
  index:
    cpu: 36
    time: 60
  deeptools:
    cpu: 8
    time: 90
  plotting:
    cpu: 2
    time: 20
```

## Configuration of Snakemake

Running Snakemake can entail quite a few command line flags. To make this easier these can be set in a global profile that is defined in a user-specific configuration directory in order to simplify this process.

For example, a profile `config.yaml` can be stored at /home/user/.config/snakemake/profile:
```yaml
printshellcmds: True
cache: True
use-conda: True
cores: 32
rerun-incomplete: True
show-failed-logs: True
```

Snakemake supports between workflow caching, so that certain resource files, such as the Bowtie2 index, can be re-used between different analyses.

To enable this append this line to your `~/.bashrc`:
```shell
export SNAKEMAKE_OUTPUT_CACHE=/path/to/snakemake-cache/
```

## Running the analysis with test data

The .test directory in the cloned GitHub repository contains a small data set that allows for a test run of the workflow:

```shell
$ snakemake --profile /home/user/.config/snakemake/profile --directory .test/
```

## Dry-run of the analysis

Before running the actual analyis with your own data, a dry-run can be performed:

```shell
$ cd path/to/analysis/directory
$ snakemake -np
```

Snakemake will create the DAG of jobs, but it will not execute anything.

## Running the analysis

Once you know that the test and/or dry run has worked, the actual analysis can be initiated as follows:
```shell
$ snakemake --profile /home/user/.config/snakemake/profile --directory .test/
```

> [!IMPORTANT]  
> Always make sure to use the absolute path (i.e. /home/user/.config/...) rather than the relative path (~/.config/...) when providing the path for the profile file.

## Report of the results

When the analysis has finished succesfully, an HTML report can be created as follows:

```shell
$ snakemake --report report.html
```

This report will contain run time information for the Snakemake rules, as well as figures generated by the workflow, and the code used to create these.

## References

TO DO