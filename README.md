# Snakemake workflow: `damid-seq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.10.6-brightgreen.svg)](https://snakemake.github.io)
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
* [Literature](https://github.com/niekwit/damid-seq?tab=readme-ov-file#references)

## Aim

`damid-seq` is a Snakemake pipeline for reproducible analysis of single/paired-end DamID-seq short read Illumina data.

The core of the pipeline is the Perl script [damidseq_pipeline](https://github.com/owenjm/damidseq_pipeline), which is a great tool for the first steps of analysing DamID-seq data. However, it does not process biological replicate data, and is not written with deployment to server, cluster, grid and cloud environments in mind.

`damid-seq` implements the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system, which overcomes the above issues. In addition, we have added many features to the DamID-seq analysis workflow.

The output of `damid-seq` is as follows:

1. Quality control of the adapter trimmed sequencing data using FastQC/MultiQC.

2. Bigwig files for visualisation of binding in genome browsers, such IGV.

4. PCA and correlation plots for checking consistency of biological replicates

5. Identified and annotated peaks using MACS2 and/or find_peaks.pl

6. Profile plot/heatmap to visualise binding around genomic features, such as transcription start sites, usingh deeptools

## DamID

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
$ mamba create -n snakemake snakemake
```

Activate the environment as follows:

```shell
$ mamba activate snakemake
```

If you want to deploy Snakemake on an HPC system using slurm also run:

```shell
$ pip install snakemake-executor-plugin-slurm
```

## Cloning `damid-seq` GitHub repository

The easiest way to obtain the workflow code is to use [snakefetch](https://pypi.org/project/snakefetch/):

```shell
$ pip install snakefetch
$ snakefetch --outdir /path/to/analysis --repo-version v0.4.0 --url https://github.com/niekwit/damid-seq
Downloading archive file for version v0.4.0 from https://github.com/niekwit/damid-seq...
Extracting config and workflow directories from tar.gz file to /home/niek/Downloads/TEST...
Done!
```

This will copy the config and workflow directories to the path set with the `--outdir` flag.

## Preparing raw sequencing data

In the directory containing config/workflow create a directory called reads:

```shell
$ cd /path/to/analysis
$ mdkir -p reads 
```

Data files from each group of biological replicates should be placed into a unique folder, e.g.:

```shell
reads
├── exp1
│   ├── Dam.fastq.gz
│   ├── HIF1A.fastq.gz
│   └── HIF2A.fastq.gz
├── exp2
│   ├── Dam.fastq.gz
│   ├── HIF1A.fastq.gz
│   └── HIF2A.fastq.gz
└── exp3
    ├── Dam.fastq.gz
    ├── HIF1A.fastq.gz
    └── HIF2A.fastq.gz
```

> [!IMPORTANT]  
> Single-end fastq files should always end with *fastq.gz*, while paired-end reads should end with *\_R1\_001.fastq.gz/\_R2\_001.fastq.gz*

> [!IMPORTANT]  
> The Dam only control should always be called Dam.*relevant_extension*

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
bowtie2:
  extra: ""
damidseq_pipeline:
  normalization: kde # kde, rpm or rawbins
  binsize: 300
  extra: "" # extra argument for damidseq_pipeline
quantile_normalisation:
  apply: True
  extra: "" # extra arguments for quantile_normalization
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
    binSize: 100
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
peak_calling_macs2:
  run: False
  mode: narrow
  qvalue: 0.05 # for narrow peaks
  broad_cutoff: 0.1 # for broad peaks
  extra: ""
consensus_peaks:
    max_size: 10 # Maximum size of peaks to be extended
    extend_by: 40 # Number of bp to extend peaks on either side
    keep: 2 # Minimum number peaks that must overlap to keep
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
    cpu: 40
    time: 60
  deeptools:
    cpu: 8
    time: 90
  plotting:
    cpu: 2
    time: 20
```

A lot of the DamID signal can come from the plasmids that were used to express the Dam-POIs, and this can skew the analysis.

To prevent this, two approaches are available:

1.  The genes (Ensembl gene IDs) fused to Dam can be set in config.yaml["fusion_genes] (separated by commas if multiple plasmids are used). This will mask the genomic locations of these genes in the fasta file that will be used to build the Bowtie2 index, hence excluding these regions from the analysis. 

> [!NOTE]
> To disable this function set the value of config.yaml["fusion_genes] to "".

2. If a plasmid is used that for example also uses an endogenous promoter besides the Dam fusion proteins, one can set a path to a fasta file containg all the plasmid sequences in config.yaml[""]. Trimmed reads are first aligned to these sequences, and the resulting non-aligning reads will then be processed as normal.

It is recommended to store this file in a directory called resources within the analysis folder (this folder will also contain all other non-experimental files such as fasta and gtf files).

> [!NOTE]
> To disable this function set the value of config.yaml["plasmid_fasta"] to none.


## Configuration of Snakemake

Running Snakemake can entail quite a few command line flags. To make this easier these can be set in a global profile that is defined in a user-specific configuration directory in order to simplify this process.

For example, a profile `config.yaml` can be stored at /home/user/.config/snakemake/profile:
```yaml
cores: 40
latency-wait: 20
use-conda: True
keep-going: False
rerun-incomplete: True
printshellcmds: True
cache: True
show-failed-logs: True
```

Snakemake supports between workflow caching, so that certain resource files, such as the Bowtie2 index, can be re-used between different analyses.

To enable this append this line to your `~/.bashrc`:
```shell
export SNAKEMAKE_OUTPUT_CACHE=/path/to/snakemake-cache/
```

## Dry-run of the analysis

Before running the actual analyis with your own data, a dry-run can be performed:

```shell
$ cd path/to/analysis/directory
$ snakemake -np
```

Snakemake will create the DAG of jobs and print the shell command, but it will not execute anything.

## Visualization of workflow

To visualize the workflow run (this command excludes the target rule):
```shell
$ mkdir -p images
$ snakemake --forceall --rulegraph | grep -v '\-> 0\|0\[label = \"all\"' | dot -Tpng > images/rule_graph.png
```

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

## Literature

:information_source: Some key DamID papers:

Van Steensel and Henikoff. Identification of in vivo DNA targets of chromatin proteins using tethered Dam methyltransferase. 2000 Nature Biotechnology.

Marshall et al. Cell-type-specific profiling of protein–DNA interactions without cell isolation using targeted DamID with next-generation sequencing. 2016 Nature Protocols.

Van den Ameele, Krautz and Brand. TaDa! Analysing cell type-specific chromatin in vivo with Targeted DamID. 2019 Current Opinion in Neurobiology.