# Snakemake workflow: `damid-seq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.25.0-brightgreen.svg)](https://snakemake.github.io)
[![Tests](https://github.com/niekwit/damid-seq/actions/workflows/main.yml/badge.svg)](https://github.com/niekwit/damid-seq/actions/workflows/main.yml)
[![DOI](https://zenodo.org/badge/708194033.svg)](https://zenodo.org/doi/10.5281/zenodo.10737672)

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

The usage of this workflow is briefly described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=niekwit%2Fdamid-seq). For a more detailed description see below.

# Aim

Snakemake pipeline for reproducible analysis or single/paired-end DamID-seq short read Illumina data.

The output of this pipeline is as follows:

1. Quality control of the raw sequencing data using FastQC/MultiQC.

2. Bigwig files for visualisation of binding in genome browsers, such IGV.

3. Identified and annotated peaks using MACS2

4. Profile plot/heatmap to visualise binding around genomic features, such as transcription start sites, usingh deeptools


# Rule graph

![Rule graph of the Snakemake damid-seq pipeline.](/images/rule_graph.png)

# Usage

## Directory structure

The analysis directory should have the following stucture

```
.
├── config
│   ├── config.yaml
│   └── samples.csv
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
├── resources
└── workflow
    ├── envs
    ├── rules
    ├── scripts
    └── Snakefile
```

In the above example, the reads directory contains a separate subdirectory for each experimental replicate, as `damidseq_pipeline` can only handle one replicate at a time. The files names (without extension) should match the sample names in `config/samples.csv` (see below).

Raw sequencing data can either be single-end (ending with .fastq.gz) or paired-end (ending with _R1_001.fastq.gz/_R2_001.fastq.gz)

## Sample meta data and analysis settings

In the `config` directory include `samples.csv` that contains the sample meta data as follows:

| sample    | genotype | treatment |
|-----------|----------|-----------|
|HIF1A      | WT       | Hypoxia   |
|HIF2A      | WT       | Hypoxia   |
|Dam        | WT       | Hypoxia   | 


`config.yaml` in the same directory contains the settings for the analysis:

```
genome: "hg38"
ensembl_genome_build: "110"
extra: "" # extra argument for damidseq_pipeline
fusion_genes: HIF1A,EPAS1 # Genes from these proteins will be removed from the analysis
deeptools:
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
peak_calling:
  iterations: 100 # N argument
  fdr: 0.01
  extra: ""
  overlapping_peaks:
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
    cpu: 8
    time: 120
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




