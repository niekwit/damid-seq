FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="bc495fc8dab6622dc357e220ee04da82e23bb5f3f8e7a47eab19bc8da4c36188"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/bowtie2/align/environment.yaml
#   prefix: /conda-envs/39f740686257be1f673b0c82e1187e65
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.4
#     - samtools =1.21
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/39f740686257be1f673b0c82e1187e65
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/bowtie2/align/environment.yaml /conda-envs/39f740686257be1f673b0c82e1187e65/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/bowtie2/build/environment.yaml
#   prefix: /conda-envs/4f56e195d46da93fbd8fc814617891ce
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.4
RUN mkdir -p /conda-envs/4f56e195d46da93fbd8fc814617891ce
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/bowtie2/build/environment.yaml /conda-envs/4f56e195d46da93fbd8fc814617891ce/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/fastqc/environment.yaml
#   prefix: /conda-envs/24b8923f8e4abe077ffe95b01bfc1652
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/24b8923f8e4abe077ffe95b01bfc1652
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/fastqc/environment.yaml /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/multiqc/environment.yaml
#   prefix: /conda-envs/75e53aec85d440a0b63166c491b3fc09
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.27.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/75e53aec85d440a0b63166c491b3fc09
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/multiqc/environment.yaml /conda-envs/75e53aec85d440a0b63166c491b3fc09/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/samtools/faidx/environment.yaml
#   prefix: /conda-envs/aca66edd82a97336414a42c80d339d33
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.21
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/aca66edd82a97336414a42c80d339d33
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/samtools/faidx/environment.yaml /conda-envs/aca66edd82a97336414a42c80d339d33/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/samtools/index/environment.yaml
#   prefix: /conda-envs/85a51d707f59b4dd6fc7102fa747fe0b
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.21
RUN mkdir -p /conda-envs/85a51d707f59b4dd6fc7102fa747fe0b
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/samtools/index/environment.yaml /conda-envs/85a51d707f59b4dd6fc7102fa747fe0b/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/trim_galore/se/environment.yaml
#   prefix: /conda-envs/975a4a78cc02737b050b54b2645b9ecc
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - trim-galore =0.6.10
RUN mkdir -p /conda-envs/975a4a78cc02737b050b54b2645b9ecc
ADD https://github.com/snakemake/snakemake-wrappers/raw/v5.8.3/bio/trim_galore/se/environment.yaml /conda-envs/975a4a78cc02737b050b54b2645b9ecc/environment.yaml

# Conda environment:
#   source: workflow/envs/R.yaml
#   prefix: /conda-envs/4554a2c6669fb2054cdfbe43af1d6154
#   name: R
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-base=4.3
#     - r-tidyverse=2.0.0
#     - r-cowplot=1.1.2
#     - r-reshape2=1.4.4
#     - r-rcolorbrewer=1.1_3
#     - r-ggrepel=0.9.4
#     - r-scales=1.3.0
#     - bioconductor-chipseeker=1.38.0
#     - bioconductor-genomicfeatures=1.54.1
#     - bioconductor-rtracklayer=1.62.0
#     - r-viridislite=0.4.2
#     - r-enrichr=3.2
#     - r-openxlsx=4.2.5.2
RUN mkdir -p /conda-envs/4554a2c6669fb2054cdfbe43af1d6154
COPY workflow/envs/R.yaml /conda-envs/4554a2c6669fb2054cdfbe43af1d6154/environment.yaml

# Conda environment:
#   source: workflow/envs/damid.yaml
#   prefix: /conda-envs/cb2a9c9685efe98f401d666d5b884a76
#   name: damid
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - pigz=2.8
#     - samtools=1.19
#     - bowtie2=2.5.3
#     - perl=5.32.1
#     - perl-inline-c=0.81
#     - r-base=4.3.2
#     - igvtools=2.16.2
#     - pandas=2.2.0
#     - pyarrow=14.0.2 #pandas dependency, somehow not in pandas recipe
RUN mkdir -p /conda-envs/cb2a9c9685efe98f401d666d5b884a76
COPY workflow/envs/damid.yaml /conda-envs/cb2a9c9685efe98f401d666d5b884a76/environment.yaml

# Conda environment:
#   source: workflow/envs/deeptools.yaml
#   prefix: /conda-envs/98f58032ba6ed531ea796d1ae05c189a
#   name: deeptools
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - deeptools=3.5.4
#     - python=3.10
#     - ucsc-bedgraphtobigwig=445
#     - wiggletools=1.2.11
#     - ucsc-wigtobigwig=447
#     - pysam=0.21.0
#     - gffpandas=1.2.2
#     - numpy=1.26.4
#     - samtools=1.13
RUN mkdir -p /conda-envs/98f58032ba6ed531ea796d1ae05c189a
COPY workflow/envs/deeptools.yaml /conda-envs/98f58032ba6ed531ea796d1ae05c189a/environment.yaml

# Conda environment:
#   source: workflow/envs/peak_calling.yaml
#   prefix: /conda-envs/bda24173dea9eee42d0a3ee88a8a19fe
#   name: peak_calling
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - python=3.9
#     - pandas=2.1.4
#     - bedtools=2.31.1
#     - bedops=2.4.41
#     - biopython=1.78
#     - bioconda::macs3=3.0.2
RUN mkdir -p /conda-envs/bda24173dea9eee42d0a3ee88a8a19fe
COPY workflow/envs/peak_calling.yaml /conda-envs/bda24173dea9eee42d0a3ee88a8a19fe/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/39f740686257be1f673b0c82e1187e65 --file /conda-envs/39f740686257be1f673b0c82e1187e65/environment.yaml && \
    conda env create --prefix /conda-envs/4f56e195d46da93fbd8fc814617891ce --file /conda-envs/4f56e195d46da93fbd8fc814617891ce/environment.yaml && \
    conda env create --prefix /conda-envs/24b8923f8e4abe077ffe95b01bfc1652 --file /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml && \
    conda env create --prefix /conda-envs/75e53aec85d440a0b63166c491b3fc09 --file /conda-envs/75e53aec85d440a0b63166c491b3fc09/environment.yaml && \
    conda env create --prefix /conda-envs/aca66edd82a97336414a42c80d339d33 --file /conda-envs/aca66edd82a97336414a42c80d339d33/environment.yaml && \
    conda env create --prefix /conda-envs/85a51d707f59b4dd6fc7102fa747fe0b --file /conda-envs/85a51d707f59b4dd6fc7102fa747fe0b/environment.yaml && \
    conda env create --prefix /conda-envs/975a4a78cc02737b050b54b2645b9ecc --file /conda-envs/975a4a78cc02737b050b54b2645b9ecc/environment.yaml && \
    conda env create --prefix /conda-envs/4554a2c6669fb2054cdfbe43af1d6154 --file /conda-envs/4554a2c6669fb2054cdfbe43af1d6154/environment.yaml && \
    conda env create --prefix /conda-envs/cb2a9c9685efe98f401d666d5b884a76 --file /conda-envs/cb2a9c9685efe98f401d666d5b884a76/environment.yaml && \
    conda env create --prefix /conda-envs/98f58032ba6ed531ea796d1ae05c189a --file /conda-envs/98f58032ba6ed531ea796d1ae05c189a/environment.yaml && \
    conda env create --prefix /conda-envs/bda24173dea9eee42d0a3ee88a8a19fe --file /conda-envs/bda24173dea9eee42d0a3ee88a8a19fe/environment.yaml && \
    conda clean --all -y