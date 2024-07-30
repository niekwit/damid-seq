FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="fb3d57cd2a2fb06ea2b46a135b44aba82c9228e52138e67a6ef7253784eda550"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/bowtie2/align/environment.yaml
#   prefix: /conda-envs/79f24222e503aa3928bb131a6d96239b
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.4
#     - samtools =1.20
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/79f24222e503aa3928bb131a6d96239b
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/bowtie2/align/environment.yaml /conda-envs/79f24222e503aa3928bb131a6d96239b/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/bowtie2/build/environment.yaml
#   prefix: /conda-envs/84012557432533881e14141abce4094f
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - bowtie2 =2.5.4
RUN mkdir -p /conda-envs/84012557432533881e14141abce4094f
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/bowtie2/build/environment.yaml /conda-envs/84012557432533881e14141abce4094f/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/fastqc/environment.yaml
#   prefix: /conda-envs/90ea9a0a5f4de478eb6aa43932774175
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/90ea9a0a5f4de478eb6aa43932774175
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/fastqc/environment.yaml /conda-envs/90ea9a0a5f4de478eb6aa43932774175/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/multiqc/environment.yaml
#   prefix: /conda-envs/df61192e96fe9bc3d27cde8b6a0fa36b
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.23
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/df61192e96fe9bc3d27cde8b6a0fa36b
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/multiqc/environment.yaml /conda-envs/df61192e96fe9bc3d27cde8b6a0fa36b/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/samtools/faidx/environment.yaml
#   prefix: /conda-envs/b29db9b4576aa24f6ca5cdabda9911b4
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.20
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/b29db9b4576aa24f6ca5cdabda9911b4
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/samtools/faidx/environment.yaml /conda-envs/b29db9b4576aa24f6ca5cdabda9911b4/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/samtools/index/environment.yaml
#   prefix: /conda-envs/adec7d399fa9492178d9caa96f93852d
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.20
RUN mkdir -p /conda-envs/adec7d399fa9492178d9caa96f93852d
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/samtools/index/environment.yaml /conda-envs/adec7d399fa9492178d9caa96f93852d/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/trim_galore/se/environment.yaml
#   prefix: /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - trim-galore =0.6.10
RUN mkdir -p /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.13.8/bio/trim_galore/se/environment.yaml /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5/environment.yaml

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
#   prefix: /conda-envs/79d668cc963ea72d97f0971437f4ce14
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
#     - pysam
#     - gffpandas=1.2.2
#     - numpy=1.26.4
#     - samtools=1.13
RUN mkdir -p /conda-envs/79d668cc963ea72d97f0971437f4ce14
COPY workflow/envs/deeptools.yaml /conda-envs/79d668cc963ea72d97f0971437f4ce14/environment.yaml

# Conda environment:
#   source: workflow/envs/peak_calling.yaml
#   prefix: /conda-envs/35dbd4f499c2bf378ca66a63855a13dc
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
#     - macs2=2.2.9.1
RUN mkdir -p /conda-envs/35dbd4f499c2bf378ca66a63855a13dc
COPY workflow/envs/peak_calling.yaml /conda-envs/35dbd4f499c2bf378ca66a63855a13dc/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/79f24222e503aa3928bb131a6d96239b --file /conda-envs/79f24222e503aa3928bb131a6d96239b/environment.yaml && \
    mamba env create --prefix /conda-envs/84012557432533881e14141abce4094f --file /conda-envs/84012557432533881e14141abce4094f/environment.yaml && \
    mamba env create --prefix /conda-envs/90ea9a0a5f4de478eb6aa43932774175 --file /conda-envs/90ea9a0a5f4de478eb6aa43932774175/environment.yaml && \
    mamba env create --prefix /conda-envs/df61192e96fe9bc3d27cde8b6a0fa36b --file /conda-envs/df61192e96fe9bc3d27cde8b6a0fa36b/environment.yaml && \
    mamba env create --prefix /conda-envs/b29db9b4576aa24f6ca5cdabda9911b4 --file /conda-envs/b29db9b4576aa24f6ca5cdabda9911b4/environment.yaml && \
    mamba env create --prefix /conda-envs/adec7d399fa9492178d9caa96f93852d --file /conda-envs/adec7d399fa9492178d9caa96f93852d/environment.yaml && \
    mamba env create --prefix /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5 --file /conda-envs/7c8f9db6001702ef40e5b666df2ea0f5/environment.yaml && \
    mamba env create --prefix /conda-envs/4554a2c6669fb2054cdfbe43af1d6154 --file /conda-envs/4554a2c6669fb2054cdfbe43af1d6154/environment.yaml && \
    mamba env create --prefix /conda-envs/cb2a9c9685efe98f401d666d5b884a76 --file /conda-envs/cb2a9c9685efe98f401d666d5b884a76/environment.yaml && \
    mamba env create --prefix /conda-envs/79d668cc963ea72d97f0971437f4ce14 --file /conda-envs/79d668cc963ea72d97f0971437f4ce14/environment.yaml && \
    mamba env create --prefix /conda-envs/35dbd4f499c2bf378ca66a63855a13dc --file /conda-envs/35dbd4f499c2bf378ca66a63855a13dc/environment.yaml && \
    mamba clean --all -y
