#!/usr/bin/env bash

set -e

# get the current working directory
WORKDIR=$(pwd)

# go to directory with fastq files
SAMPLEDIR=${snakemake_wildcards[dir]}
cd "reads/${SAMPLEDIR}"

# Check if data is paired-end or single-end
END=${snakemake_params[paired]}

if [ "$END" == "True" ]; then
    ARG="--paired"
else
    ARG=""
fi

# run DamID-seq pipeline
damidseq_pipeline $ARG --gatc_frag_file=../${snakemake_input[gatc]} --bowtie2_genome_dir=../${snakemake_params[idxdir]} > ../${snakemake_log[0]} 2>&1

# go back to working directory
cd ${WORKDIR}

# move output files to results/sample/bedgraph
mkdir -p ${snakemake_output[dir]}
mv *.bedgraph ${snakemake_output[dir]}
#mkdir -p logs/damid_pipeline/
#mv pipeline-*.log ${WORKDIR}/logs/damid_pipeline/


