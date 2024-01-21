#!/usr/bin/env bash

set -e


#TEMP_DIR=${snakemake_output[temp_dir]}
#TEMP_DIR=temp/${snakemake_wildcards[dir]}/${snakemake_wildcards[sample]}
DEST_DIR=results/trimmed/${snakemake_wildcards[dir]}
BASENAME=${snakemake_wildcards[sample]}

#mkdir -p $TEMP_DIR
mkdir -p $DEST_DIR

if [ ${snakemake_params[paired]} == "YES" ]
    then
    INPUT="${snakemake_input[r1]} ${snakemake_input[r2]}"
else 
    INPUT=${snakemake_input[r1]}
fi

trim_galore ${snakemake_params[extra]} --cores ${snakemake[threads]} --output_dir $DEST_DIR --basename $BASENAME $INPUT > ${snakemake_log[0]} 2>&1

#mv $TEMP_DIR/*_trimmed.fq.gz $DEST_DIR
#mv $TEMP_DIR/*trimming_report.txt $DEST_DIR

#rm -r $TEMP_DIR

