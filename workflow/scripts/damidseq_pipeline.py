"""
damidseq_pipeline will be run for one replicate directory, split by condition (genotype_treatment), where the appropriate Dam only control is matched to non-Dam samples for the same condition.

Note: other replicate(s) will be run in parallel by Snakemake.

Problem: parallel analysis of replicate(s) will generate same output files in same directory.
Solution: work in a temporary directory and move files to appropriate locations at the end.
"""

import os
from pathlib import Path
import tempfile
import pandas as pd
import glob
from snakemake.shell import shell

def move_files(extension):
    """Move files to appropriate locations
    """
    # Create output directory for requested file type
    out_dir = os.path.join(cwd, f"results/{extension}/{directory}")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # All output has been made in main directory, so move it to appropriate locations
    shell(
        "mv *.{extension} {out_dir}"
        )

def allfastq(directory):
    """Returns absolute file paths for all fastq.gz files in a directory.
    Based on https://stackoverflow.com/questions/9816816/get-absolute-paths-of-all-files-in-a-directory
    """
    path = os.path.abspath(directory)
    all_files= [entry.path for entry in os.scandir(path) if entry.is_file()]
    all_fastq = [x for x in all_files if x.endswith(".fastq.gz")]
    return all_fastq

# Get current working directory
cwd = os.getcwd()

# Load Snakemake variables
trim_dir = snakemake.params["trim_dir"]
flags = snakemake.input["flag"]
gatc = os.path.join(cwd, snakemake.input["gatc"])
bowtie2_idx = os.path.join(cwd, snakemake.params["idxdir"])
paired = snakemake.params["paired"]
threads = snakemake.threads
bins = snakemake.params["binsize"]
extra = snakemake.params["extra"]

# Get sample directory
directory = list(set([os.path.basename(os.path.dirname(x)) for x in flags]))
assert len(directory) == 1, "Too many replicate directories used..."
directory = directory[0]

# Path to damidseq_pipeline script
damidseq_pipeline = os.path.join(cwd,"resources/damidseq_pipeline/damidseq_pipeline")

# Load sample table
csv = pd.read_csv("config/samples.csv")

# Check if treatment column contains any NaN values, if so replace with "none"
if csv["treatment"].isnull().values.any():
    csv.fillna({"treatment": "none"}, inplace=True)

# Combine genotypes and treatments into one condition column
csv["condition"] = csv["genotype"] + "_" + csv["treatment"]

# Create sample table without Dam only samples
csv_no_dam = csv[~csv["sample"].str.contains("Dam")]

# For each unique condition, find Dam only control file
conditions = csv["condition"].unique().tolist()
dam_controls = {} # Has dam control for each condition
for c in conditions:
    # Get Dam sample that matches condition in csv and add to dam_controls, i.e. Dam in sample column
    dam_controls[c] = csv[csv["sample"].str.contains("Dam") & csv["condition"].str.contains(c)]["sample"].tolist()[0]

# Get all fastq files
all_fastq = allfastq(f"results/{trim_dir}/{directory}/")

# Run damidseq_pipeline for each condition
for condition, dam_control in dam_controls.items():
    print(f"Analysing data in results/trimmed/{directory}...\nUsing {dam_control} as control")
    
    # Create temporary directory and move there
    temp_dir = tempfile.TemporaryDirectory()
    print(f"Creating temporary directory {temp_dir.name}")
    os.chdir(temp_dir.name)
       
    # Define Dam only control file(s)
    if paired:
        arg = "--paired"
        dam = sorted([x for x in all_fastq if f"{dam_control}_1.fastq.gz" or f"{dam_control}_2.fastq.gz" in x])
        dam = " ".join(dam)
    else:
        arg = ""
        dam = [x for x in all_fastq if f"{dam_control}.fastq.gz" in x][0]
        
    # Get all non-Dam only files that match condition
    INPUT = csv_no_dam[csv_no_dam["condition"].str.contains(condition)]["sample"].tolist()
    INPUT = sorted([x for x in all_fastq if any([y in x for y in INPUT])])
    INPUT = " ".join(INPUT)
    assert len(INPUT) > len(all_fastq), "No trimmed Dam only files found..."
    
    # Run damidseq_pipeline
    shell(
        "perl {damidseq_pipeline} "
        "{arg} "
        "--threads={threads} "
        "--dam={dam} "
        "--bins={bins}"
        "{extra} "
        "--gatc_frag_file={gatc} "
        "--bowtie2_genome_dir={bowtie2_idx} "
        "{INPUT} "
        )

    print("Moving output files from temporary directory to appropriate locations")
    # Move log file to logs directory
    target = os.path.join(cwd, f"logs/damidseq_pipeline/{directory}")
    shell(
        "mv pipeline-*.log {target}"
        )

    # Move bedgraph files to output directory
    move_files("bedgraph")

    # Move bam files to output directory
    move_files("bam")
    
    if not paired:
        # Rename bam files so that they end with .bam and not -ext300.bam
        bam_files = glob.glob(f"results/bam/{directory}/*-ext300.bam")
        assert len(bam_files) > 0, "No single-end bam files found to remove -ext300 from file name..."
        for bam in bam_files:
            new_name = bam.replace("-ext300", "")
            os.rename(bam, new_name)
    
    # Destroy temporary directory
    print(f"Cleaning up temporary directory {temp_dir.name}")
    temp_dir.cleanup()
