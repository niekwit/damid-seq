"""
damidseq_pipeline will be run for one replicate directory, split by condition (genotype_treatment), where the appropriate Dam only control is matched to non-Dam samples for the same condition.

Note: other replicate(s) will be run in parallel by Snakemake.

Problem: parallel analysis of replicate(s) will generate same output files in same directory.
Solution: work in a temporary directory and move files to appropriate locations at the end.
"""

import logging
import os
from pathlib import Path
import tempfile
import re
import glob
from snakemake.shell import shell


def move_files(extension):
    """Move files to appropriate locations"""
    # Create output directory for requested file type
    out_dir = os.path.join(cwd, f"results/{extension}/{directory}")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # All output has been made in main directory, so move it to appropriate locations
    shell("mv *.{extension} {out_dir}")


# Get current working directory
cwd = os.getcwd()

# Load Snakemake variables
bams = snakemake.input["bam"]
bams = [os.path.join(cwd, x) for x in bams]
gatc = os.path.join(cwd, snakemake.input["gatc"])
threads = snakemake.threads
bins = snakemake.params["binsize"]
normalization_method = snakemake.params["normalization_method"]
idx = os.path.join(cwd, snakemake.params["idx"])
extra = snakemake.params["extra"]
log = os.path.join(cwd, snakemake.log[0])

# Setup log file
logging.basicConfig(
    format="%(levelname)s:%(message)s",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
)

# Get sample directory
directory = list(set([os.path.basename(os.path.dirname(x)) for x in bams]))
assert len(directory) == 1, "Too many replicate directories used..."
directory = directory[0]

# Path to damidseq_pipeline script
damidseq_pipeline = os.path.join(cwd, "resources/damidseq_pipeline/damidseq_pipeline")

# Get dam bam
dam_bam = [x for x in bams if x.lower().endswith(f"dam.sorted.bam")]
assert len(dam_bam) == 1, "No Dam only bam file found..."
# TO DO: allow for multiple Dam only bam files (i.e. multiple conditions)

# Get non-Dam bam files
non_dam_bams = [x for x in bams if x != dam_bam[0]]

logging.info(f"Analysing data in results/bam/{directory}...")

# Create temporary directory and move there
temp_dir = tempfile.TemporaryDirectory()
logging.info(f"Creating temporary directory {temp_dir.name}")
os.chdir(temp_dir.name)

# Run damidseq_pipeline
command = [
    f"perl {damidseq_pipeline} "
    f"--threads={threads} "
    f"--dam={dam_bam[0]} "
    f"--bins={bins} "
    f"--norm_method={normalization_method} "
    f"{extra} "
    f"--gatc_frag_file={gatc} "
    f"--bowtie2_genome_dir={idx} "  # error if not parsed
    f"{' '.join(non_dam_bams)} "
]
logging.info(f"Running damidseq_pipeline with command: {' '.join(command)}")
shell(" ".join(command))

logging.info("Moving output files from temporary directory to appropriate locations")
# Move log file to logs directory
target = os.path.join(cwd, f"logs/damidseq_pipeline/{directory}")
shell("mv pipeline-*.log {target}")

# Move bedgraph files to output directory
move_files("bedgraph")

# Go to parent directory
os.chdir(cwd)

# Rename bedgraph files so that normalization method is not included in file name
bedgraphs = glob.glob(f"results/bedgraph/{directory}/*.bedgraph")
assert len(bedgraphs) > 0, "No bedgraph files found to rename..."
for bedgraph in bedgraphs:
    new_name = bedgraph.replace(
        f".{normalization_method}-norm.gatc.bedgraph", "-norm.gatc.bedgraph"
    )
    os.rename(bedgraph, new_name)

# Destroy temporary directory
logging.info(f"Cleaning up temporary directory {temp_dir.name}")
temp_dir.cleanup()
