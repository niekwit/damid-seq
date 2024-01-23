import os
from pathlib import Path
import glob
from snakemake.shell import shell

# Load Snakemake variables
flags = snakemake.input["flag"]
gatc = snakemake.input["gatc"]
bowtie2_idx = snakemake.params["idxdir"]
paired = snakemake.params["paired"]
threads = snakemake.threads
extra = snakemake.params["extra"]

# Get sample directories
dirs = list(set([os.path.basename(os.path.dirname(x)) for x in flags]))

assert len(dirs) != 0, "No data found in results/trimmed/"

# Path to damidseq_pipeline script
damidseq_pipeline = "resources/damidseq_pipeline/damidseq_pipeline"

def move_files(extension):
    # Create output directory for requested file type
    out_dir = f"results/{extension}/{d}"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # All output has been made in main directory, so move it to appropriate locations
    # Move files to output directory
    shell(
        "mv *.{extension} "
        "{out_dir}"
        )

#TO DO: if multiple treatments/conditions exist, find appropriate Dam only for these individual conditions/treatments

# Run damidseq_pipeline for each set of trimmed fastq files
for d in dirs:
    print(f"Analysing data in results/trimmed/{d}...")
    
    # Define Dam only control file(s)
    if paired:
        arg = "--paired"
        dam = f"results/trimmed/{d}/Dam_R1.fastq.gz,results/trimmed/{d}/Dam_R2.fastq.gz"
    else:
        arg = ""
        dam = f"results/trimmed/{d}/Dam.fastq.gz"
     
    # Specify input files
    all_fastq = glob.glob(f"results/trimmed/{d}/*.fastq.gz")
    INPUT = " ".join([x for x in all_fastq if not "Dam" in x])
    assert len(INPUT) > len(all_fastq), "No trimmed Dam only files found..."
     
    # Run damidseq_pipeline
    shell(
        "perl {damidseq_pipeline} "
        "{arg} "
        "--threads={threads} "
        "--dam={dam} "
        "{extra} "
        "--gatc_frag_file={gatc} "
        "--bowtie2_genome_dir={bowtie2_idx} "
        "{INPUT} "
        )
    print("Moving output files to appropriate locations...")
    
    # Move log file to logs directory
    shell(
        "mv pipeline-*.log "
        "logs/damidseq_pipeline/{d}"
        )

    # Move bedgraph files to output directory
    move_files("bedgraph")
    
    # Move bam files to output directory
    move_files("bam")
    

