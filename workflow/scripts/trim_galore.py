import os
from pathlib import Path
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Load Snakemake variables
flag = snakemake.output["flag"]
outdir = os.path.dirname(flag)
basename = snakemake.wildcards["sample"]
threads = snakemake.threads
extra = snakemake.params["extra"]

if snakemake.params["paired"]:
    INPUT = f"{snakemake.input['r1']} {snakemake.input['r2']}"
    extra = f"{extra} --paired"
else:
    INPUT = snakemake.input["r1"]

# Create output directory for trimmed fastq files
Path(outdir).mkdir(parents=True, exist_ok=True)

# Run trim_galore
shell(
    "trim_galore "
    "{extra} "
    "--cores {threads} "
    "--output_dir {outdir} "
    "--basename {basename} "
    "{INPUT} "
    "{log} "
)

# Rename output files
if snakemake.params["paired"]:
    os.rename(os.path.join(outdir, f"{basename}_val_1.fq.gz"), snakemake.output["r1"])
    os.rename(os.path.join(outdir, f"{basename}_val_2.fq.gz"), snakemake.output["r2"])
else:
    os.rename(os.path.join(outdir, f"{basename}_trimmed.fq.gz"), snakemake.output["r1"])
