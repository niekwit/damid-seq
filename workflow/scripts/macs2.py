import os
import re
import glob
import pandas as pd
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Load Snakemake variables
bam = snakemake.input["bam"]
fdr = snakemake.params["fdr"]
extra = snakemake.params["extra"]
mode = snakemake.params["mode"]
paired_end = snakemake.params["paired_end"]
genome = snakemake.params["genome"]
data_dir = snakemake.params["data_dir"]
out_dir = snakemake.params["outdir"]

#dam = os.path.join(data_dir, "Dam.bam")
name = re.sub(".sorted.bam$", "", os.path.basename(bam))

# Get Dam only bam file
dam = [x for x in glob.glob(os.path.join(data_dir, "*.sorted.bam")) if "dam" in x.lower()]
if len(dam) == 0:
    raise ValueError("No Dam only bam file found...")
elif len(dam) > 1:
    raise ValueError("Too many Dam only bam files found...")
else:
    dam = dam[0]

# Construct MACS2 arguments
if paired_end:
        data_format = "BAMPE"
else:
    data_format = "BAM"

if "hg" in genome:
    genome = "hs"
elif "mm" in genome:
    genome = "mm"
elif "dm" in genome:
    genome = "dm"

if mode == "broad":
    broad = f"--broad --broad-cutoff {fdr} "
    qvalue= ""
else:
    broad = ""
    qvalue = fdr
    qvalue = f"-q {fdr}"
    
# Run MACS2
shell(
    "macs2 callpeak "
    "--treatment {bam} "
    "--control {dam} "
    "--outdir {out_dir} "
    "--format {data_format} "
    "--name {name} "
    "-g {genome} "
    "{qvalue} "
    "{broad} "
    "{extra} "
    "{log}"
)