import os
import re
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
name = re.sub(".bam$", "", os.path.basename(bam))

# Load sample table
csv = pd.read_csv("config/samples.csv")

# Check if treatment column contains any NaN values, if so replace with "none"
if csv["treatment"].isnull().values.any():
    csv.fillna({"treatment": "none"}, inplace=True)

# Combine genotypes and treatments into one condition column
csv["condition"] = csv["genotype"] + "_" + csv["treatment"]

# Get condition for name
condition = csv[csv["sample"] == name]["condition"].tolist()[0]

# Get dam sample that matches bam sample (name) condition in csv
dam = csv[csv["sample"].str.contains("Dam")]

if len(dam) == 1:
    dam = dam["sample"].tolist()[0]
else:
    dam = dam[csv["condition"].str.contains(condition)]["sample"].tolist()[0]
dam = os.path.join(data_dir, dam + ".bam")

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