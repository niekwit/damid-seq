import glob
import os
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Get arguments from snakemake
find_peaks = f"{snakemake.input['fp']}/find_peaks"
bedgraph = snakemake.input["bg"]
gff = snakemake.output["gff"]
data = snakemake.output["data"]
fdr = snakemake.params["fdr"]
frac = snakemake.params["frac"]
min_count = snakemake.params["mc"]
min_quant = snakemake.params["mq"]
n = snakemake.params["n"]
step = snakemake.params["step"]
up = snakemake.params["up"]
outdir = snakemake.params["out"]

# Run find_peaks
shell(
    "perl {find_peaks} "
    "--fdr={fdr} "
    "--frac={frac} "
    "--min_count={min_count} "
    "--min_quant={min_quant} "
    "--n={n} "
    "--step={step} "
    "--unified_peaks={up} "
    "{bedgraph} "
    "{log}"
    )

# Files are created in bedgraph directory, so move to peak directory
# Problem: folder containing files contains date/time
# Solution: use glob to find files
bg_dir = os.path.dirname(bedgraph)
gff_temp = glob.glob(f"{bg_dir}/*.gff")[0] # peak file
data_temp = glob.glob(f"{bg_dir}/*-data")[0] # data file

# Rename and move files
os.makedirs(outdir, exist_ok=True)
os.replace(gff_temp, gff)
os.replace(data_temp, data)

# Remove empty directory
os.rmdir(bg_dir)

