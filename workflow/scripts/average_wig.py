from snakemake.shell import shell

# Load Snakemake variables
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
threads = snakemake.threads

all_bw = snakemake.input
condition = snakemake.wildcards["bg_sample"]
wig = snakemake.output["wig"]

# Get all samples in condition
bw = [x for x in all_bw if condition in x] 

# Create average bigwig file
shell(
    "wiggletools write {wig} mean {bw} {log}"
    )
