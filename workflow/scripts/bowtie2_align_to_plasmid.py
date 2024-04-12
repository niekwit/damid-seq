from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Load Snakemake variables
flag = snakemake.input["flag"]
idx = snakemake.params["idxdir"]
paired = snakemake.params["paired"]
extra = snakemake.params["extra"]
out_base = snakemake.params["out_base"]

# Prepare input and output fastq file flags
base = flag.replace(".flag", "")
if paired:
    r1 = base + "_1.fastq.gz"
    r2 = base + "_2.fastq.gz"
    fastq_in = f"-1 {r1} -2 {r2}"
    fastq_out = f"--un-conc-gz {out_base}"
else:
    fastq_in = f"-U {base + '.fastq.gz'}"
    fastq_out = f"--un-gz {out_base + '.fastq.gz'}"
    
# Align to plasmid index and only keep non-aligned reads in new fastq file(s)
shell(
    "bowtie2 -x {idx} {fastq_in} -p {snakemake.threads} {extra} {fastq_out} > /dev/null {log}"
    )