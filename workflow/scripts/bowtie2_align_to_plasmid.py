"""
Align to plasmid index and only keep non-aligned reads in new fastq file(s)
"""

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

# Load Snakemake variables
fastq = snakemake.input["fastq"]
idx = snakemake.params["idxdir"]
paired = snakemake.params["paired"]
extra = snakemake.params["extra"]
bam = snakemake.output["bam"]
out_fastq = snakemake.output["fastq"]

# Prepare input and output fastq file names
if paired:
    r1 = fastq[0]
    r2 = fastq[1]
    fastq_in = f"-1 {r1} -2 {r2}"
    out_base = snakemake.params["out_base"]
    fastq_out = f"--un-conc-gz {out_base}"
else:
    fastq_in = f"-U {fastq[0]}"
    fastq_out = f"--un-gz {out_fastq}"

shell(
    "(bowtie2 "
    "-x {idx} "
    "{fastq_in} "
    "-p {snakemake.threads} "
    "{extra} "
    "{fastq_out} "
    "| samtools view -F 4 --with-header "  # exclude unmapped reads
    "> {bam} )"
    "{log}"
)
