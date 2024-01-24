from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

gff_files = snakemake.input["gff"]
samples = snakemake.wildcards["bg_sample"]
out_file = snakemake.output[0]
min_overlap = snakemake.params["min_overlap"]
extra = snakemake.params["extra"]

# For each GFF file that matches a value for bg_sample wildcard get common peaks
for sample in samples:
    # get all GFF files that match sample
    gffs = sorted([x for x in gff_files if sample in x])
    
    # Get first GFF file
    a_gff = gffs[0]
    
    # Get all other GFF files
    b_gffs = gffs[1:]
    
    # Get common peaks with bedtools intersect
    shell(
        "bedtools intersect "
        "-a {a_gff} "
        "-b {b_gffs} "
        "{extra} "
        "-r "
        "-f {min_overlap} > {out_file} "
        "{log} "
    )

