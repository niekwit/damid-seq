from snakemake.shell import shell

log = snakemake.log[0]

genome = snakemake.params["genome"]

if "hg" in genome:
    shell(
        "perl $CONDA_PREFIX/share/homer/configureHomer.pl "
        "-install human-p > {log} 2>&1"
    )
elif "mm" in genome:
    shell(
        "perl $CONDA_PREFIX/share/homer/configureHomer.pl "
        "-install mouse-p > {log} 2>&1"
    )
else:
    print(f"Genome {genome} not supported")
    exit(1)

shell(
    "perl $CONDA_PREFIX/share/homer/configureHomer.pl "
    "-install {genome} >> {log} 2>&1"
)

