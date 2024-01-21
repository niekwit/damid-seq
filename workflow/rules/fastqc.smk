rule fastqc:
    input:
        "reads/{sample}{end}.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}{end}.html",
        zip="results/qc/fastqc/{sample}{end}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}{end}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    wrapper:
        "v3.3.4/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}{end}_fastqc.zip", sample=SAMPLES, end=["_R1","_R2"])
    output:
        "results/qc/multiqc.html",
        "results/qc/multiqc_data/multiqc_general_stats.txt"
    params:
        extra="",  # Optional: extra parameters for multiqc
    threads: config["resources"]["fastqc"]["cpu"]
    log:
        "logs/multiqc/multiqc.log"
    wrapper:
        "v3.3.4/bio/multiqc"