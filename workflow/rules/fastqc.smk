if paired_end:
    rule fastqc:
        input:
            "reads/{dir}/{sample}{end}.fastq.gz",
        output:
            html="results/qc/fastqc/{dir}/{sample}{end}.html",
            zip="results/qc/fastqc/{dir}/{sample}{end}_fastqc.zip"
        params:
            extra = "--quiet"
        log:
            "logs/fastqc/{dir}/{sample}{end}.log"
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        wrapper:
            "v3.3.4/bio/fastqc"


    rule multiqc:
        input:
            expand("results/qc/fastqc/{dir}/{sample}{end}_fastqc.zip", dir=DIRS,sample=SAMPLES, end=["_R1_001","_R2_001"])
        output:
            r="results/qc/multiqc/multiqc.html",
            d=directory("results/qc/multiqc/"),
        params:
            extra="",  # Optional: extra parameters for multiqc
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        log:
            "logs/multiqc/multiqc.log"
        conda:
            "../envs/trim.yaml"
        shell:
            "multiqc " 
            "--force "
            "--outdir {output.d} "
            "-n multiqc.html "
            "{params.extra} "
            "{input} "
            "> {log} 2>&1"

else:
    rule fastqc:
        input:
            "reads/{dir}/{sample}.fastq.gz"
        output:
            html="results/qc/fastqc/{dir}/{sample}.html",
            zip="results/qc/fastqc/{dir}/{sample}_fastqc.zip"
        params:
            extra = "--quiet"
        log:
            "logs/fastqc/{dir}/{sample}.log"
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        wrapper:
            "v3.3.4/bio/fastqc"


    rule multiqc:
        input:
            expand("results/qc/fastqc/{dir}/{sample}_fastqc.zip", dir=DIRS, sample=SAMPLES)
        output:
            r="results/qc/multiqc/multiqc.html",
            d=directory("results/qc/multiqc/"),
        params:
            extra="", # Optional: extra parameters for multiqc
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        log:
            "logs/multiqc/multiqc.log"
        conda:
            "../envs/trim.yaml"
        shell:
            "multiqc " 
            "--force "
            "--outdir {output.d} "
            "-n multiqc.html "
            "{params.extra} "
            "{input} "
            "> {log} 2>&1"
            