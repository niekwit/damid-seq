if paired_end:
    rule post_trim_fastqc:
        input:
            "results/trimmed/{dir}/{sample}{end}.fastq.gz",
        output:
            html="results/qc/fastqc/{dir}/{sample}{end}_post_trim.html",
            zip="results/qc/fastqc/{dir}/{sample}{end}_post_trim_fastqc.zip"
        params:
            extra = "--quiet"
        log:
            "logs/fastqc/{dir}/{sample}{end}.log"
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        wrapper:
            f"{wrapper_version}/bio/fastqc"


    rule pre_trim_fastqc:
        input:
            "reads/{dir}/{sample}{end}.fastq.gz",
        output:
            html="results/qc/fastqc/{dir}/{sample}{end}_pre_trim.html",
            zip="results/qc/fastqc/{dir}/{sample}{end}_pre_trim_fastqc.zip"
        params:
            extra = "--quiet"
        log:
            "logs/fastqc/{dir}/{sample}{end}.log"
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        wrapper:
            f"{wrapper_version}/bio/fastqc"


    rule multiqc:
        input:
            expand("results/qc/fastqc/{dir}/{sample}{end}{trim}_fastqc.zip", dir=DIRS, sample=SAMPLES, end=["_1","_2"], trim=["_pre_trim","_post_trim"])
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
            "--dirs " # Prepend directory to sample names
            "-n multiqc.html "
            "{params.extra} "
            "{input} "
            "> {log} 2>&1"
else:
    rule post_trim_fastqc:
        input:
            "results/trimmed/{dir}/{sample}.fastq.gz"
        output:
            html="results/qc/fastqc/{dir}/post_trim/{sample}.html",
            zip="results/qc/fastqc/{dir}/post_trim/{sample}_fastqc.zip"
        params:
            extra = "--quiet"
        log:
            "logs/fastqc/{dir}/{sample}.log"
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        wrapper:
            f"{wrapper_version}/bio/fastqc"


    rule pre_trim_fastqc:
        input:
            "reads/{dir}/{sample}.fastq.gz",
        output:
            html="results/qc/fastqc/{dir}/pre_trim/{sample}.html",
            zip="results/qc/fastqc/{dir}/pre_trim/{sample}_fastqc.zip"
        params:
            extra = "--quiet"
        log:
            "logs/fastqc/{dir}/{sample}.log"
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        wrapper:
            f"{wrapper_version}/bio/fastqc"


    rule multiqc:
        input:
            expand("results/qc/fastqc/{dir}/{trim}/{sample}_fastqc.zip", dir=DIRS, sample=SAMPLES, trim=["pre_trim","post_trim"])
        output:
            "results/qc/multiqc/multiqc.html",
        params:
            extra="--dirs --dirs-depth 2", # Optional: extra parameters for multiqc
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        log:
            "logs/multiqc/multiqc.log"
        wrapper:
            f"{wrapper_version}/bio/multiqc"