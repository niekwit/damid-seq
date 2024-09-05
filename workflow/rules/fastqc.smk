if paired_end:
    rule post_trim_fastqc:
        input:
            "results/trimmed/{dir}/{sample}_{end}.fastq.gz",
        output:
            html="results/qc/fastqc/{dir}/post_trim/{sample}_R{end}_001.html",
            zip="results/qc/fastqc/{dir}/post_trim/{sample}_R{end}_001_fastqc.zip"
        params:
            extra = "--quiet"
        log:
            "logs/fastqc/post_trim/{dir}/{sample}_R{end}_001.log"
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        wrapper:
            "v3.13.8/bio/fastqc"


    rule pre_trim_fastqc:
        input:
            "reads/{dir}/{sample}_R{end}_001.fastq.gz",
        output:
            html="results/qc/fastqc/{dir}/pre_trim/{sample}_R{end}_001.html",
            zip="results/qc/fastqc/{dir}/pre_trim/{sample}_R{end}_001_fastqc.zip"
        params:
            extra = "--quiet"
        log:
            "logs/fastqc/pre_trim/{dir}/{sample}_R{end}_001.log"
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        wrapper:
            "v3.13.8/bio/fastqc"


    rule multiqc:
        input:
            expand("results/qc/fastqc/{dir}/{trim}/{sample}_R{end}_001_fastqc.zip", dir=DIRS, sample=SAMPLES, end=["1","2"], trim=["pre_trim","post_trim"])
        output:
            r="results/qc/multiqc/multiqc.html",
        params:
            extra="--dirs --dirs-depth 2",  # Optional: extra parameters for multiqc
        threads: config["resources"]["fastqc"]["cpu"]
        resources:
            runtime=config["resources"]["fastqc"]["time"],
            mem_mb = 2048,
        log:
            "logs/multiqc/multiqc.log"
        wrapper:
            "v3.13.8/bio/multiqc"
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
            "v3.13.8/bio/fastqc"


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
            "v3.13.8/bio/fastqc"


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
            "v3.13.8/bio/multiqc"