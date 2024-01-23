if paired_end:
    rule trim_galore_pe:
        input:
            r1="reads/{dir}/{sample}_R1_001.fastq.gz", 
            r2="reads/{dir}/{sample}_R2_001.fastq.gz",
        output:
            r1="results/trimmed/{dir}/{sample}_1.fastq.gz",
            #report_fwd="results/trimmed/{dir}/{sample}_R1_trimming_report.txt",
            r2="results/trimmed/{dir}/{sample}_2.fastq.gz",
            #report_rev="results/trimmed/{dir}/{sample}_R2_trimming_report.txt",
            #temp_dir=temp(directory("temp/{dir}/{sample}")),
            #dest_dir=directory("results/trimmed/{dir}/{sample}"),
            flag=touch("results/trimmed/{dir}/{sample}.flag"),
        threads: config["resources"]["trim"]["cpu"],
        resources:
            runtime=config["resources"]["trim"]["time"],
        params:
            paired=True,
            extra="--illumina -q 20",
        log:
            "logs/trim_galore/{dir}/{sample}.log",
        conda:
            "../envs/trim.yaml"
        script:
            "../scripts/trim_galore.py"
else:
    rule trim_galore_se:
        input:
            r1="reads/{dir}/{sample}.fastq.gz",
        output:
            r1="results/trimmed/{dir}/{sample}.fastq.gz",
            #report_fwd="results/trimmed/{dir}/{sample}_trimming_report.txt",
            #temp_dir=temp(directory("temp/{dir}/{sample}")),
            #dest_dir=directory("results/trimmed/{dir}/{sample}"),
            flag=touch("results/trimmed/{dir}/{sample}.flag"),
        threads: config["resources"]["trim"]["cpu"],
        resources:
            runtime=config["resources"]["trim"]["time"],
        params:
            paired=False,
            extra="--illumina -q 20",
        log:
            "logs/trim_galore/{dir}/{sample}.log",
        conda:
            "../envs/trim.yaml"
        script:
            "../scripts/trim_galore.py"

        