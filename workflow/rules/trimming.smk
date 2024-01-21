if config["paired_end"]:
    rule trim_galore_pe:
        input:
            r1="reads/{dir}/{sample}_R1_001.fastq.gz", 
            r2="reads/{dir}/{sample}_R2_001.fastq.gz",
        output:
            fastq_fwd="results/trimmed/{dir}/{sample}_val_1.fq.gz",
            #report_fwd="results/trimmed/{dir}/{sample}_R1_trimming_report.txt",
            fastq_rev="results/trimmed/{dir}/{sample}_val_2.fq.gz",
            #report_rev="results/trimmed/{dir}/{sample}_R2_trimming_report.txt",
            #temp_dir=temp(directory("temp/{dir}/{sample}")),
            #dest_dir=directory("results/trimmed/{dir}/{sample}"),
            flag=touch("results/trimmed/{dir}/{sample}.flag"),
        threads: config["resources"]["trim"]["cpu"],
        resources:
            runtime=config["resources"]["trim"]["time"],
        params:
            paired="YES",
            extra="--illumina -q 20",
        log:
            "logs/trim_galore/{dir}/{sample}.log",
        conda:
            "../envs/damid.yaml"
        script:
            "../scripts/trim_galore.sh"
else:
    rule trim_galore_se:
        input:
            r1="reads/{dir}/{sample}.fastq.gz",
        output:
            fastq="results/trimmed/{dir}/{sample}_trimmed.fq.gz",
            #report_fwd="results/trimmed/{dir}/{sample}_trimming_report.txt",
            #temp_dir=temp(directory("temp/{dir}/{sample}")),
            #dest_dir=directory("results/trimmed/{dir}/{sample}"),
            flag=touch("results/trimmed/{dir}/{sample}.flag"),
        threads: config["resources"]["trim"]["cpu"],
        resources:
            runtime=config["resources"]["trim"]["time"],
        params:
            paired="NO",
            extra="--illumina -q 20",
        log:
            "logs/trim_galore/{dir}/{sample}.log",
        conda:
            "../envs/damid.yaml"
        script:
            "../scripts/trim_galore.sh"

        