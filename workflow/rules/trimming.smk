if paired_end:

    rule trim_galore_pe:
        input:
            [
                "reads/{dir}/{sample}_R1_001.fastq.gz",
                "reads/{dir}/{sample}_R2_001.fastq.gz",
            ],
        output:
            fasta_fwd="results/trimmed/{dir}/{sample}_1.fastq.gz",
            fasta_rev="results/trimmed/{dir}/{sample}_2.fastq.gz",
            report_fwd="results/trimmed/{dir}/{sample}_1.fastq.gz_trimming_report.txt",
            report_rev="results/trimmed/{dir}/{sample}_2.fastq.gz_trimming_report.txt",
        params:
            extra="--illumina -q 20",
        threads: config["resources"]["trim"]["cpu"]
        resources:
            runtime=config["resources"]["trim"]["time"],
        log:
            "logs/trim_galore/{dir}/{sample}.log",
        wrapper:
            "v5.8.3/bio/trim_galore/pe"

else:

    rule trim_galore_se:
        input:
            "reads/{dir}/{sample}.fastq.gz",
        output:
            fasta="results/trimmed/{dir}/{sample}.fastq.gz",
            report="results/trimmed/{dir}/{sample}.fastq.gz_trimming_report.txt",
        params:
            extra="--illumina -q 20",
        threads: config["resources"]["trim"]["cpu"]
        resources:
            runtime=config["resources"]["trim"]["time"],
        log:
            "logs/trim_galore/{dir}/{sample}.log",
        wrapper:
            "v5.8.3/bio/trim_galore/se"
