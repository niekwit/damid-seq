if config["plasmid_fasta"] == "none":
    if paired_end:

        rule bowtie2_align_pe:
            input:
                sample=[
                    "results/trimmed/{dir}/{sample}_1.fastq.gz",
                    "results/trimmed/{dir}/{sample}_2.fastq.gz",
                ],
                idx=multiext(
                    f"resources/bowtie2_index/{resources.genome}_{resources.build}_{maskedgenes}.masked/index",
                    ".1.bt2",
                    ".2.bt2",
                    ".3.bt2",
                    ".4.bt2",
                    ".rev.1.bt2",
                    ".rev.2.bt2",
                ),
            output:
                bam="results/bam/{dir}/{sample}.bam",
            log:
                "logs/bowtie2_align/{dir}/{sample}.log",
            threads: config["resources"]["bowtie2"]["cpu"]
            resources:
                runtime=config["resources"]["bowtie2"]["time"],
            params:
                extra=config["bowtie2"]["extra"],
            wrapper:
                "v5.8.3/bio/bowtie2/align"

        rule sort_bowtie2_bam:
            input:
                "results/bam/{dir}/{sample}.bam",
            output:
                "results/bam/{dir}/{sample}.sorted.bam",
            log:
                "logs/samtools/sort/{dir}/{sample}.log",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"],
            wrapper:
                "v5.8.3/bio/samtools/sort"

        rule index_bowtie2_bam:
            input:
                "results/bam/{dir}/{sample}.sorted.bam",
            output:
                "results/bam/{dir}/{sample}.sorted.bam.bai",
            log:
                "logs/samtools/index/{dir}/{sample}.log",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"],
            wrapper:
                "v5.8.3/bio/samtools/index"

    else:

        rule bowtie2_align_se:
            input:
                sample=["results/trimmed/{dir}/{sample}.fastq.gz"],
                idx=multiext(
                    f"resources/bowtie2_index/{resources.genome}_{resources.build}_{maskedgenes}.masked/index",
                    ".1.bt2",
                    ".2.bt2",
                    ".3.bt2",
                    ".4.bt2",
                    ".rev.1.bt2",
                    ".rev.2.bt2",
                ),
            output:
                bam="results/bam/{dir}/{sample}.bt2.bam",
            log:
                "logs/bowtie2_align/{dir}/{sample}.log",
            threads: config["resources"]["bowtie2"]["cpu"]
            resources:
                runtime=config["resources"]["bowtie2"]["time"],
            params:
                extra=config["bowtie2"]["extra"],
            wrapper:
                "v5.8.3/bio/bowtie2/align"

        rule sort_bowtie2_bam:
            input:
                "results/bam/{dir}/{sample}.bt2.bam",
            output:
                temp("results/bam/{dir}/{sample}.bt2_sorted.bam"),
            log:
                "logs/samtools/index/{dir}/{sample}.log",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"],
            wrapper:
                "v5.8.3/bio/samtools/sort"

        rule index_bowtie2_bam:
            input:
                "results/bam/{dir}/{sample}.bt2_sorted.bam",
            output:
                temp("results/bam/{dir}/{sample}.bt2_sorted.bam.bai"),
            log:
                "logs/samtools/index/{dir}/{sample}.log",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"],
            wrapper:
                "v5.8.3/bio/samtools/index"

        # Single-end fragments are extended to n bp
        # or to the closest GATC site
        # Fragments are also sorted
        rule extend_fragments:
            input:
                bam="results/bam/{dir}/{sample}.bt2_sorted.bam",
                bai="results/bam/{dir}/{sample}.bt2_sorted.bam.bai",
                gatc_gff=f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.GATC.gff",
                fai=f"{resources.fasta}.fai",
            output:
                bam="results/bam/{dir}/{sample}.sorted.bam",
            log:
                "logs/extend_reads/{dir}/{sample}.log",
            conda:
                "../envs/deeptools.yaml"
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"],
            params:
                n=300,
                q=20,
            shell:
                f"perl {os.path.join(workflow.basedir, 'scripts/extend_reads.pl')} "
                "{input.bam} "
                "{output.bam} "
                "{input.gatc_gff} "
                "{params.n} "
                "{params.q} "
                "{log}"

        rule index_bam:
            input:
                "results/bam/{dir}/{sample}.sorted.bam",
            output:
                "results/bam/{dir}/{sample}.sorted.bam.bai",
            log:
                "logs/samtools/index/{dir}/{sample}.log",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"],
            wrapper:
                "v5.8.3/bio/samtools/index"

else:
    logger.info(
        f"Removing reads aligning to {config["plasmid_fasta"]} from trimmed reads..."
    )
    check_plasmid()
    if paired_end:

        rule bowtie2_remove_plasmid_reads:
            input:
                fastq=[
                    "results/trimmed/{dir}/{sample}_1.fastq.gz",
                    "results/trimmed/{dir}/{sample}_2.fastq.gz",
                ],
                idx=multiext(
                    f"resources/bowtie2_index/{plasmid_name}/index",
                    ".1.bt2",
                    ".2.bt2",
                    ".3.bt2",
                    ".4.bt2",
                    ".rev.1.bt2",
                    ".rev.2.bt2",
                ),
            output:
                fastq=[
                    "results/trimmed_no_plasmid/{dir}/{sample}_1.fastq.gz",
                    "results/trimmed_no_plasmid/{dir}/{sample}_2.fastq.gz",
                ],
                bam="results/bam/plasmid_reads/{dir}/{sample}.bam",
            log:
                "logs/bowtie2_align_to_plasmid/{dir}/{sample}.log",
            conda:
                "../envs/damid.yaml"
            threads: config["resources"]["damid"]["cpu"]
            resources:
                runtime=config["resources"]["damid"]["time"],
            params:
                idxdir=lambda wildcards, input: input["idx"][0][:-6],
                paired=paired_end,
                out_base=lambda wildcards, output: output["r1_fastq"].replace(
                    "_1.fastq.gz", "_%.fastq.gz"
                ),
                extra=config["bowtie2"]["extra"],
            script:
                "../scripts/bowtie2_align_to_plasmid.py"

        rule bowtie2_align:
            input:
                sample=[
                    "results/trimmed_no_plasmid/{dir}/{sample}_1.fastq.gz",
                    "results/trimmed_no_plasmid/{dir}/{sample}_2.fastq.gz",
                ],
                idx=multiext(
                    f"resources/bowtie2_index/{resources.genome}_{resources.build}_{maskedgenes}.masked/index",
                    ".1.bt2",
                    ".2.bt2",
                    ".3.bt2",
                    ".4.bt2",
                    ".rev.1.bt2",
                    ".rev.2.bt2",
                ),
            output:
                bam=temp(expand("results/bam/{{dir}}/{sample}.bam", sample=SAMPLES)),
            log:
                "logs/bowtie2_align/{dir}/{sample}.log",
            threads: config["resources"]["bowtie2"]["cpu"]
            resources:
                runtime=config["resources"]["bowtie2"]["time"],
            params:
                extra=config["bowtie2"]["extra"],
            wrapper:
                "v5.8.3/bio/bowtie2/align"

    else:

        rule bowtie2_remove_plasmid_reads:
            input:
                fastq=["results/trimmed/{dir}/{sample}.fastq.gz"],
                idx=multiext(
                    f"resources/bowtie2_index/{plasmid_name}/index",
                    ".1.bt2",
                    ".2.bt2",
                    ".3.bt2",
                    ".4.bt2",
                    ".rev.1.bt2",
                    ".rev.2.bt2",
                ),
            output:
                bam="results/bam/plasmid_reads/{dir}/{sample}.bam",
                fastq=temp("results/trimmed_no_plasmid/{dir}/{sample}.fastq.gz"),
            log:
                "logs/bowtie2_align_to_plasmid/{dir}/{sample}.log",
            conda:
                "../envs/damid.yaml"
            threads: config["resources"]["damid"]["cpu"]
            resources:
                runtime=config["resources"]["damid"]["time"],
            params:
                idxdir=lambda wildcards, input: input["idx"][0][:-6],
                paired=paired_end,
                out_base=lambda wildcards, output: output["fastq"].replace(
                    ".fastq.gz", ""
                ),
                extra=config["bowtie2"]["extra"],
            script:
                "../scripts/bowtie2_align_to_plasmid.py"

        rule bowtie2_align_se:
            input:
                sample=["results/trimmed_no_plasmid/{dir}/{sample}.fastq.gz"],
                idx=multiext(
                    f"resources/bowtie2_index/{resources.genome}_{resources.build}_{maskedgenes}.masked/index",
                    ".1.bt2",
                    ".2.bt2",
                    ".3.bt2",
                    ".4.bt2",
                    ".rev.1.bt2",
                    ".rev.2.bt2",
                ),
            output:
                bam=temp("results/bam/{dir}/{sample}.bt2.bam"),
            log:
                "logs/bowtie2_align/{dir}/{sample}.log",
            threads: config["resources"]["bowtie2"]["cpu"]
            resources:
                runtime=config["resources"]["bowtie2"]["time"],
            params:
                extra=config["bowtie2"]["extra"],
            wrapper:
                "v5.8.3/bio/bowtie2/align"

        rule sort_bowtie2_bam:
            input:
                "results/bam/{dir}/{sample}.bt2.bam",
            output:
                temp("results/bam/{dir}/{sample}.bt2_sorted.bam"),
            log:
                "logs/samtools/index/{dir}/{sample}.log",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"],
            wrapper:
                "v5.8.3/bio/samtools/sort"

        rule index_bowtie2_bam:
            input:
                "results/bam/{dir}/{sample}.bt2_sorted.bam",
            output:
                temp("results/bam/{dir}/{sample}.bt2_sorted.bam.bai"),
            log:
                "logs/samtools/index/{dir}/{sample}.log",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"],
            wrapper:
                "v5.8.3/bio/samtools/index"

        # Single-end fragments are extended to n bp
        # or to the closest GATC site
        # Fragments are also sorted
        rule extend_fragments:
            input:
                bam="results/bam/{dir}/{sample}.bt2_sorted.bam",
                bai="results/bam/{dir}/{sample}.bt2_sorted.bam.bai",
                gatc_gff=f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.GATC.gff",
                fai=f"{resources.fasta}.fai",
            output:
                bam="results/bam/{dir}/{sample}.sorted.bam",
            log:
                "logs/extend_reads/{dir}/{sample}.log",
            conda:
                "../envs/deeptools.yaml"
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"],
            params:
                n=300,
                q=20,
            shell:
                f"perl {os.path.join(workflow.basedir, 'scripts/extend_reads.pl')} "
                "{input.bam} "
                "{output.bam} "
                "{input.gatc_gff} "
                "{params.n} "
                "{params.q} "
                "{log}"

        rule index_bam:
            input:
                "results/bam/{dir}/{sample}.sorted.bam",
            output:
                "results/bam/{dir}/{sample}.sorted.bam.bai",
            log:
                "logs/samtools/index/{dir}/{sample}.log",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"],
            wrapper:
                "v5.8.3/bio/samtools/index"


# dir wildcard is escaped so that multiple dirs can be run in parallel
rule damidseq_pipeline:  # Ignore dir wildcard in expand statement (double braces)
    input:
        git="resources/damidseq_pipeline",
        bam=expand("results/bam/{{dir}}/{sample}.sorted.bam", sample=SAMPLES),
        bai=expand("results/bam/{{dir}}/{sample}.sorted.bam.bai", sample=SAMPLES),
        gatc=f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.GATC.gff",
        idx=multiext(
            f"resources/bowtie2_index/{resources.genome}_{resources.build}_{maskedgenes}.masked/index",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        bf=expand(
            "results/bedgraph/{{dir}}/{bg_sample}-vs-Dam-norm.gatc.bedgraph",
            bg_sample=BG_SAMPLES,
        ),
    log:
        "logs/damidseq_pipeline/{dir}/damidseq_pipeline.log",
    conda:
        "../envs/damid.yaml"
    threads: config["resources"]["damid"]["cpu"]
    resources:
        runtime=config["resources"]["damid"]["time"],
        tmpdir=config["resources"]["damid"]["tmpdir"],
    params:
        idx=lambda wildcards, input: input["idx"][0][:-6],
        binsize=config["damidseq_pipeline"]["binsize"],
        normalization_method=config["damidseq_pipeline"]["normalization"],
        extra=config["damidseq_pipeline"]["extra"],
    script:
        "../scripts/damidseq_pipeline.py"


rule bam2bigwig:
    input:
        bam="results/bam/{dir}/{sample}.sorted.bam",
        bai="results/bam/{dir}/{sample}.sorted.bam.bai",
    output:
        "results/bigwig/bam2bigwig/{dir}/{sample}.bw",
    log:
        "logs/bam2bigwig/{dir}/{sample}.log",
    conda:
        "../envs/deeptools.yaml"
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"],
    params:
        bs=config["deeptools"]["bamCoverage"]["binSize"],
        n=config["deeptools"]["bamCoverage"]["normalizeUsing"],
        extra=config["deeptools"]["bamCoverage"]["extra"],
    shell:
        "bamCoverage "
        "-b {input.bam} "
        "-o {output} "
        "--binSize {params.bs} "
        "-p {threads} "
        "--normalizeUsing {params.n} "
        "{params.extra} "
        "> {log} 2>&1 "


rule damidbind:
    input:
        bg=expand(
            "results/bedgraph/{dir}/{bg_sample}-vs-Dam-norm.gatc.bedgraph",
            bg_sample=BG_SAMPLES,
            dir=DIRS,
        ),
        gff=expand(
            "results/peaks/fdr{fdr}/{dir}/{bg_sample}.peaks.gff",
            bg_sample=BG_SAMPLES,
            dir=DIRS,
            fdr=fdr,
        ),
    output:
        diff_diagn_plot="results/damidbind/{comparison}/diagnostic_plots_diff.pdf",
        venn="results/damidbind/{comparison}/venn.pdf",
        volcano="results/damidbind/{comparison}/volcano.pdf",
        csv="results/damidbind/{comparison}/peaks.csv",
    log:
        "logs/damidbind/{comparison}/damidbind.log",
    conda:
        "../envs/damidbind.yaml"
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"],
    params:
        outdir_bg=lambda wildcards, output: os.path.join(
            os.path.dirname(output.csv), "bedgraph"
        ),
        outdir_peaks=lambda wildcards, output: os.path.join(
            os.path.dirname(output.csv), "peaks"
        ),
        genome=config["genome"],
        norm_method=config["differential_peaks"]["normalization"],
        fdr=config["differential_peaks"]["fdr"],
        filter_occupancy=config["differential_peaks"]["filter_occupancy"],
    script:
        "../scripts/damidbind.R"
