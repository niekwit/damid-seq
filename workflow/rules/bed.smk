rule gff2bed:
    input:
        "results/peaks/{dir}/{bg_sample}.peaks.gff"
    output:
        "results/peaks/{dir}/{bg_sample}.bed"
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/gff2bed/{dir}_{bg_sample}.log"
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "gff2bed < {input} > {output}"


rule sort_peak_bed:
    input:
        "results/peaks/{dir}/{bg_sample}.bed"
    output:
        "results/peaks/{dir}/{bg_sample}.sorted.bed"
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/sort_peak_bed/{dir}_{bg_sample}.log"
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"


# create bed file of overlapping peaks between replicate conditions
rule overlapping_peaks: # escape bg_sample wildcard to get all replicate bg_samples
    input:
        beds=expand("results/peaks/{dir}/{{bg_sample}}.sorted.bed", dir=DIRS)
    output:
        "results/peaks/overlapping_peaks/{bg_sample}.overlap.bed"
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/overlapping_peaks/{bg_sample}.log"
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "bedtools multiinter {params.extra} -i {input.beds} > {output}"


rule filter_overlapping_peaks:
    input:
        bed="results/peaks/overlapping_peaks/{bg_sample}.overlap.bed",
        cs=f"resources/{resources.genome}_chrom.sizes",
    output:
        "results/peaks/overlapping_peaks/{bg_sample}.filtered.bed",
    params:
        m=config["peak_calling_perl"]["overlapping_peaks"]["max_size"],
        e=config["peak_calling_perl"]["overlapping_peaks"]["extend_by"],
        k=config["peak_calling_perl"]["overlapping_peaks"]["keep"],
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/filter_bed_file/{bg_sample}.log"
    conda:
        "../envs/peak_calling.yaml"
    script:
        workflow.source_path("../scripts/filter_overlapping_peaks.py")


rule annotate_peaks:
    input:
        bed="results/peaks/overlapping_peaks/{bg_sample}.filtered.bed",
        adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
        gtf=resources.gtf,
    output:
        txt=report("results/peaks/overlapping_peaks/{bg_sample}.annotated.txt",  caption="../report/annotated_peaks.rst", category="Annotated peaks"),
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/annotate_peaks/{bg_sample}.log"
    conda:
        "../envs/R.yaml"
    script:
        workflow.source_path("../scripts/annotate_peaks.R")
