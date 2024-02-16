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


rule remove_plasmid_loci:
    input:
        bed="results/peaks/{dir}/{bg_sample}.sorted.bed",
        bl="resources/blacklist.bed",
        txt="resources/genes.txt",
    output:
        bed="results/peaks/{dir}/{bg_sample}.no_plasmid.bed"
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/remove_plasmid_loci/{dir}_{bg_sample}.log"
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "bedtools intersect -v -a {input.bed} -b {input.bl} > {output.bed}"


# create bed file of overlapping peaks between replicate conditions
rule overlapping_peaks: # escape bg_sample wildcard to get all replicate bg_samples
    input:
        beds=expand("results/peaks/{dir}/{{bg_sample}}.no_plasmid.bed", dir=DIRS)
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
        m=config["peak_calling"]["overlapping_peaks"]["max_size"],
        e=config["peak_calling"]["overlapping_peaks"]["extend_by"],
        k=config["peak_calling"]["overlapping_peaks"]["keep"],
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/extend_bed_regions/{bg_sample}.log"
    conda:
        "../envs/peak_calling.yaml"
    script:
        "../scripts/filter_overlapping_peaks.py"


rule annotate_peaks:
    input:
        bed="results/peaks/overlapping_peaks/{bg_sample}.filtered.bed",
        adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
        gtf=resources.gtf,
    output:
        txt="results/peaks/overlapping_peaks/{bg_sample}.annotated.txt",
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
        "../scripts/annotate_peaks.R"


