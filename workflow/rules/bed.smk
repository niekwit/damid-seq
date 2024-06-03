rule gff2bed:
    input:
        "results/peaks/fdr{fdr}/{dir}/{bg_sample}.peaks.gff"
    output:
        "results/peaks/fdr{fdr}/{dir}/{bg_sample}.bed"
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/gff2bed/{dir}_{bg_sample}_fdr{fdr}.log"
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "gff2bed < {input} > {output}"


rule sort_peak_bed:
    input:
        "results/peaks/fdr{fdr}/{dir}/{bg_sample}.bed"
    output:
        "results/peaks/fdr{fdr}/{dir}/{bg_sample}.sorted.bed"
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/sort_peak_bed/{dir}_{bg_sample}_fdr{fdr}.log"
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "sort -k1,1 -k2,2n {input} > {output}"


# Create bed file of consensus peaks between replicate conditions
rule consensus_peaks: # Escape bg_sample wildcard to get all replicate bg_samples
    input:
        beds=expand("results/peaks/fdr{fdr}/{dir}/{{bg_sample}}.sorted.bed", fdr=fdr , dir=DIRS)
    output:
        "results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.bed"
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/consensus_peaks/fdr{fdr}/{bg_sample}.log"
    conda:
        "../envs/peak_calling.yaml"
    shell:
        "bedtools multiinter {params.extra} -i {input.beds} > {output}"


rule filter_consensus_peaks:
    input:
        bed="results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.bed",
        peaks=expand("results/peaks/fdr{fdr}/{dir}/{{bg_sample}}.sorted.bed", fdr=fdr , dir=DIRS),
        cs=f"resources/{resources.genome}_chrom.sizes",
    output:
        "results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.filtered.bed",
    params:
        k=config["consensus_peaks"]["keep"],
        max_size=config["consensus_peaks"]["max_size"],
        e=config["consensus_peaks"]["extend_by"],
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/filter_bed_file_fdr{fdr}/{bg_sample}.log"
    conda:
        "../envs/peak_calling.yaml"
    script:
        "../scripts/filter_consensus_peaks.py"


rule annotate_peaks:
    input:
        bed="results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.filtered.bed",
        adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
        gtf=resources.gtf,
    output:
        txt=report("results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.annotated.txt",  caption="../report/annotated_peaks.rst", category="Annotated peaks"),
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/annotate_peaks/fdr{fdr}/{bg_sample}.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/annotate_peaks.R"


rule get_gene_names:
    input:
        txt="results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.annotated.txt"
    output:
        ids="results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.geneIDs.txt"
    threads: 1
    resources:
        runtime=5
    log:
        "logs/peaks_fdr{fdr}_geneIDs/{bg_sample}.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "sed '1d' {input.txt} | "
        "awk '{{print $(NF-4),$(NF-1)}}' |"
        " sort | "
        "uniq > {output.ids}"