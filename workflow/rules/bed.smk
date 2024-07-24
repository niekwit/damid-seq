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
        in_file="results/peaks/fdr{fdr}/{dir}/{bg_sample}.bed",
        genome=f"resources/{resources.genome}_chrom_order.txt",
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
    wrapper:
        f"{wrapper_version}/bio/bedtools/sort"


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

if config["peak_calling_perl"]["run"]:
    if config["consensus_peaks"]["enrichment_analysis"]["run"]:
        rule enrichment_analysis:
            input:
                txt="results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.geneIDs.txt",
            output:
                xlsx="results/peaks/fdr{fdr}/consensus_peaks/enrichment_analysis/{bg_sample}.xlsx",
            params:
                extra="",
                genome=resources.genome,
                dbs=config["consensus_peaks"]["enrichment_analysis"]["dbs"]
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/enrichment_analysis/peaks/fdr{fdr}/{bg_sample}.log"
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/enrichment_analysis.R"

        rule plot_enrichment:
            input:
                xlsx="results/peaks/fdr{fdr}/consensus_peaks/enrichment_analysis/{bg_sample}.xlsx",
            output:
                plots=report(expand("results/plots/peaks/fdr{{fdr}}/enrichment_analysis/{{bg_sample}}/{db}.pdf", db=DBS), caption="../report/enrichment_analysis.rst", category="Pathway enrichment analysis"),
            params:
                terms=config["consensus_peaks"]["enrichment_analysis"]["terms"],
                dirname=lambda w, output: os.path.dirname(output[0]),
            threads: config["resources"]["plotting"]["cpu"]
            resources:
                runtime=config["resources"]["plotting"]["time"]
            log:
                "logs/plot_enrichment/peaks/fdr{fdr}/{bg_sample}.log"
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/plot_enrichment.R"

    rule count_reads_in_peaks:
        # Adapted from https://www.biostars.org/p/337872/#337890
        input:
            bam="results/bam/{dir}/{bg_sample}.extended.sorted.bam",
            bai="results/bam/{dir}/{bg_sample}.extended.sorted.bam.bai",
            b="results/peaks/fdr{fdr}/{dir}/{bg_sample}.sorted.bed",
            order=f"resources/{resources.genome}_chrom_order.txt",
        output:
            total_read_count="results/peaks/fdr{fdr}/read_counts/{dir}/{bg_sample}.total.count",
            peak_read_count="results/peaks/fdr{fdr}/read_counts/{dir}/{bg_sample}.peak.count",
        params:
            extra="",
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        log:
            "logs/bedtools_intersect/fdr{fdr}/{dir}/{bg_sample}.log"
        conda:
            "../envs/peak_calling.yaml"
        shell:
            "bedtools bamtobed "
            "{params.extra} "
            "-i {input.bam} | "
            "bedtools sort -g {input.order} "
            "tee >(wc -l > {output.total_read_count}) | "
            "bedtools intersect "
            "{params.extra} "
            "-sorted "
            "-g {input.order} "
            "-c "
            "-a {input.b} "
            "-b stdin | "
            "awk '{{i+=$NF}}END{{print i}}' > "
            "{output.peak_read_count} "
                

    rule plot_fraction_of_reads_in_peaks:
        input:
            total_read_count=expand("results/peaks/fdr{fdr}/read_counts/{dir}/{bg_sample}.total.count", dir=DIRS, fdr=fdr, bg_sample=BG_SAMPLES),
            peak_read_count=expand("results/peaks/fdr{fdr}/read_counts/{dir}/{bg_sample}.peak.count", dir=DIRS, fdr=fdr, bg_sample=BG_SAMPLES),
        output:
            plot=report("results/plots/peaks/fdr{fdr}/frip.pdf", caption="../report/frip.rst", category="Fraction of reads in peaks"),
            csv="results/peaks/fdr{fdr}/frip.csv",
        params:
            extra="",
        threads: config["resources"]["plotting"]["cpu"]
        resources:
            runtime=config["resources"]["plotting"]["time"]
        log:
            "logs/plot_frip/fdr{fdr}.log"
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/plot_frip.R"