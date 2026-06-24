if config["peak_calling_perl"]["run"]:
    fdr = config["peak_calling_perl"]["fdr"]

    rule peak_calling_perl:
        input:
            fp="resources/find_peaks",
            bg="results/bedgraph/{dir}/{bg_sample}-vs-Dam-norm.gatc.bedgraph",
        output:
            gff="results/peaks/fdr{fdr}/{dir}/{bg_sample}.peaks.gff",
            data="results/peaks/fdr{fdr}/{dir}/{bg_sample}.data",
        log:
            "logs/find_peaks/fdr{fdr}/{dir}/{bg_sample}.log",
        conda:
            "../envs/peak_calling.yaml"
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        params:
            outdir=lambda w, output: os.path.dirname(output["gff"]),
            n=config["peak_calling_perl"]["iterations"],
            fdr=config["peak_calling_perl"]["fdr"],
            frac=config["peak_calling_perl"]["fraction"],
            mc=config["peak_calling_perl"]["min_count"],
            mq=config["peak_calling_perl"]["min_quantile"],
            step=config["peak_calling_perl"]["step"],
            up=config["peak_calling_perl"]["unified_peaks"],
        script:
            "../scripts/run_find_peaks.py"

    rule gff2bed:
        input:
            "results/peaks/fdr{fdr}/{dir}/{bg_sample}.peaks.gff",
        output:
            "results/peaks/fdr{fdr}/{dir}/{bg_sample}.bed",
        log:
            "logs/gff2bed/{dir}_{bg_sample}_fdr{fdr}.log",
        conda:
            "../envs/peak_calling.yaml"
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        params:
            extra="",
        shell:
            "gff2bed < {input} > {output}"

    rule sort_peak_bed:
        input:
            in_file="results/peaks/fdr{fdr}/{dir}/{bg_sample}.bed",
            genome=f"resources/{resources.genome}_chrom_order.txt",
        output:
            "results/peaks/fdr{fdr}/{dir}/{bg_sample}.sorted.bed",
        log:
            "logs/sort_peak_bed/{dir}_{bg_sample}_fdr{fdr}.log",
        conda:
            "../envs/peak_calling.yaml"
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        params:
            extra="",
        wrapper:
            "v5.8.3/bio/bedtools/sort"

    # Create bed file of consensus peaks between replicate conditions
    rule consensus_peaks:  # Escape bg_sample wildcard to get all replicate bg_samples
        input:
            beds=expand(
                "results/peaks/fdr{fdr}/{dir}/{{bg_sample}}.sorted.bed",
                fdr=fdr,
                dir=DIRS,
            ),
        output:
            "results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.bed",
        log:
            "logs/consensus_peaks/fdr{fdr}/{bg_sample}.log",
        conda:
            "../envs/peak_calling.yaml"
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        params:
            extra="",
        shell:
            "bedtools multiinter {params.extra} -i {input.beds} > {output}"

    rule filter_consensus_peaks:
        input:
            bed="results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.bed",
            peaks=expand(
                "results/peaks/fdr{fdr}/{dir}/{{bg_sample}}.sorted.bed",
                fdr=fdr,
                dir=DIRS,
            ),
            cs=f"resources/{resources.genome}_chrom.sizes",
        output:
            "results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.filtered.bed",
        log:
            "logs/filter_bed_file_fdr{fdr}/{bg_sample}.log",
        conda:
            "../envs/peak_calling.yaml"
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        params:
            k=config["consensus_peaks"]["keep"],
            max_size=config["consensus_peaks"]["max_size"],
            e=config["consensus_peaks"]["extend_by"],
            extra="",
        script:
            "../scripts/filter_consensus_peaks.py"

    rule annotate_peaks:
        input:
            bed="results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.filtered.bed",
            adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
            gtf=resources.gtf,
        output:
            txt=report(
                "results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.annotated.txt",
                caption="../report/annotated_peaks.rst",
                category="Annotated peaks",
            ),
        log:
            "logs/annotate_peaks/fdr{fdr}/{bg_sample}.log",
        conda:
            "../envs/R.yaml"
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        params:
            extra="",
        script:
            "../scripts/annotate_peaks.R"

    rule get_gene_names:
        input:
            txt="results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.annotated.txt",
        output:
            ids="results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.geneIDs.txt",
        log:
            "logs/peaks_fdr{fdr}_geneIDs/{bg_sample}.log",
        conda:
            "../envs/deeptools.yaml"
        threads: 1
        resources:
            runtime=5,
        shell:
            "sed '1d' {input.txt} | "
            "awk '{{print $(NF-4),$(NF-1)}}' |"
            " sort | "
            "uniq > {output.ids}"

    if config["consensus_peaks"]["enrichment_analysis"]["run"]:

        rule enrichment_analysis:
            input:
                txt="results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.geneIDs.txt",
            output:
                xlsx="results/peaks/fdr{fdr}/consensus_peaks/enrichment_analysis/{bg_sample}.xlsx",
            log:
                "logs/enrichment_analysis/peaks/fdr{fdr}/{bg_sample}.log",
            conda:
                "../envs/R.yaml"
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"],
            params:
                extra="",
                genome=resources.genome,
                dbs=config["consensus_peaks"]["enrichment_analysis"]["dbs"],
            script:
                "../scripts/enrichment_analysis.R"

        rule plot_enrichment:
            input:
                xlsx="results/peaks/fdr{fdr}/consensus_peaks/enrichment_analysis/{bg_sample}.xlsx",
            output:
                plots=report(
                    expand(
                        "results/plots/peaks/fdr{{fdr}}/enrichment_analysis/{{bg_sample}}/{db}.pdf",
                        db=DBS,
                    ),
                    caption="../report/enrichment_analysis.rst",
                    category="Pathway enrichment analysis",
                ),
            log:
                "logs/plot_enrichment/peaks/fdr{fdr}/{bg_sample}.log",
            conda:
                "../envs/R.yaml"
            threads: config["resources"]["plotting"]["cpu"]
            resources:
                runtime=config["resources"]["plotting"]["time"],
            params:
                terms=config["consensus_peaks"]["enrichment_analysis"]["terms"],
                dirname=lambda w, output: os.path.dirname(output[0]),
            script:
                "../scripts/plot_enrichment.R"

    rule count_reads_in_peaks:
        # Adapted from https://www.biostars.org/p/337872/#337890
        input:
            bams=expand(
                "results/bam/{dir}/{bg_sample}.sorted.bam",
                dir=DIRS,
                bg_sample=BG_SAMPLES,
            ),
            bai=expand(
                "results/bam/{dir}/{bg_sample}.sorted.bam.bai",
                dir=DIRS,
                bg_sample=BG_SAMPLES,
            ),
            peaks=expand(
                "results/peaks/fdr{fdr}/{dir}/{bg_sample}.sorted.bed",
                fdr=fdr,
                dir=DIRS,
                bg_sample=BG_SAMPLES,
            ),
        output:
            csv="results/peaks/fdr{fdr}/read_counts/reads_in_peaks.csv",
        log:
            "logs/frip/perl_fdr{fdr}.log",
        conda:
            "../envs/deeptools.yaml"
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        params:
            extra="",
        script:
            "../scripts/count_reads_in_peaks.py"

    rule plot_fraction_of_reads_in_peaks:
        input:
            csv="results/peaks/fdr{fdr}/read_counts/reads_in_peaks.csv",
        output:
            plot=report(
                "results/plots/peaks/fdr{fdr}/frip.pdf",
                caption="../report/frip.rst",
                category="Fraction of reads in peaks",
            ),
        log:
            "logs/plot_frip/fdr{fdr}.log",
        conda:
            "../envs/R.yaml"
        threads: config["resources"]["plotting"]["cpu"]
        resources:
            runtime=config["resources"]["plotting"]["time"],
        params:
            extra="",
        script:
            "../scripts/plot_frip.R"
