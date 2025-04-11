if config["peak_calling_macs3"]["run"]:
    fdr = get_fdr()

    rule peak_calling_macs3:
        input:
            treatment=f"results/bam/{{dir}}/{{bg_sample}}.sorted.bam",
            control=f"results/bam/{{dir}}/Dam.sorted.bam",
        output:
            multiext(
                "results/macs3_{mode}/fdr{fdr}/{dir}/{bg_sample}",
                "_peaks.xls",
                "_peaks.{mode}Peak",
                "_cutoff_analysis.txt",
            ),
        params:
            name=lambda w, input: os.path.basename(input["treatment"]).replace(
                ".sorted.bam", ""
            ),
            outdir=lambda w, output: os.path.dirname(output[0]),
            args=macs3_params(),  # Contains format, genome, cutoffs, mode and extra
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        log:
            "logs/macs3_{mode}/fdr{fdr}/{dir}/{bg_sample}.log",
        conda:
            "../envs/peak_calling.yaml"
        shell:
            "macs3 callpeak "
            "--treatment {input.treatment} "
            "--control {input.control} "
            "--name {params.name} "
            "--outdir {params.outdir} "
            "--keep-dup all "
            "--cutoff-analysis "
            "{params.args} "
            "2> {log}"

    rule consensus_peaks_macs3:
        input:
            peaks=expand(
                "results/macs3_{mode}/fdr{fdr}/{dir}/{{bg_sample}}_peaks.{mode}Peak",
                dir=DIRS,
                fdr=fdr,
                mode=PEAK_MODE,
            ),
        output:
            "results/macs3_{mode}/fdr{fdr}/consensus_peaks/{bg_sample}.intersect.bed",
        params:
            extra="",
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        log:
            "logs/consensus_peaks_macs3/{mode}/fdr{fdr}/{bg_sample}.log",
        conda:
            "../envs/peak_calling.yaml"
        shell:
            "bedtools multiinter {params.extra} -i {input.peaks} > {output}"

    rule filter_consensus_peaks_macs3:
        input:
            bed="results/macs3_{mode}/fdr{fdr}/consensus_peaks/{bg_sample}.intersect.bed",
            peaks=expand(
                "results/macs3_{mode}/fdr{fdr}/{dir}/{{bg_sample}}_peaks.{mode}Peak",
                mode=PEAK_MODE,
                dir=DIRS,
                fdr=fdr,
            ),
            cs=f"resources/{resources.genome}_chrom.sizes",
        output:
            ext_bed="results/macs3_{mode}/fdr{fdr}/consensus_peaks/{bg_sample}.intersect.filtered.bed",
        params:
            k=config["consensus_peaks"]["keep"],
            max_size=config["consensus_peaks"]["max_size"],
            e=config["consensus_peaks"]["extend_by"],
        threads: 1
        resources:
            runtime=15,
        log:
            "logs/filter_consensus_peaks/macs3_{mode}/fdr{fdr}/{bg_sample}.log",
        conda:
            "../envs/peak_calling.yaml"
        script:
            "../scripts/filter_consensus_peaks.py"

    rule peak_annotation_plots_macs3:
        input:
            gtf=resources.gtf,
            bed=expand(
                "results/macs3_{mode}/fdr{fdr}/consensus_peaks/{bg_sample}.intersect.filtered.bed",
                mode=PEAK_MODE,
                fdr=fdr,
                bg_sample=BG_SAMPLES,
            ),
        output:
            fd=report(
                "results/plots/macs3_{mode}/fdr{fdr}/feature_distributions.pdf",
                caption="../report/feature_distributions.rst",
                category="Peak annotation",
            ),
            dt=report(
                "results/plots/macs3_{mode}/fdr{fdr}/distance_to_tss.pdf",
                caption="../report/distance_to_tss.rst",
                category="Peak annotation",
            ),
        params:
            extra="",
        threads: config["resources"]["plotting"]["cpu"]
        resources:
            runtime=config["resources"]["plotting"]["time"],
        log:
            "logs/plotting/macs3_{mode}_peak_annotation_plots/fdr{fdr}.log",
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/peak_annotation_plots.R"

    rule annotate_peaks_macs3:
        input:
            bed="results/macs3_{mode}/fdr{fdr}/consensus_peaks/{bg_sample}.intersect.filtered.bed",
            adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
            gtf=resources.gtf,
        output:
            txt=report(
                "results/macs3_{mode}/fdr{fdr}/{bg_sample}.annotated.txt",
                caption="../report/annotated_peaks.rst",
                category="Annotated peaks",
            ),
        params:
            extra="",
        log:
            "logs/annotate_peaks_macs3_{mode}/fdr{fdr}/{bg_sample}.log",
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/annotate_peaks.R"

    rule get_gene_names_macs3:
        input:
            txt="results/macs3_{mode}/fdr{fdr}/{bg_sample}.annotated.txt",
        output:
            ids="results/macs3_{mode}/fdr{fdr}/{bg_sample}.geneIDs.txt",
        threads: 1
        resources:
            runtime=5,
        log:
            "logs/geneIDs_peaks_macs3_{mode}_fdr{fdr}/{bg_sample}.log",
        conda:
            "../envs/deeptools.yaml"
        shell:
            "sed '1d' {input.txt} | "
            "awk '{{print $(NF-4),$(NF-1)}}' | "
            "sort | "
            "uniq > {output.ids}"

    rule count_reads_in_peaks_macs3:
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
            csv="results/macs3_{mode}/fdr{fdr}/read_counts/reads_in_peaks.csv",
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        log:
            "logs/frip/{mode}_fdr{fdr}.log",
        conda:
            "../envs/deeptools.yaml"
        script:
            "../scripts/count_reads_in_peaks.py"

    rule plot_fraction_of_reads_in_peaks_macs3:
        input:
            csv="results/macs3_{mode}/fdr{fdr}/read_counts/reads_in_peaks.csv",
        output:
            plot=report(
                "results/plots/macs3_{mode}/fdr{fdr}/frip.pdf",
                caption="../report/frip.rst",
                category="Fraction of reads in peaks",
            ),
        params:
            extra="",
        threads: config["resources"]["plotting"]["cpu"]
        resources:
            runtime=config["resources"]["plotting"]["time"],
        log:
            "logs/plot_frip/{mode}/fdr{fdr}.log",
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/plot_frip.R"

    if config["consensus_peaks"]["enrichment_analysis"]["run"]:

        rule enrichment_analysis_macs3:
            input:
                txt="results/macs3_{mode}/fdr{fdr}/{bg_sample}.geneIDs.txt",
            output:
                xlsx="results/macs3_{mode}/fdr{fdr}/enrichment_analysis/{bg_sample}.xlsx",
            params:
                extra="",
                genome=resources.genome,
                dbs=config["consensus_peaks"]["enrichment_analysis"]["dbs"],
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"],
            log:
                "logs/enrichment_analysis/macs3_{mode}/fdr{fdr}/{bg_sample}.log",
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/enrichment_analysis.R"

        rule plot_enrichment_macs3:
            input:
                xlsx="results/macs3_{mode}/fdr{fdr}/enrichment_analysis/{bg_sample}.xlsx",
            output:
                plots=report(
                    expand(
                        "results/plots/macs3_{{mode}}/fdr{{fdr}}/enrichment_analysis/{{bg_sample}}/{db}.pdf",
                        db=DBS,
                    ),
                    caption="../report/enrichment_analysis.rst",
                    category="Pathway enrichment analysis",
                ),
            params:
                terms=config["consensus_peaks"]["enrichment_analysis"]["terms"],
                dirname=lambda w, output: os.path.dirname(output[0]),
            threads: config["resources"]["plotting"]["cpu"]
            resources:
                runtime=config["resources"]["plotting"]["time"],
            log:
                "logs/plot_enrichment/macs3_{mode}/fdr{fdr}/{bg_sample}.log",
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/plot_enrichment.R"
