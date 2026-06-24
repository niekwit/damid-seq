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
        log:
            "logs/macs3_{mode}/fdr{fdr}/{dir}/{bg_sample}.log",
        conda:
            "../envs/peak_calling.yaml"
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        params:
            name=lambda w, input: os.path.basename(input["treatment"]).replace(
                ".sorted.bam", ""
            ),
            outdir=lambda w, output: os.path.dirname(output[0]),
            args=macs3_params(),  # Contains format, genome, cutoffs, mode and extra
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
        log:
            "logs/consensus_peaks_macs3/{mode}/fdr{fdr}/{bg_sample}.log",
        conda:
            "../envs/peak_calling.yaml"
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        params:
            extra="",
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
        log:
            "logs/filter_consensus_peaks/macs3_{mode}/fdr{fdr}/{bg_sample}.log",
        conda:
            "../envs/peak_calling.yaml"
        threads: 1
        resources:
            runtime=15,
        params:
            k=config["consensus_peaks"]["keep"],
            max_size=config["consensus_peaks"]["max_size"],
            e=config["consensus_peaks"]["extend_by"],
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
        log:
            "logs/plotting/macs3_{mode}_peak_annotation_plots/fdr{fdr}.log",
        conda:
            "../envs/R.yaml"
        threads: config["resources"]["plotting"]["cpu"]
        resources:
            runtime=config["resources"]["plotting"]["time"],
        params:
            extra="",
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
        log:
            "logs/annotate_peaks_macs3_{mode}/fdr{fdr}/{bg_sample}.log",
        conda:
            "../envs/R.yaml"
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
        params:
            extra="",
        script:
            "../scripts/annotate_peaks.R"

    rule get_gene_names_macs3:
        input:
            txt="results/macs3_{mode}/fdr{fdr}/{bg_sample}.annotated.txt",
        output:
            ids="results/macs3_{mode}/fdr{fdr}/{bg_sample}.geneIDs.txt",
        log:
            "logs/geneIDs_peaks_macs3_{mode}_fdr{fdr}/{bg_sample}.log",
        conda:
            "../envs/deeptools.yaml"
        threads: 1
        resources:
            runtime=5,
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
                "results/macs3_{mode}/fdr{fdr}/{dir}/{bg_sample}_peaks.{mode}Peak",
                fdr=fdr,
                dir=DIRS,
                bg_sample=BG_SAMPLES,
                mode=PEAK_MODE,
            ),
        output:
            csv="results/macs3_{mode}/fdr{fdr}/read_counts/reads_in_peaks.csv",
        log:
            "logs/frip/{mode}_fdr{fdr}.log",
        conda:
            "../envs/deeptools.yaml"
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"],
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
        log:
            "logs/plot_frip/{mode}/fdr{fdr}.log",
        conda:
            "../envs/R.yaml"
        threads: config["resources"]["plotting"]["cpu"]
        resources:
            runtime=config["resources"]["plotting"]["time"],
        params:
            extra="",
        script:
            "../scripts/plot_frip.R"

    if config["consensus_peaks"]["enrichment_analysis"]["run"]:

        rule enrichment_analysis_macs3:
            input:
                txt="results/macs3_{mode}/fdr{fdr}/{bg_sample}.geneIDs.txt",
            output:
                xlsx="results/macs3_{mode}/fdr{fdr}/enrichment_analysis/{bg_sample}.xlsx",
            log:
                "logs/enrichment_analysis/macs3_{mode}/fdr{fdr}/{bg_sample}.log",
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
            log:
                "logs/plot_enrichment/macs3_{mode}/fdr{fdr}/{bg_sample}.log",
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
