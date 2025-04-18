rule plotPCA:
    input:
        "results/deeptools/PCA.tab",
    output:
        pca=report("results/plots/PCA.pdf", caption="../report/pca.rst", category="PCA"),
        scree=report(
            "results/plots/scree.pdf", caption="../report/scree.rst", category="PCA"
        ),
    params:
        extra="",
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"],
    log:
        "logs/plotting/plotPCA.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_PCA.R"


rule plot_correlation:
    input:
        "results/deeptools/scores_per_bin.npz",
    output:
        tab="results/deeptools/correlation.tab",
        pdf=report(
            "results/plots/sample_correlation.pdf",
            caption="../report/correlation.rst",
            category="Sample correlation",
        ),
    params:
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"],
    log:
        "logs/plotting/plotCorrelation.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotCorrelation "
        "--corData {input} "
        "--corMethod spearman "
        "--whatToPlot heatmap "
        "--removeOutliers "
        "--plotFile {output.pdf} "
        "--outFileCorMatrix {output.tab} "
        "--colorMap viridis "
        "--skipZeros "
        "{params.extra} "
        "> {log} 2>&1"


rule plot_heatmap:
    input:
        mat="results/deeptools/average_bw_matrix.gz",
    output:
        pdf=report(
            "results/plots/heatmap.pdf",
            caption="../report/heatmap.rst",
            category="Heatmap",
        ),
        mat="results/deeptools/heatmap_matrix.gz",
    params:
        im=config["deeptools"]["plotHeatmap"]["interpolationMethod"],
        pt=config["deeptools"]["plotHeatmap"]["plotType"],
        cm=config["deeptools"]["plotHeatmap"]["colorMap"],
        a=config["deeptools"]["plotHeatmap"]["alpha"],
        extra=config["deeptools"]["plotHeatmap"]["extra"],
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"],
    log:
        "logs/deeptools/plotHeatmap.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotHeatmap "
        "--matrixFile {input.mat} "
        "--outFileNameMatrix {output.mat} "
        "--outFileName {output.pdf} "
        "--colorMap {params.cm} "
        "--alpha {params.a} "
        "{params.extra} "
        "> {log} 2>&1"


rule plot_profile:
    input:
        mat="results/deeptools/average_bw_matrix.gz",
    output:
        pdf=report(
            "results/plots/profile_plot.pdf",
            caption="../report/profile_plot.rst",
            category="Profile plot",
        ),
    params:
        rl=computematrix_args(region_labels=True),
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"],
    log:
        "logs/deeptools/plotProfile.log",
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotProfile "
        "--matrixFile {input.mat} "
        "--outFileName {output.pdf} "
        "--perGroup "
        "--plotType=lines "
        "--legendLocation=upper-right "
        "{params.extra} "
        "> {log} 2>&1"


if config["peak_calling_perl"]["run"]:
    fdr = config["peak_calling_perl"]["fdr"]

    rule peak_annotation_plots:
        input:
            gtf=resources.gtf,
            bed=expand(
                "results/peaks/fdr{fdr}/consensus_peaks/{bg_sample}.filtered.bed",
                fdr=fdr,
                bg_sample=BG_SAMPLES,
            ),
        output:
            fd=report(
                "results/plots/peaks/fdr{fdr}/feature_distributions.pdf",
                caption="../report/feature_distributions.rst",
                category="Peak annotation",
            ),
            dt=report(
                "results/plots/peaks/fdr{fdr}/distance_to_tss.pdf",
                caption="../report/distance_to_tss.rst",
                category="Peak annotation",
            ),
        params:
            extra="",
        threads: config["resources"]["plotting"]["cpu"]
        resources:
            runtime=config["resources"]["plotting"]["time"],
        log:
            "logs/plotting/peak_annotation_plots_fdr{fdr}.log",
        conda:
            "../envs/R.yaml"
        script:
            "../scripts/peak_annotation_plots.R"

elif config["peak_calling_macs3"]["run"]:
    if config["peak_calling_macs3"]["mode"] == "narrow":
        fdr = config["peak_calling_macs3"]["qvalue"]

        rule peak_annotation_plots:
            input:
                gtf=resources.gtf,
                bed=expand(
                    "results/macs2_narrow/fdr{fdr}/consensus_peaks/{bg_sample}.filtered.bed",
                    fdr=fdr,
                    bg_sample=BG_SAMPLES,
                ),
            output:
                fd=report(
                    "results/plots/macs2_narrow/fdr{fdr}/feature_distributions.pdf",
                    caption="../report/feature_distributions.rst",
                    category="Peak annotation",
                ),
                dt=report(
                    "results/plots/macs2_narrow/fdr{fdr}/distance_to_tss.pdf",
                    caption="../report/distance_to_tss.rst",
                    category="Peak annotation",
                ),
            params:
                extra="",
            threads: config["resources"]["plotting"]["cpu"]
            resources:
                runtime=config["resources"]["plotting"]["time"],
            log:
                "logs/plotting/peak_annotation_plots_fdr{fdr}_macs2_narrow.log",
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/peak_annotation_plots.R"

    elif config["peak_calling_macs3"]["mode"] == "broad":
        fdr = config["peak_calling_macs3"]["broad_cutoff"]

        rule peak_annotation_plots:
            input:
                gtf=resources.gtf,
                bed=expand(
                    "results/macs2_broad/fdr{fdr}/consensus_peaks/{bg_sample}.filtered.bed",
                    fdr=fdr,
                    bg_sample=BG_SAMPLES,
                ),
            output:
                fd=report(
                    "results/plots/macs2_broad/fdr{fdr}/feature_distributions.pdf",
                    caption="../report/feature_distributions.rst",
                    category="Peak annotation",
                ),
                dt=report(
                    "results/plots/macs2_broad/fdr{fdr}/distance_to_tss.pdf",
                    caption="../report/distance_to_tss.rst",
                    category="Peak annotation",
                ),
            params:
                extra="",
            threads: config["resources"]["plotting"]["cpu"]
            resources:
                runtime=config["resources"]["plotting"]["time"],
            log:
                "logs/plotting/peak_annotation_plots_fdr{fdr}_macs2_broad.log",
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/peak_annotation_plots.R"


rule plot_mapping_rates:
    input:
        log=expand("logs/bowtie2_align/{dir}/{sample}.log", dir=DIRS, sample=SAMPLES),
    output:
        pdf=report(
            "results/plots/mapping_rates.pdf",
            caption="../report/mapping_rates.rst",
            category="Mapping rates",
        ),
    params:
        extra="",
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"],
    log:
        "logs/plotting/plot_mapping_rates.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_mapping_rates.R"
