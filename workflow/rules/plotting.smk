rule plotPCA:
    input:
        "results/deeptools/PCA.tab",
    output:
        pca="results/plots/PCA.pdf",
        scree="results/plots/scree.pdf",
    params:
        extra=""
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plotting/plotPCA.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_PCA.R"


rule plotCorrelation:
    input:
        "results/deeptools/scores_per_bin.npz"
    output:
        tab="results/deeptools/correlation.tab",
        pdf="results/plots/sample_correlation.pdf",
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/plotting/plotCorrelation.log"
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
        "> {log} 2>&1"


rule plot_heatmap:
    input:
        mat="results/deeptools/average_bw_matrix.gz",
    output:
        pdf="results/plots/heatmap.pdf",
        mat="results/deeptools/heatmap_matrix.gz",
    params:
        im=config["deeptools"]["plotHeatmap"]["interpolationMethod"],
        pt=config["deeptools"]["plotHeatmap"]["plotType"],
        cm=config["deeptools"]["plotHeatmap"]["colorMap"],
        a=config["deeptools"]["plotHeatmap"]["alpha"],
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/deeptools/plotHeatmap.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotHeatmap "
        "--matrixFile {input.mat} "
        "--outFileNameMatrix {output.mat} "
        "--outFileName {output.pdf} "
        "--perGroup "
        "--colorMap {params.cm} "
        "--alpha {params.a} "
        "{params.extra} "
        "> {log} 2>&1"


rule plot_profile:
    input:
        mat="results/deeptools/average_bw_matrix.gz",
    output:
        pdf="results/plots/profile_plot.pdf",
    params:
        rl = computematrix_args(region_labels=True),
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/deeptools/plotProfile.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotProfile "
        "--matrixFile {input.mat} "
        "--outFileName {output.pdf} "
        "--perGroup "
        "--plotType=fill "
        "--legendLocation=upper-right "
        "{params.extra} "
        "> {log} 2>&1"


rule peak_annotation_plots:
    input:
        gtf=resources.gtf,
        bed=expand("results/peaks/overlapping_peaks/{bg_sample}.extended.bed", bg_sample=BG_SAMPLES),
    output:
        fd="results/plots/peaks/feature_distributions.pdf",
        dt="results/plots/peaks/distance_to_tss.pdf",
        pa="results/plots/peaks/pathway_enrichment.pdf",
        v="results/plots/peaks/venn_overlap_conditions.pdf",
    params:
        extra="",
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plotting/peak_annotation_plots.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/peak_annotation_plots.R"

