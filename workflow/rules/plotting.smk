rule plotPCA:
    input:
        "results/deeptools/PCA.tab",
    output:
        pca="results/deeptools/PCA.pdf",
        scree="results/deeptools/scree.pdf",
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
        pdf="results/deeptools/correlation.pdf"
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/plotting/plotCorrelation.log"
    conda:
        "../envs/R.yaml"
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
        im = config["deeptools"]["plotHeatmap"]["interpolationMethod"],
        pt = config["deeptools"]["plotHeatmap"]["plotType"],
        cm = config["deeptools"]["plotHeatmap"]["colorMap"],
        a = config["deeptools"]["plotHeatmap"]["alpha"],
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