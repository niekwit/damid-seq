rule multiBigwigSummary:
    input:
        expand("results/bigwig/{dir}/{bg_sample}.bw", dir=DIRS , bg_sample=BG_SAMPLES),
    output:
        "results/deeptools/scores_per_bin.npz",
    params:
        labels=lambda wildcards, input: [x.replace("results/bigwig/", "").replace(".bw","") for x in input],
        extra=""
    threads: config["resources"]["deeptools"]["cpu"] * 4
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/deeptools/multiBigwigSummary.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "multiBigwigSummary bins "
        "--bwfiles {input} "
        "--labels {params.labels} "
        "--outFileName {output} "
        "--numberOfProcessors {threads} "
        "{params.extra} "
        "> {log} 2>&1"


rule PCA:
    input:
        "results/deeptools/scores_per_bin.npz",
    output:
        "results/deeptools/PCA.tab",
    params:
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/deeptools/PCA.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotPCA "
        "--corData {input} "
        "--outFileNameData {output} "
        "--transpose "
        "{params.extra} "
        "> {log} 2>&1"


rule computeMatrix:
    input:
        bw=expand("results/bigwig/average_bw/{bg_sample}.bw", bg_sample=BG_SAMPLES),
        gtf=resources.gtf,
    output:
        mat="results/deeptools/average_bw_matrix.gz",
    params:
        args=computematrix_args(),
    threads: config["resources"]["deeptools"]["cpu"] * 4 # Otherwise it will take very long
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/deeptools/computeMatrix.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "computeMatrix "
        "{params.args} "
        "--numberOfProcessors {threads} "
        "--smartLabels "
        "--scoreFileName {input.bw} "
        "--outFileName {output.mat} "
        "> {log} 2>&1"