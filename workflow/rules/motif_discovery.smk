# create bed file of overlapping peaks between replicate conditions
rule overlapping_peaks:
    input:
        gff=expand("results/peaks/{dir}/{bg_sample}.gff", dir=DIRS, bg_sample=BG_SAMPLES),
    output:
        "results/peaks/overlapping_peaks/{bg_sample}.bed"
    params:
        min_overlap=config["peak_calling"]["overlapping_peaks"]["min_overlap"],
        extra=""
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/overlapping_peaks/{bg_sample}.log"
    conda:
        "../envs/peak_calling.yaml"
    script:
        "../scripts/overlapping_peaks.py"

        