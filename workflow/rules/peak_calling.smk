rule peak_calling:
    input:
        fp="resources/find_peaks_py",
        bg="results/bedgraph/{dir}/{bg_sample}-vs-Dam.kde-norm.gatc.bedgraph",
    output:
        gff="results/peaks/{dir}/{bg_sample}.gff",
    params:
        dir=lambda w, input: os.path.dirname(input),
        seed=config["peak_calling"]["seed"],
        it=config["peak_calling"]["iterations"],
        fdr=config["peak_calling"]["fdr"],
        extra=config["peak_calling"]["extra"],
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    conda:
        "../envs/damid.yaml"
    log:
        "logs/find_peaks_py/{dir}/{bg_sample}.log"
    shell:
        "python {input.fp}/find_peaks.py "
        "{params.extra} "
        "--n {params.it} "
        "--fdr {params.fdr} "
        "--seed {params.seed} "
        "--outdir {params.dir} "
        "{input.bg} "
        "> {log} 2>&1 "