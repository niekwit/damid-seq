rule peak_calling:
    input:
        fp="resources/find_peaks",
        bg="results/bedgraph/{dir}/{bg_sample}-vs-Dam.kde-norm.gatc.bedgraph",
    output:
        gff="results/peaks/{dir}/{bg_sample}.peaks.gff",
        data="results/peaks/{dir}/{bg_sample}.data",
    params:
        out=lambda w, output: os.path.dirname(output["gff"]),
        #seed=config["peak_calling"]["seed"],
        n=config["peak_calling"]["iterations"],
        fdr=config["peak_calling"]["fdr"],
        frac=config["peak_calling"]["fraction"],
        mc=config["peak_calling"]["min_count"],
        mq=config["peak_calling"]["min_quantile"],
        step=config["peak_calling"]["step"],
        up=config["peak_calling"]["unified_peaks"],
        #extra=config["peak_calling"]["extra"],
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    conda:
        "../envs/peak_calling.yaml"
    log:
        "logs/find_peaks/{dir}/{bg_sample}.log"
    script:
        "../scripts/run_find_peaks.py"
        
'''
    shell:
        "python {input.fp}/find_peaks.py "
        "{params.extra} "
        "--n {params.it} "
        "--fdr {params.fdr} "
        "--frac {params.frac} "
        "--min_count {params.mc} "
        "--min_quant {params.mq} "
        "--step {params.step} "
        "--unified_peaks {params.up} "
        "--seed {params.seed} "
        "--outdir {params.out} "
        "{input.bg} "
        "> {log} 2>&1 "
'''