rule bedgraph2bigwig:
    input:
        cs=f"resources/{resources.genome}_chrom.sizes",
        bg="results/bedgraph/{dir}/{sample}-vs-Dam.gatc.bedgraph" # exclude dam sample from sample wildcard here!!!!
    output:
        bw="results/bigwig/{dir}/{sample}.bw"
    params:
        extra=""
    threads: config["resources"]["fastqc"]["cpu"]
    resources: 
        runtime=config["resources"]["fastqc"]["time"],
    conda: 
        "../envs/deeptools.yaml"
    shell:
        "bedGraphToBigWig "
        "{params.extra} "
        "{input.bg} "
        "{input.cs} "
        "{output} > {log} 2>&1"


rule average_bigwigs:
    input:
        expand("results/bigwig/{dir}/{sample}.bw", dir=DIRS, sample=SAMPLES),
    output:
        bw="results/bigwig/average_bw/{sample}.bw",
    params:
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/deeptools/bw_average_{condition}.log"
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/average_bigwig.py"