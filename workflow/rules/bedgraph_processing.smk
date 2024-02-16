rule bedgraph2bigwig:
    input:
        cs=f"resources/{resources.genome}_chrom.sizes",
        bg="results/bedgraph/{dir}/{bg_sample}-vs-Dam.kde-norm.gatc.bedgraph" # exclude Dam sample from sample wildcard here!!!!
    output:
        bw="results/bigwig/{dir}/{bg_sample}.bw",
    params:
        extra=""
    #wildcard_constraints:
    #    sample="^((?!Dam).)*$"
    threads: config["resources"]["fastqc"]["cpu"]
    resources: 
        runtime=config["resources"]["fastqc"]["time"],
    conda: 
        "../envs/deeptools.yaml"
    log:
        "logs/bedgraph2bigwig/{dir}/bw_{bg_sample}.log"
    shell:
        "bedGraphToBigWig "
        "{params.extra} "
        "{input.bg} "
        "{input.cs} "
        "{output.bw} > {log} 2>&1"


rule average_wig:
    input:
        expand("results/bigwig/{dir}/{bg_sample}.bw", dir=DIRS, bg_sample=BG_SAMPLES),
    output:
        wig=temp("results/bigwig/average_bw/{bg_sample}.wig"),
    params:
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/wiggletools/wig_average_{bg_sample}.log"
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/average_wig.py"


rule wig2bigwig:
    input:
        wig="results/bigwig/average_bw/{bg_sample}.wig",
        cs=f"resources/{resources.genome}_chrom.sizes",
    output:
        "results/bigwig/average_bw/{bg_sample}.bw",
    params:
        extra="",
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/wigToBigWig/{bg_sample}.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "wigToBigWig {input.wig} {input.cs} {output}"

        