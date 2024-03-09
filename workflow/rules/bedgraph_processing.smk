rule quantile_normalisation:
    input:
        bg=expand("results/bedgraph/{dir}/{bg_sample}-vs-Dam.kde-norm.gatc.bedgraph", dir=DIRS, bg_sample=BG_SAMPLES)
    output:
        bg=expand("results/bedgraph/{dir}/{bg_sample}-vs-Dam.quantile-norm.gatc.bedgraph", dir=DIRS, bg_sample=BG_SAMPLES)
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    #log:
    #    "logs/quantile_norm/{dir}/{bg_sample}.log"
    conda:
        "../envs/damid.yaml"
    shell:
        "python workflow/scripts/quantile_norm_bedgraph.py {input.bg}"

rule bedgraph2bigwig:
    input:
        cs=f"resources/{resources.genome}_chrom.sizes",
        bg="results/bedgraph/{dir}/{bg_sample}-vs-Dam.quantile-norm.gatc.bedgraph" # exclude Dam sample from sample wildcard here!!!!
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


rule reverse_log2:
    input:
        bg="results/bedgraph/{dir}/{bg_sample}-vs-Dam.kde-norm.gatc.bedgraph"
    output:
        bg="results/bedgraph/{dir}/{bg_sample}-vs-Dam.rev_log2.bedgraph",
    threads: 1
    resources:
        runtime=30
    log:
        "logs/rev_log2/{dir}/{bg_sample}.log"
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/reverse_log2.py"


use rule bedgraph2bigwig as bedgraph2bigwig_rev_log2 with:
    input:
        bg="results/bedgraph/{dir}/{bg_sample}-vs-Dam.rev_log2.bedgraph",
        cs=f"resources/{resources.genome}_chrom.sizes",
    output:
        bw="results/bigwig_rev_log2/{dir}/{bg_sample}.bw",
    log:
        "logs/bedgraph2bigwig_rev_log2/{dir}/bw_{bg_sample}.log"


use rule average_wig as average_wig_rev_log2 with:
    input:
        expand("results/bigwig_rev_log2/{dir}/{bg_sample}.bw", dir=DIRS, bg_sample=BG_SAMPLES),
    output:
        wig=temp("results/bigwig_rev_log2/average_bw/{bg_sample}.wig"),
    log:
        "logs/wiggletools_rev_log2/wig_average_{bg_sample}.log"


use rule wig2bigwig as wig2bigwig_rev_log2 with:
    input:
        wig="results/bigwig_rev_log2/average_bw/{bg_sample}.wig",
        cs=f"resources/{resources.genome}_chrom.sizes",
    output:
        "results/bigwig_rev_log2/average_bw/{bg_sample}.bw",
    log:
        "logs/wigToBigWig_rev_log2/{bg_sample}.log"
