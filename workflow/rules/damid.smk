# dir wildcard is escaped so that multiple dirs can be run in parallel
rule damidseq_pipeline: # ignore dir wildcard in expand statement (double braces)
    input:
        git="resources/damidseq_pipeline",
        flag=expand("results/trimmed/{{dir}}/{sample}.flag", sample=SAMPLES),
        gatc=f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.GATC.gff",
        idx=multiext(
            f"resources/bowtie2_index/{resources.genome}_{resources.build}_{maskedgenes}.masked/index",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        bf=expand("results/bedgraph/{{dir}}/{bg_sample}-vs-Dam.kde-norm.gatc.bedgraph", bg_sample=BG_SAMPLES),
        bam=expand("results/bam/{{dir}}/{sample}-ext300.bam", sample=SAMPLES),
    params:
        idxdir=lambda wildcards, input: input["idx"][0][:-6],
        paired=paired_end,
        extra=config["extra"],
    conda:
        "../envs/damid.yaml"
    threads: config["resources"]["damid"]["cpu"]
    resources:
        runtime=config["resources"]["damid"]["time"]
    log:
        "logs/damidseq_pipeline/{dir}/damidseq_pipeline.log"
    script:
        "../scripts/damidseq_pipeline.py"
