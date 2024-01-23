rule damidseq_pipeline: # ignore dir wildcard in expand statement (double braces)
    input:
        git="resources/damidseq_pipeline",
        flag=expand("results/trimmed/{dir}/{sample}.flag", dir=DIRS, sample=SAMPLES),
        gatc=f"resources/{resources.genome}.GATC.gff",
        idx=multiext(
            f"resources/bowtie2_index/{resources.genome}/index",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    params:
        idxdir=f"resources/bowtie2_index/{resources.genome}/index",
        paired=config["paired_end"],
        extra=config["extra"],
    output:
        bg=directory(expand("results/bedgraph/{dir}", dir=DIRS)),
        bam=directory(expand("results/bam/{dir}", dir=DIRS)),
        bf=expand("results/bedgraph/{dir}/{bg_sample}-vs-Dam.kde-norm.gatc.bedgraph", dir=DIRS, bg_sample=BG_SAMPLES),
        #bg=expand("results/bedgraph/{dir}/{sample,^((?!Dam).)*$}-vs-Dam.gatc.bedgraph", dir=DIRS, sample=SAMPLES), # exclude Dam sample from sample wildcard 
    #wildcard_constraints:
    #    sample="^((?!Dam).)*$" # exclude Dam sample from sample wildcard 
    conda:
        "../envs/damid.yaml"
    threads: config["resources"]["damid"]["cpu"]
    resources:
        runtime=config["resources"]["damid"]["time"]
    log:
        expand("logs/damidseq_pipeline/{dir}/damidseq_pipeline.log", dir=DIRS)
    script:
        "../scripts/damidseq_pipeline.py"
