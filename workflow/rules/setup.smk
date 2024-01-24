rule make_gatc_tracks:
    input:
        sw="resources/damidseq_pipeline",
        fa=resources.fasta,
    output:
        "resources/genome.GATC.gff",
    params:
        genome=lambda wildcards, output: output[0][:-9]
    threads: config["resources"]["fastqc"]["cpu"],
    resources:
        time=config["resources"]["fastqc"]["time"],
    conda:
        "../envs/damid.yaml",
    log:
        "logs/make_gatc_tracks/tracks.log",
    shell:
        "perl resources/damidseq_pipeline/gatc.track.maker.pl "
        "--name={params.genome} "
        "{input.fa} > {log} 2>&1 "


rule bowtie2_build_index:
    input:
        ref=resources.fasta,
    output:
        multiext(
            "resources/bowtie2_index/index",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build_index/build.log",
    params:
        extra="",  # optional parameters
    threads: config["resources"]["index"]["cpu"]
    resources:
        runtime=config["resources"]["index"]["time"],
    wrapper:
        "v3.3.3/bio/bowtie2/build"

