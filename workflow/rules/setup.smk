rule make_gatc_tracks:
    input:
        fa=resources.fasta,
    output:
        gatc=f"resources/{resources.genome}.GATC.gff",
    params:
        genome=f"resources/{resources.genome}",
    threads: config["resources"]["fastqc"]["cpu"],
    resources:
        time=config["resources"]["fastqc"]["time"],
    conda:
        "../envs/damid.yaml",
    log:
        "logs/make_gatc_tracks/{params.genome}.log",
    shell:
        "gatc.track.maker.pl "
        "--name={params.genome} "
        "{input.fa} > {log} 2>&1 "


rule bowtie2_build_index:
    input:
        ref=resources.fasta,
    output:
        multiext(
            f"resources/bowtie2_index/{resources.genome}/index",
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

