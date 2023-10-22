rule bowtie2_build:
    input:
        ref=resources.fasta,
    output:
        multiext(
            f"resources/bowtie2_index/{resources.genome}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        "logs/bowtie2_build/build.log",
    params:
        extra="",  # optional parameters
    threads: config["damid"]["threads"]
    wrapper:
        "v2.6.0/bio/bowtie2/build"


rule create_gatc_fragments:
    input:
        resources.fasta
    output:
        f"resources/gatc_fragment_file_{resources.genome}.gff.gz"
    conda:
        "envs/damid.yml"
    threads: config["damid"]["threads"]
    log:
        "logs/create_gatc_fragment_file/gatc.log"
    shell:
        "perl gatc.track.maker.pl "
        f"--name={resources.genome} "
        f"{resources.fasta} 2> {log}"


rule damidseq_pipeline:
    input:
        gatc=f"resources/gatc_fragment_file_{resources.genome}.gff.gz",
        b2dir=f"resources/bowtie2_index/{resources.genome}",
    output:
        directory("results/bedgraph"),
    conda:
        "envs/damid.yml"
    threads: config["damid"]["threads"]
    log:
        "logs/damidseq_pipeline/damidseq_pipeline.log"
    shell:
        "cd reads/ && "
        "damidseq_pipeline "
        "--paired "
        "--gatc_frag_file={input.gatc "
        "--bowtie2_genome_dir={input.b2dir} 2> {log} && "
        "cd .."

