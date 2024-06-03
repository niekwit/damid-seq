rule get_fasta:
    output:
        resources.fasta,
    retries: 3
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log"
    threads: config["resources"]["plotting"]["cpu"]
    resources: 
        runtime=config["resources"]["plotting"]["time"]
    conda:
        "../envs/damid.yaml"
    script:
        "../scripts/get_resource.sh"


use rule get_fasta as get_gtf with:
        output:
            resources.gtf,
        params:
            url=resources.gtf_url,
        log:
            "logs/resources/get_gtf.log"


rule install_damidseq_pipeline_software:
    output:
        directory("resources/damidseq_pipeline"),
    params:
        url="https://github.com/owenjm/damidseq_pipeline.git",
        version="-b v1.5.3",
    retries: 3
    log:
        "logs/resources/install_find_peaks_software.log"
    threads: 1
    resources: 
        runtime=5
    conda:
        "../envs/damid.yaml"
    shell:
        "git clone "
        "{params.url} "
        "{params.version} "
        "{output} > {log} 2>&1"


use rule install_damidseq_pipeline_software as install_find_peak_software with:
    output:
        directory("resources/find_peaks"),
    params:
        url="https://github.com/owenjm/find_peaks.git",
        version=""
    log:
        "logs/resources/install_find_peaks_software.log"


rule create_annotation_file:
    input:
        gtf=resources.gtf,
    output:
        rdata=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
    log:
        "logs/resources/create_annotation_file.log"
    threads: config["resources"]["plotting"]["cpu"]
    resources: 
        runtime=config["resources"]["plotting"]["time"]
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/create_annotation_file.R"


rule masked_fasta:
    input:
        fa=resources.fasta,
        gtf=resources.gtf,
    output:
        out=f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.fa",
    params:
        g2m=maskedgenes,
        genome=resources.genome,
    log:
        "logs/resources/masked_fasta.log"
    threads: config["resources"]["plotting"]["cpu"]
    resources: 
        runtime=config["resources"]["plotting"]["time"]
    conda:
        "../envs/peak_calling.yaml"
    script:
        "../scripts/mask_fasta.py"


rule index_fasta:
    input:
        f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.fa",
    output:
        f"{resources.fasta}.fai",
    log:
        "logs/resources/index_fasta.log"
    threads: config["resources"]["plotting"]["cpu"]
    resources: 
        runtime=config["resources"]["plotting"]["time"]
    wrapper:
        f"{wrapper_version}/bio/samtools/faidx"


rule chrom_sizes:
    input:
        fa=resources.fasta,
        fai=f"{resources.fasta}.fai",
    output:
        f"resources/{resources.genome}_chrom.sizes",
    log:
        "logs/resources/chrom_sizes.log"
    threads: config["resources"]["plotting"]["cpu"]
    resources: 
        runtime=config["resources"]["plotting"]["time"]
    conda:
        "../envs/damid.yaml"
    shell:
        "awk '{{print $1,$2}}' {input.fai} | "
        r"sed 's/ /\t/'  > {output}"


rule make_gatc_tracks:
    input:
        sw="resources/damidseq_pipeline",
        fa=f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.fa",
    output:
        f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.GATC.gff",
    params:
        genome=lambda wildcards, output: output[0][:-9]
    threads: config["resources"]["fastqc"]["cpu"],
    resources:
        runtime=config["resources"]["fastqc"]["time"],
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
        ref=f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.fa",
    output:
        multiext(
            f"resources/bowtie2_index/{resources.genome}_{resources.build}_{maskedgenes}.masked/index",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        f"logs/bowtie2_build_index/{resources.genome}_{resources.build}.log",
    params:
        extra="",  # optional parameters
    threads: config["resources"]["index"]["cpu"]
    resources:
        runtime=config["resources"]["index"]["time"],
    wrapper:
        f"{wrapper_version}/bio/bowtie2/build"


if config["plasmid_fasta"] != "none":
    rule bowtie2_build_index_plasmid:
        input:
            ref=config["plasmid_fasta"],
        output:
            multiext(
                f"resources/bowtie2_index/{plasmid_name}/index",
                ".1.bt2",
                ".2.bt2",
                ".3.bt2",
                ".4.bt2",
                ".rev.1.bt2",
                ".rev.2.bt2",
            ),
        log:
            f"logs/bowtie2_build_index/{plasmid_name}.log",
        params:
            extra="",  # optional parameters
        threads: config["resources"]["index"]["cpu"]
        resources:
            runtime=config["resources"]["index"]["time"],
        wrapper:
            f"{wrapper_version}/bio/bowtie2/build"    
