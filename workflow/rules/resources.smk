rule get_fasta:
    output:
        resources.fasta,
    log:
        "logs/resources/get_fasta.log",
    retries: 3
    conda:
        "../envs/damid.yaml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"],
    params:
        url=resources.fasta_url,
    script:
        "../scripts/get_resource.sh"


use rule get_fasta as get_gtf with:
    output:
        resources.gtf,
    log:
        "logs/resources/get_gtf.log",
    params:
        url=resources.gtf_url,


rule install_damidseq_pipeline_software:
    output:
        directory("resources/damidseq_pipeline"),
    log:
        "logs/resources/install_find_peaks_software.log",
    retries: 3
    conda:
        "../envs/damid.yaml"
    threads: 1
    resources:
        runtime=5,
    params:
        url="https://github.com/owenjm/damidseq_pipeline.git",
        version="-b v1.5.3",
    shell:
        "git clone " "{params.url} " "{params.version} " "{output} > {log} 2>&1"


rule install_find_peak_software:
    output:
        directory("resources/find_peaks"),
    log:
        "logs/resources/install_find_peaks_software.log",
    retries: 3
    conda:
        "../envs/damid.yaml"
    threads: 1
    resources:
        runtime=5,
    params:
        url="https://github.com/owenjm/find_peaks.git",
        commit="2259915d08c7d5ab0151d7b15e2292603d7a05c7",
    shell:
        "git clone {params.url} {output} && "
        "git -C {output} checkout {params.commit} > {log} 2>&1"


rule create_annotation_file:
    input:
        gtf=resources.gtf,
    output:
        rdata=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
    log:
        "logs/resources/create_annotation_file.log",
    conda:
        "../envs/R.yaml"
    threads: config["resources"]["trim"]["cpu"]
    resources:
        runtime=config["resources"]["trim"]["time"],
    script:
        "../scripts/create_annotation_file.R"


rule mask_fasta:
    input:
        fa=resources.fasta,
        gtf=resources.gtf,
    output:
        out=f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.fa",
    log:
        "logs/resources/masked_fasta.log",
    conda:
        "../envs/peak_calling.yaml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"],
    params:
        g2m=maskedgenes,
        genome=resources.genome,
        f2m=config["fusion_genes"]["feature_to_mask"],
    script:
        "../scripts/mask_fasta.py"


rule index_fasta:
    input:
        f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.fa",
    output:
        f"{resources.fasta}.fai",
    log:
        "logs/resources/index_fasta.log",
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"],
    wrapper:
        "v5.8.3/bio/samtools/faidx"


rule chrom_sizes:
    input:
        fa=resources.fasta,
        fai=f"{resources.fasta}.fai",
    output:
        f"resources/{resources.genome}_chrom.sizes",
    log:
        "logs/resources/chrom_sizes.log",
    conda:
        "../envs/damid.yaml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"],
    shell:
        "awk '{{print $1,$2}}' {input.fai} | " r"sed 's/\s/\t/'  > {output}"


rule chrom_order:
    input:
        f"resources/{resources.genome}_chrom.sizes",
    output:
        f"resources/{resources.genome}_chrom_order.txt",
    log:
        "logs/resources/chrom_order.log",
    conda:
        "../envs/damid.yaml"
    threads: 1
    resources:
        runtime=5,
    shell:
        r"sort -k 1,1V {input} > {output}"


rule make_gatc_tracks:
    input:
        fa=f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.fa",
    output:
        out=temp(
            f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.GATC.unsorted.gff"
        ),
    log:
        "logs/make_gatc_tracks/tracks.log",
    conda:
        "../envs/peak_calling.yaml"
    threads: config["resources"]["damid"]["cpu"]
    resources:
        runtime=45,
    params:
        genome=resources.genome,
        motif="GATC",
    script:
        "../scripts/gatc.track.maker.py"


rule sort_gatc_tracks:
    input:
        f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.GATC.unsorted.gff",
    output:
        f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.GATC.gff",
    log:
        "logs/make_gatc_tracks/sort_tracks.log",
    conda:
        "../envs/damid.yaml"
    threads: 1
    resources:
        runtime=10,
    shell:
        "sort -k1,1V -k4,4n {input} > {output} 2> {log}"


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
    threads: config["resources"]["index"]["cpu"]
    resources:
        runtime=config["resources"]["index"]["time"],
    params:
        extra="",  # optional parameters
    wrapper:
        "v5.8.3/bio/bowtie2/build"


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
        threads: config["resources"]["index"]["cpu"]
        resources:
            runtime=config["resources"]["index"]["time"],
        params:
            extra="",  # optional parameters
        wrapper:
            "v5.8.3/bio/bowtie2/build"
