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


rule index_fasta:
    input:
        resources.fasta,
    output:
        f"{resources.fasta}.fai",
    log:
        "logs/resources/index_fasta.log"
    threads: config["resources"]["plotting"]["cpu"]
    resources: 
        runtime=config["resources"]["plotting"]["time"]
    wrapper:
        "v3.3.3/bio/samtools/faidx"


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
        "{output} > {log} 2>&1"


use rule install_damidseq_pipeline_software as install_find_peak_software with:
    output:
        directory("resources/find_peaks"),
    params:
        url="https://github.com/owenjm/find_peaks.git",
    log:
        "logs/resources/install_find_peaks_software.log"


rule create_blacklist:
    input:
        gtf=resources.gtf,
    output:
        bed="resources/blacklist.bed",
        txt="resources/genes.txt",
    params:
        genes=config["fusion_genes"],
    log:
        "logs/resources/create_blacklist.log"
    threads: config["resources"]["deeptools"]["cpu"]
    resources: 
        runtime=config["resources"]["deeptools"]["time"]
    conda:
        "../envs/peak_calling.yaml"
    script:
        "../scripts/create_blacklist.py"
    