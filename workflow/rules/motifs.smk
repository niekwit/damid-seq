if config["motif_analysis"]["run_analysis"]:
    rule convert_bed2fasta:
        input:
            bed="results/peaks/overlapping_peaks/{bg_sample}.filtered.bed",
            fasta=resources.fasta
        output:
            out="results/motifs/fasta/{bg_sample}.fa"
        params:
            extra="",
        threads: config["resources"]["deeptools"]["cpu"]
        resources: 
            runtime=config["resources"]["deeptools"]["time"]
        conda:
            "../envs/peak_calling.yaml"
        log:
            "logs/motifs/convert_bed2fasta/{bg_sample}.log"
        script:
            "../scripts/convert_bed2fasta.py"


    rule create_background_fasta:
        input:
            bed="results/peaks/overlapping_peaks/{bg_sample}.filtered.bed",
            fasta=resources.fasta,
            cs=f"resources/{resources.genome}_chrom.sizes",
        output:
            out="results/motifs/fasta/background_{bg_sample}.fa"
        params:
            extra="",
        threads: config["resources"]["deeptools"]["cpu"]
        resources: 
            runtime=config["resources"]["deeptools"]["time"]
        conda:
            "../envs/peak_calling.yaml"
        log:
            "logs/motifs/create_background_fasta/{bg_sample}.log"
        script:
            "../scripts/create_background_fasta.py"


    rule detect_motifs:
        input:
            fa="results/motifs/fasta/{bg_sample}.fa",
            bgfa="results/motifs/fasta/background_{bg_sample}.fa",
        output:
            outdir=directory("results/motifs/{bg_sample}/"),
        params:
            _len=config["motif_analysis"]["len"],
            mask=config["motif_analysis"]["mask"],
            extra="",
        threads: config["resources"]["deeptools"]["cpu"]
        resources: 
            runtime=config["resources"]["deeptools"]["time"]
        conda:
            "../envs/peak_calling.yaml"
        log:
            "logs/homer/{bg_sample}.log"
        shell:
            "findMotifsGenome.pl "
            "{input.fa} "
            "fasta "
            "{output.outdir} "
            "-fasta {input.bgfa} "
            "-len {params._len} "
            "{params.mask} "
            "-p {threads} "
            "{params.extra} "
            "> {log} 2>&1"