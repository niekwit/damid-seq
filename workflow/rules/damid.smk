if config["plasmid_fasta"] == "none":
    if paired_end:
        rule bowtie2_align_pe:
            input:
                sample=["results/trimmed/{dir}/{sample}_1.fastq.gz",
                        "results/trimmed/{dir}/{sample}_2.fastq.gz"],
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
                bam=temp("results/bam/{dir}/{sample}.bam"),
            params:
                extra=config["bowtie2"]["extra"],
            threads: config["resources"]["bowtie2"]["cpu"]
            resources:
                runtime=config["resources"]["bowtie2"]["time"]
            log:
                "logs/bowtie2_align/{dir}/{sample}.log"
            wrapper:
                f"{wrapper_version}/bio/bowtie2/align"


        rule sort_bowtie2_bam:
            input:
                "results/bam/{dir}/{sample}.bam",
            output:
                "results/bam/{dir}/{sample}.sorted.bam",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/samtools/sort/{dir}/{sample}.log"
            wrapper:
                f"{wrapper_version}/bio/samtools/sort"

        
        rule index_bowtie2_bam:
            input:
                "results/bam/{dir}/{sample}.sorted.bam",
            output:
                "results/bam/{dir}/{sample}.sorted.bam.bai",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/samtools/index/{dir}/{sample}.log"
            wrapper:
                f"{wrapper_version}/bio/samtools/index"

        
        rule damidseq_pipeline: # Ignore dir wildcard in expand statement (double braces)
            input:
                git="resources/damidseq_pipeline",
                bam=expand("results/bam/{{dir}}/{sample}.sorted.bam", sample=SAMPLES),
                bai=expand("results/bam/{{dir}}/{sample}.sorted.bam.bai", sample=SAMPLES),
                gatc=f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.GATC.gff",
            output:
                bg=expand("results/bedgraph/{{dir}}/{bg_sample}-vs-Dam-norm.gatc.bedgraph", bg_sample=BG_SAMPLES),
            params:
                idx=f"resources/bowtie2_index/{resources.genome}_{resources.build}_{maskedgenes}.masked/index",
                binsize=config["damidseq_pipeline"]["binsize"],
                normalization_method=config["damidseq_pipeline"]["normalization"],
                extra=config["damidseq_pipeline"]["extra"],
            conda:
                "../envs/damid.yaml"
            threads: config["resources"]["damid"]["cpu"]
            resources:
                runtime=config["resources"]["damid"]["time"],
                tmpdir=config["resources"]["damid"]["tmpdir"],
            log:
                "logs/damidseq_pipeline/{dir}/damidseq_pipeline.log"
            script:
                "../scripts/damidseq_pipeline.py"
    else:
        rule bowtie2_align_se:
            input:
                sample=["results/trimmed/{dir}/{sample}.fastq.gz"],
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
                bam=temp("results/bam/{dir}/{sample}.bam"),
            params:
                extra=config["bowtie2"]["extra"]
            threads: config["resources"]["bowtie2"]["cpu"]
            resources:
                runtime=config["resources"]["bowtie2"]["time"]
            log:
                "logs/bowtie2_align/{dir}/{sample}.log"
            wrapper:
                f"{wrapper_version}/bio/bowtie2/align"


        rule sort_bowtie2_bam:
            input:
                "results/bam/{dir}/{sample}.bam",
            output:
                temp("results/bam/{dir}/{sample}.sorted.bam"),
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/samtools/index/{dir}/{sample}.log"
            wrapper:
                f"{wrapper_version}/bio/samtools/sort"


        rule index_bowtie2_bam:
            input:
                "results/bam/{dir}/{sample}.sorted.bam",
            output:
                temp("results/bam/{dir}/{sample}.sorted.bam.bai"),
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/samtools/index/{dir}/{sample}.log"
            wrapper:
                f"{wrapper_version}/bio/samtools/index"


        # Single-end fragments are extended to n bp
        # or to the closest GATC site
        rule extend_fragments:
            input:
                bam="results/bam/{dir}/{sample}.sorted.bam",
                bai="results/bam/{dir}/{sample}.sorted.bam.bai",
                gatc_gff=f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.GATC.gff",
                fai=f"{resources.fasta}.fai",
            output:
                bam=temp("results/bam/{dir}/{sample}.extended.bam"), 
            params:
                n=300,
                q=20,
            conda:
                "../envs/deeptools.yaml"
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/extend_reads/{dir}/{sample}.log"
            script:
                "../scripts/extend_fragments.py"

        '''
        rule sort_extended_bam:
            input:
                "results/bam/{dir}/{sample}.extended.bam",
            output:
                "results/bam/{dir}/{sample}.sorted.bam",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/samtools/sort/{dir}/{sample}.log"
            wrapper:
                f"{wrapper_version}/bio/samtools/sort"
        '''
        
        rule index_bam:
            input:
                "results/bam/{dir}/{sample}.extended.bam",
            output:
                "results/bam/{dir}/{sample}.extended.bam.bai",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/samtools/index/{dir}/{sample}.log"
            wrapper:
                f"{wrapper_version}/bio/samtools/index"
        
    
        # dir wildcard is escaped so that multiple dirs can be run in parallel
        rule damidseq_pipeline: # Ignore dir wildcard in expand statement (double braces)
            input:
                git="resources/damidseq_pipeline",
                bam=expand("results/bam/{{dir}}/{sample}.extended.bam", sample=SAMPLES),
                bai=expand("results/bam/{{dir}}/{sample}.extended.bam.bai", sample=SAMPLES),
                gatc=f"resources/{resources.genome}_{resources.build}_{maskedgenes}.masked.GATC.gff",
            output:
                bg=expand("results/bedgraph/{{dir}}/{bg_sample}-vs-Dam-norm.gatc.bedgraph", bg_sample=BG_SAMPLES),
            params:
                idx=f"resources/bowtie2_index/{resources.genome}_{resources.build}_{maskedgenes}.masked/index",
                binsize=config["damidseq_pipeline"]["binsize"],
                normalization_method=config["damidseq_pipeline"]["normalization"],
                extra=config["damidseq_pipeline"]["extra"],
            conda:
                "../envs/damid.yaml"
            threads: config["resources"]["damid"]["cpu"]
            resources:
                runtime=config["resources"]["damid"]["time"],
                tmpdir=config["resources"]["damid"]["tmpdir"],
            log:
                "logs/damidseq_pipeline/{dir}/damidseq_pipeline.log"
            script:
                "../scripts/damidseq_pipeline.py"
    '''
    old
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
            bf=expand("results/bedgraph/{{dir}}/{bg_sample}-vs-Dam-norm.gatc.bedgraph", bg_sample=BG_SAMPLES),
            bam=temp(expand("results/bam/{{dir}}/{sample}.bam", sample=SAMPLES)),
        params:
            idxdir=lambda wildcards, input: input["idx"][0][:-6],
            paired=paired_end,
            binsize=config["damidseq_pipeline"]["binsize"],
            trim_dir="trimmed",
            normalization_method=config["damidseq_pipeline"]["normalization"],
            extra=config["damidseq_pipeline"]["extra"],
        conda:
            "../envs/damid.yaml"
        threads: config["resources"]["damid"]["cpu"]
        resources:
            runtime=config["resources"]["damid"]["time"],
            tmpdir=config["resources"]["damid"]["tmpdir"],
        log:
            "logs/damidseq_pipeline/{dir}/damidseq_pipeline.log"
        script:
            "../scripts/damidseq_pipeline.py"
    '''
else:
    if paired_end:
        rule bowtie2_remove_plasmid_reads:
            input:
                flag="results/trimmed/{dir}/{sample}.flag",
                idx=multiext(
                    f"resources/bowtie2_index/{plasmid_name}/index",
                    ".1.bt2",
                    ".2.bt2",
                    ".3.bt2",
                    ".4.bt2",
                    ".rev.1.bt2",
                    ".rev.2.bt2",
                ),
            output:
                r1_fastq="results/trimmed_no_plasmid/{dir}/{sample}_1.fastq.gz",
                r2_fastq="results/trimmed_no_plasmid/{dir}/{sample}_2.fastq.gz",
                bam="results/bam/plasmid_reads/{dir}/{sample}.bam",
            params:
                idxdir=lambda wildcards, input: input["idx"][0][:-6],
                paired=paired_end,
                trim_dir="trimmed",
                out_base= lambda wildcards, output: output["r1_fastq"].replace("_1.fastq.gz", "_%.fastq.gz"),
                extra=config["bowtie2"]["extra"],
            conda:
                "../envs/damid.yaml"
            threads: config["resources"]["damid"]["cpu"]
            resources:
                runtime=config["resources"]["damid"]["time"]
            log:
                "logs/bowtie2_align_to_plasmid/{dir}/{sample}.log"
            script:
                "../scripts/bowtie2_align_to_plasmid.py"

        '''
        rule bowtie2_align:
            input:
                sample=["results/trimmed_no_plasmid/{dir}/{sample}_1.fastq.gz", "results/trimmed_no_plasmid/{dir}/{sample}_2.fastq.gz"],
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
                bam=temp(expand("results/bam/{{dir}}/{sample}.bam", sample=SAMPLES)),
            params:
                extra=config["bowtie2"]["extra"],
            threads: config["resources"]["bowtie2"]["cpu"]
            resources:
                runtime=config["resources"]["bowtie2"]["time"]
            log:
                "logs/bowtie2_align/{dir}/{sample}.log"
            wrapper:
                f"{wrapper_version}/bio/bowtie2/align"
        '''
        rule damidseq_pipeline: # ignore dir wildcard in expand statement (double braces)
            input:
                git="resources/damidseq_pipeline",
                r1_fastq=expand("results/trimmed_no_plasmid/{{dir}}/{sample}_1.fastq.gz", sample=SAMPLES),
                r2_fastq=expand("results/trimmed_no_plasmid/{{dir}}/{sample}_2.fastq.gz", sample=SAMPLES),
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
                bf=expand("results/bedgraph/{{dir}}/{bg_sample}-vs-Dam-norm.gatc.bedgraph", bg_sample=BG_SAMPLES),
                bam=temp(expand("results/bam/{{dir}}/{sample}.bam", sample=SAMPLES)),
            params:
                idxdir=lambda wildcards, input: input["idx"][0][:-6],
                paired=paired_end,
                binsize=config["damidseq_pipeline"]["binsize"],
                normalization_method=config["damidseq_pipeline"]["normalization"],
                extra=config["damidseq_pipeline"]["extra"],
            conda:
                "../envs/damid.yaml"
            threads: config["resources"]["damid"]["cpu"]
            resources:
                runtime=config["resources"]["damid"]["time"],
                tmpdir=config["resources"]["damid"]["tmpdir"],
            log:
                "logs/damidseq_pipeline/{dir}/damidseq_pipeline.log"
            script:
                "../scripts/damidseq_pipeline.py"
    else:
        rule bowtie2_remove_plasmid_reads:
            input:
                flag="results/trimmed/{dir}/{sample}.flag",
                idx=multiext(
                    f"resources/bowtie2_index/{plasmid_name}/index",
                    ".1.bt2",
                    ".2.bt2",
                    ".3.bt2",
                    ".4.bt2",
                    ".rev.1.bt2",
                    ".rev.2.bt2",
                ),
            output:
                bam="results/bam/plasmid_reads/{dir}/{sample}.bam",
                fastq="results/trimmed_no_plasmid/{dir}/{sample}.fastq.gz",
            params:
                idxdir=lambda wildcards, input: input["idx"][0][:-6],
                paired=paired_end,
                out_base= lambda wildcards, output: output["fastq"].replace(".fastq.gz", ""),
                extra=config["bowtie2"]["extra"],
            conda:
                "../envs/damid.yaml"
            threads: config["resources"]["damid"]["cpu"]
            resources:
                runtime=config["resources"]["damid"]["time"]
            log:
                "logs/bowtie2_align_to_plasmid/{dir}/{sample}.log"
            script:
                "../scripts/bowtie2_align_to_plasmid.py"


        rule damidseq_pipeline: # Ignore dir wildcard in expand statement (double braces)
            input:
                git="resources/damidseq_pipeline",
                fastq=expand("results/trimmed_no_plasmid/{{dir}}/{sample}.fastq.gz", sample=SAMPLES),
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
                bf=expand("results/bedgraph/{{dir}}/{bg_sample}-vs-Dam.norm.gatc.bedgraph", bg_sample=BG_SAMPLES),
                bam=temp(expand("results/bam/{{dir}}/{sample}.bam", sample=SAMPLES)),
            params:
                idxdir=lambda wildcards, input: input["idx"][0][:-6],
                paired=paired_end,
                trim_dir="trimmed_no_plasmid",
                binsize=config["damidseq_pipeline"]["binsize"],
                normalization_method=config["damidseq_pipeline"]["normalization"],
                extra=config["damidseq_pipeline"]["extra"],
            conda:
                "../envs/damid.yaml"
            threads: config["resources"]["damid"]["cpu"]
            resources:
                runtime=config["resources"]["damid"]["time"],
                tmpdir=config["resources"]["damid"]["tmpdir"],
            log:
                "logs/damidseq_pipeline/{dir}/damidseq_pipeline.log"
            script:
                "../scripts/damidseq_pipeline.py"


rule bam2bigwig:
    input:
        bam="results/bam/{dir}/{sample}.sorted.bam",
        bai="results/bam/{dir}/{sample}.sorted.bam.bai",
    output:
        "results/bigwig/bam2bigwig/{dir}/{sample}.bw"
    params:
        bs=config["deeptools"]["bamCoverage"]["binSize"],
        n=config["deeptools"]["bamCoverage"]["normalizeUsing"],
        extra=config["deeptools"]["bamCoverage"]["extra"]
    conda:
        "../envs/deeptools.yaml"
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"]
    log:
        "logs/bam2bigwig/{dir}/{sample}.log"
    shell:
        "bamCoverage "
        "-b {input.bam} "
        "-o {output} "
        "--binSize {params.bs} "
        "-p {threads} "
        "--normalizeUsing {params.n} "
        "{params.extra} "
        "> {log} 2>&1 "