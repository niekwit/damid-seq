if config["peak_calling_perl"]["run"]:
    fdr = config["peak_calling_perl"]["fdr"]
    rule peak_calling_perl:
        input:
            fp="resources/find_peaks",
            bg="results/bedgraph/{dir}/{bg_sample}-vs-Dam-norm.gatc.bedgraph",
        output:
            gff="results/peaks/fdr{fdr}/{dir}/{bg_sample}.peaks.gff",
            data="results/peaks/fdr{fdr}/{dir}/{bg_sample}.data",
        params:
            outdir=lambda w, output: os.path.dirname(output["gff"]),
            n=config["peak_calling_perl"]["iterations"],
            fdr=config["peak_calling_perl"]["fdr"],
            frac=config["peak_calling_perl"]["fraction"],
            mc=config["peak_calling_perl"]["min_count"],
            mq=config["peak_calling_perl"]["min_quantile"],
            step=config["peak_calling_perl"]["step"],
            up=config["peak_calling_perl"]["unified_peaks"],
        threads: config["resources"]["deeptools"]["cpu"]
        resources:
            runtime=config["resources"]["deeptools"]["time"]
        conda:
            "../envs/peak_calling.yaml"
        log:
            "logs/find_peaks/fdr{fdr}/{dir}/{bg_sample}.log"
        script:
            "../scripts/run_find_peaks.py"

if config["peak_calling_macs2"]["run"]:
    if config["peak_calling_macs2"]["mode"] == "narrow":
        fdr = config["peak_calling_macs2"]["qvalue"]
        
        rule peak_calling_MACS2_narrow:
            input:
                bam="results/bam/{dir}/{bg_sample}.bam",
            output:
                multiext("results/macs2_broad/fdr{fdr}/{dir}/{bg_sample}",
                        "_peaks.xls",
                        "_peaks.narrowPeak",
                        "_summits.bed"
                        )
            params:
                outdir=lambda w, output: os.path.dirname(output[0]),
                paired_end=paired_end,
                mode="narrow",
                fdr=fdr,
                genome=resources.genome,
                data_dir=lambda w, input: os.path.dirname(input[0]),
                extra=config["peak_calling_macs2"]["extra"]
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/macs2_narrow/fdr{fdr}/{dir}/{bg_sample}.log"
            conda:
                "../envs/peak_calling.yaml"
            script:
                "../scripts/macs2.py"


        rule consensus_peaks_macs2_narrow:
            input: 
                peaks=expand("results/macs2_narrow/fdr{fdr}/{dir}/{{bg_sample}}_peaks.narrowPeak", dir=DIRS, fdr=fdr),
            output:
                "results/macs2_narrow/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.bed"
            params:
                extra=""
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/consensus_peaks_macs2_narrow/fdr{fdr}/{bg_sample}.log"
            conda:
                "../envs/peak_calling.yaml"
            shell:
                "bedtools multiinter {params.extra} -i {input.peaks} > {output}"


        rule filter_consensus_peaks_macs2_narrow:
            input:
                bed="results/macs2_narrow/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.bed",
                peaks=expand("results/macs2_narrow/fdr{fdr}/{dir}/{{bg_sample}}_peaks.narrowPeak", dir=DIRS, fdr=fdr),
                cs=f"resources/{resources.genome}_chrom.sizes",
            output:
                ext_bed="results/macs2_narrow/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.filtered.bed",
            params:
                k=config["consensus_peaks"]["keep"],
                max_size=config["consensus_peaks"]["max_size"],
                e=config["consensus_peaks"]["extend_by"],
            threads: 1
            resources:
                runtime=15
            log:
                "logs/filter_consensus_peaks/macs2_narrow/fdr{fdr}/{bg_sample}.log"
            conda:
                "../envs/peak_calling.yaml"
            script:
                "../scripts/filter_consensus_peaks.py"


        rule peak_annotation_plots_macs2_narrow:
            input:
                gtf=resources.gtf,
                bed=expand("results/macs2_narrow/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.filtered.bed", fdr=fdr,bg_sample=BG_SAMPLES),
            output:
                fd=report("results/plots/macs2_narrow/fdr{fdr}/feature_distributions.pdf", caption="../report/feature_distributions.rst", category="Peak annotation MACS2 narrow"),
                dt=report("results/plots/macs2_narrow/fdr{fdr}/distance_to_tss.pdf", caption="../report/distance_to_tss.rst", category="Peak annotation MACS2 narrow"),
            params:
                extra="",
            threads: config["resources"]["plotting"]["cpu"]
            resources:
                runtime=config["resources"]["plotting"]["time"]
            log:
                "logs/plotting/macs2_narrow_peak_annotation_plots/fdr{fdr}.log"
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/peak_annotation_plots.R"


        rule annotate_peaks_macs_narrow:
            input:
                bed="results/macs2_narrow/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.filtered.bed",
                adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
                gtf=resources.gtf,
            output:
                txt=report("results/macs2_narrow/fdr{fdr}/{bg_sample}.annotated.txt", caption="../report/annotated_peaks.rst", category="Annotated peaks MACS2 narrow"),
            params:
                extra=""
            log:
                "logs/annotate_peaks_macs2_narrow/fdr{fdr}/{bg_sample}.log"
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/annotate_peaks.R"


        rule get_gene_names_macs2:
            input:
                txt="results/macs2_narrow/fdr{fdr}/{bg_sample}.annotated.txt"
            output:
                ids="results/macs2_narrow/fdr{fdr}/{bg_sample}.geneIDs.txt"
            threads: 1
            resources:
                runtime=5
            log:
                "logs/geneIDs_peaks_macs2_narrow_fdr{fdr}/{bg_sample}.log"
            conda:
                "../envs/deeptools.yaml"
            shell:
                "sed '1d' {input.txt} | "
                "awk '{{print $(NF-4),$(NF-1)}}' | "
                "sort | "
                "uniq > {output.ids}"


        if config["consensus_peaks"]["enrichment_analysis"]["run"]:
            rule enrichment_analysis:
                input:
                    txt="results/macs2_narrow/fdr{fdr}/{bg_sample}.geneIDs.txt",
                output:
                    xlsx="results/macs2_narrow/fdr{fdr}/enrichment_analysis/{bg_sample}.xlsx",
                params:
                    extra="",
                    genome=resources.genome,
                    dbs=config["consensus_peaks"]["enrichment_analysis"]["dbs"]
                threads: config["resources"]["deeptools"]["cpu"]
                resources:
                    runtime=config["resources"]["deeptools"]["time"]
                log:
                    "logs/enrichment_analysis/macs2_narrow/fdr{fdr}/{bg_sample}.log"
                conda:
                    "../envs/R.yaml"
                script:
                    "../scripts/enrichment_analysis.R"

            rule plot_enrichment:
                input:
                    xlsx="results/macs2_narrow/fdr{fdr}/enrichment_analysis/{bg_sample}.xlsx",
                output:
                    plots=expand("results/plots/macs2_narrow/fdr{{fdr}}/enrichment_analysis/{{bg_sample}}/{db}.pdf", db=DBS),
                params:
                    terms=config["consensus_peaks"]["enrichment_analysis"]["terms"],
                    dirname=lambda w, output: os.path.dirname(output[0]),
                threads: config["resources"]["plotting"]["cpu"]
                resources:
                    runtime=config["resources"]["plotting"]["time"]
                log:
                    "logs/plot_enrichment/macs2_narrow/fdr{fdr}/{bg_sample}.log"
                conda:
                    "../envs/R.yaml"
                script:
                    "../scripts/plot_enrichment.R"


        rule count_reads_in_peaks:
        # Adapted from https://www.biostars.org/p/337872/#337890
            input:
                bam="results/bam/{dir}/{bg_sample}.bam",
                b="results/macs2_narrow/fdr{fdr}/{dir}/{bg_sample}_peaks.narrowPeak",
            output:
                total_read_count="results/macs2_narrow/fdr{fdr}/read_counts/{dir}/{bg_sample}.total.count",
                peak_read_count="results/macs2_narrow/fdr{fdr}/read_counts/{dir}/{bg_sample}.peak.count",
            params:
                extra="",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/bedtools_intersect/fdr{fdr}/{dir}/{bg_sample}.log"
            conda:
                "../envs/peak_calling.yaml"
            shell:
                "bedtools bamtobed "
                "{params.extra} "
                "-i {input.bam} | "
                "sort -k1,1 -k2,2n | "
                "tee >(wc -l > {output.total_read_count}) | "
                "bedtools intersect "
                "{params.extra} "
                "-sorted "
                "-c "
                "-a {input.b} "
                "-b stdin | "
                "awk '{{i+=$NF}}END{{print i}}' > "
                "{output.peak_read_count} "
                "{log}"
            

        rule plot_fraction_of_reads_in_peaks:
            input:
                total_read_count=expand("results/macs2_broad/fdr{fdr}/read_counts/{dir}/{bg_sample}.total.count", dir=DIRS, fdr=fdr, bg_sample=BG_SAMPLES),
                peak_read_count=expand("results/macs2_broad/fdr{fdr}/read_counts/{dir}/{bg_sample}.peak.count", dir=DIRS, fdr=fdr, bg_sample=BG_SAMPLES),
            output:
                plot="results/plots/macs2_broad/fdr{fdr}/frip.pdf",
                csv="results/macs2_broad/fdr{fdr}/frip.csv",
            params:
                extra="",
            threads: config["resources"]["plotting"]["cpu"]
            resources:
                runtime=config["resources"]["plotting"]["time"]
            log:
                "logs/plot_frip/fdr{fdr}.log"
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/plot_frip.R"
    
    elif config["peak_calling_macs2"]["mode"] == "broad":
        fdr = config["peak_calling_macs2"]["broad_cutoff"]

        rule peak_calling_MACS2_broad:
            input:
                bam="results/bam/{dir}/{bg_sample}.bam",
            output:
                multiext("results/macs2_broad/fdr{fdr}/{dir}/{bg_sample}",
                        "_peaks.xls",
                        "_peaks.broadPeak",
                        "_peaks.gappedPeak"
                        )
            params:
                outdir=lambda w, output: os.path.dirname(output[0]),
                paired_end=paired_end,
                mode="broad",
                fdr=fdr,
                genome=resources.genome,
                data_dir=lambda w, input: os.path.dirname(input[0]),
                extra=config["peak_calling_macs2"]["extra"]
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/macs2_broad/fdr{fdr}/{dir}/{bg_sample}.log"
            conda:
                "../envs/peak_calling.yaml"
            script:
                "../scripts/macs2.py"


        rule consensus_peaks_macs2_broad:
            input: 
                peaks=expand("results/macs2_broad/fdr{fdr}/{dir}/{{bg_sample}}_peaks.broadPeak", dir=DIRS, fdr=fdr),
            output:
                "results/macs2_broad/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.bed"
            params:
                extra=""
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/consensus_peaks_macs2_broad/fdr{fdr}/{bg_sample}.log"
            conda:
                "../envs/peak_calling.yaml"
            shell:
                "bedtools multiinter {params.extra} -i {input.peaks} > {output}"

        
        rule filter_consensus_peaks_macs2_broad:
            input:
                bed="results/macs2_broad/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.bed",
                peaks=expand("results/macs2_broad/fdr{fdr}/{dir}/{{bg_sample}}_peaks.broadPeak", dir=DIRS, fdr=fdr),
                cs=f"resources/{resources.genome}_chrom.sizes",
            output:
                ext_bed="results/macs2_broad/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.filtered.bed",
            params:
                k=config["consensus_peaks"]["keep"],
                max_size=config["consensus_peaks"]["max_size"],
                e=config["consensus_peaks"]["extend_by"],
            threads: 1
            resources:
                runtime=15
            log:
                "logs/filter_consensus_peaks/macs2_broad/fdr{fdr}/{bg_sample}.log"
            conda:
                "../envs/peak_calling.yaml"
            script:
                "../scripts/filter_consensus_peaks.py"

                    
        rule peak_annotation_plots_macs2_broad:
            input:
                gtf=resources.gtf,
                bed=expand("results/macs2_broad/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.filtered.bed", fdr=fdr,bg_sample=BG_SAMPLES),
            output:
                fd=report("results/plots/macs2_broad/fdr{fdr}/feature_distributions.pdf", caption="../report/feature_distributions.rst", category="Peak annotation MACS2 broad"),
                dt=report("results/plots/macs2_broad/fdr{fdr}/distance_to_tss.pdf", caption="../report/distance_to_tss.rst", category="Peak annotation MACS2 broad"),
            params:
                extra="",
            threads: config["resources"]["plotting"]["cpu"]
            resources:
                runtime=config["resources"]["plotting"]["time"]
            log:
                "logs/plotting/macs2_broad_peak_annotation_plots_fdr{fdr}.log"
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/peak_annotation_plots.R"

        
        rule annotate_peaks_macs_broad:
            input:
                bed="results/macs2_broad/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.filtered.bed",
                adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
                gtf=resources.gtf,
            output:
                txt=report("results/macs2_broad/fdr{fdr}/{bg_sample}.annotated.txt", caption="../report/annotated_peaks.rst", category="Annotated peaks MACS2 broad"),
            params:
                extra=""
            log:
                "logs/annotate_peaks_macs2_broad_fdr{fdr}/{bg_sample}.log"
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/annotate_peaks.R" 


        rule get_gene_names_macs2:
            input:
                txt="results/macs2_broad/fdr{fdr}/{bg_sample}.annotated.txt"
            output:
                ids="results/macs2_broad/fdr{fdr}/{bg_sample}.geneIDs.txt"
            threads: 1
            resources:
                runtime=5
            log:
                "logs/geneIDs_peaks_macs2_broad_fdr{fdr}/{bg_sample}.log"
            conda:
                "../envs/deeptools.yaml"
            shell:
                "sed '1d' {input.txt} | "
                "awk '{{print $(NF-4),$(NF-1)}}' | "
                "sort | "
                "uniq > {output.ids}"


        if config["consensus_peaks"]["enrichment_analysis"]["run"]:
            rule enrichment_analysis:
                input:
                    txt="results/macs2_broad/fdr{fdr}/{bg_sample}.geneIDs.txt",
                output:
                    xlsx="results/macs2_broad/fdr{fdr}/enrichment_analysis/{bg_sample}.xlsx",
                params:
                    extra="",
                    genome=resources.genome,
                    dbs=config["consensus_peaks"]["enrichment_analysis"]["dbs"]
                threads: config["resources"]["deeptools"]["cpu"]
                resources:
                    runtime=config["resources"]["deeptools"]["time"]
                log:
                    "logs/enrichment_analysis/macs2_broad/fdr{fdr}/{bg_sample}.log"
                conda:
                    "../envs/R.yaml"
                script:
                    "../scripts/enrichment_analysis.R"

            
            rule plot_enrichment:
                input:
                    xlsx="results/macs2_broad/fdr{fdr}/enrichment_analysis/{bg_sample}.xlsx",
                output:
                    plots=expand("results/plots/macs2_broad/fdr{{fdr}}/enrichment_analysis/{{bg_sample}}/{db}.pdf", db=DBS),
                params:
                    terms=config["consensus_peaks"]["enrichment_analysis"]["terms"],
                    dirname=lambda w, output: os.path.dirname(output[0]),
                threads: config["resources"]["plotting"]["cpu"]
                resources:
                    runtime=config["resources"]["plotting"]["time"]
                log:
                    "logs/plot_enrichment/macs2_broad/fdr{fdr}/{bg_sample}.log"
                conda:
                    "../envs/R.yaml"
                script:
                    "../scripts/plot_enrichment.R"


        rule count_reads_in_peaks:
        # Adapted from https://www.biostars.org/p/337872/#337890
            input:
                bam="results/bam/{dir}/{bg_sample}.bam",
                b="results/macs2_broad/fdr{fdr}/{dir}/{bg_sample}_peaks.broadPeak",
            output:
                total_read_count="results/macs2_broad/fdr{fdr}/read_counts/{dir}/{bg_sample}.total.count",
                peak_read_count="results/macs2_broad/fdr{fdr}/read_counts/{dir}/{bg_sample}.peak.count",
            params:
                extra="",
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/bedtools_intersect/fdr{fdr}/{dir}/{bg_sample}.log"
            conda:
                "../envs/peak_calling.yaml"
            shell:
                "bedtools bamtobed "
                "{params.extra} "
                "-i {input.bam} | "
                "sort -k1,1 -k2,2n | "
                "tee >(wc -l > {output.total_read_count}) | "
                "bedtools intersect "
                "{params.extra} "
                "-sorted "
                "-c "
                "-a {input.b} "
                "-b stdin | "
                "awk '{{i+=$NF}}END{{print i}}' > "
                "{output.peak_read_count} "
                #"{log}"
            

        rule plot_fraction_of_reads_in_peaks:
            input:
                total_read_count=expand("results/macs2_broad/fdr{fdr}/read_counts/{dir}/{bg_sample}.total.count", dir=DIRS, fdr=fdr, bg_sample=BG_SAMPLES),
                peak_read_count=expand("results/macs2_broad/fdr{fdr}/read_counts/{dir}/{bg_sample}.peak.count", dir=DIRS, fdr=fdr, bg_sample=BG_SAMPLES),
            output:
                plot="results/plots/macs2_broad/fdr{fdr}/frip.pdf",
                csv="results/macs2_broad/fdr{fdr}/frip.csv",
            params:
                extra="",
            threads: config["resources"]["plotting"]["cpu"]
            resources:
                runtime=config["resources"]["plotting"]["time"]
            log:
                "logs/plot_frip/fdr{fdr}.log"
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/plot_frip.R"