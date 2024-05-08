from snakemake.logging import logger

if config["peak_calling_perl"]["run"]:
    fdr = config["peak_calling_perl"]["fdr"]
    rule peak_calling_perl:
        input:
            fp="resources/find_peaks",
            bg="results/bedgraph/{dir}/{bg_sample}-vs-Dam.kde-norm.gatc.bedgraph",
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
        fdr = peak_fdr("macs2_narrow")
        
        rule peak_calling_MACS2_narrow:
            input: 
                treatment=expand("results/bam/{dir}/{{bg_sample}}.bam", bg_sample=BG_SAMPLES),
                control=expand("results/bam/{dir}/{dam_sample}.bam", dam_sample=DAM_SAMPLES),
            output:
                multiext(f"results/macs2_narrow/fdr{fdr}/{{bg_sample}}",
                        "_peaks.xls",
                        "_peaks.narrowPeak",
                        "_summits.bed"
                        )
            params:
                macs2_params(),
            log:
                f"logs/macs2_narrow/fdr{fdr}/{{bg_sample}}.log"
            wrapper:
                f"{wrapper_version}/bio/macs2/callpeak"


        rule peak_annotation_plots_macs2_narrow:
            input:
                gtf=resources.gtf,
                bed=expand(f"results/macs2_narrow/fdr{fdr}/{{bg_sample}}_peaks.narrowPeak", bg_sample=BG_SAMPLES),
            output:
                fd=report(f"results/plots/macs2_narrow/fdr{fdr}/feature_distributions.pdf", caption="../report/feature_distributions.rst", category="Peak annotation MACS2 narrow"),
                dt=report(f"results/plots/macs2_narrow/fdr{fdr}/distance_to_tss.pdf", caption="../report/distance_to_tss.rst", category="Peak annotation MACS2 narrow"),
            params:
                extra="",
            threads: config["resources"]["plotting"]["cpu"]
            resources:
                runtime=config["resources"]["plotting"]["time"]
            log:
                f"logs/plotting/macs2_narrow_peak_annotation_plots_fdr{fdr}.log"
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/peak_annotation_plots.R"

        
        rule annotate_peaks_macs_narrow:
            input:
                bed=f"results/macs2_narrow/fdr{fdr}/{{bg_sample}}_peaks.narrowPeak",
                adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
                gtf=resources.gtf,
            output:
                txt=report(f"results/macs2_narrow/fdr{fdr}/{{bg_sample}}.annotated.txt", caption="../report/annotated_peaks.rst", category="Annotated peaks MACS2 narrow"),
            params:
                extra=""
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                f"logs/annotate_peaks_macs2_narrow/fdr{fdr}/{{bg_sample}}.log"
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/annotate_peaks.R"

        
        rule get_gene_names_macs2:
            input:
                txt=f"results/macs2_narrow/fdr{fdr}/{{bg_sample}}.annotated.txt"
            output:
                ids=f"results/macs2_narrow/fdr{fdr}/{{bg_sample}}.geneIDs.txt"
            threads: 1
            resources:
                runtime=5
            log:
                f"logs/geneIDs_peaks_macs2_narrow_fdr{fdr}/{{bg_sample}}.log"
            conda:
                "../envs/deeptools.yaml"
            shell:
                "sed '1d' {input.txt} | "
                "awk '{{print $(NF-4),$(NF-1)}}' |"
                " sort | "
                "uniq > {output.ids}"

        
        rule peak_calling_MACS2_narrow_single:
            input: 
                treatment=expand("results/bam/{dir}/{bg_sample}.bam", dir=DIRS, bg_sample=BG_SAMPLES),
                control=expand("results/bam/{dir}/{dam_sample}.bam", dam_sample=DAM_SAMPLES, dir=DIRS),
            output:
                multiext(f"results/macs2_narrow_single/fdr{fdr}/{{bg_sample}}",
                        "_peaks.xls",
                        "_peaks.narrowPeak",
                        "_summits.bed"
                        )
            params:
                macs2_params(),
            log:
                f"logs/macs2_narrow_single/fdr{fdr}/{{bg_sample}}.log"
            wrapper:
                f"{wrapper_version}/bio/macs2/callpeak"

        
        rule overlapping_peaks_macs2_narrow:
            input: 
                peaks=expand(f"results/macs2_narrow_single/fdr{fdr}/{{bg_sample}}_peaks.narrowPeak", fdr=peak_fdr("macs2_narrow"), bg_sample=BG_SAMPLES),
            output:
                f"results/macs2_narrow_single/fdr{fdr}/overlapping_peaks/{{bg_sample}}.overlap.bed"
            params:
                extra=""
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                f"logs/overlapping_peaks_macs2_narrow_fdr{fdr}/{{bg_sample}}.log"
            conda:
                "../envs/peak_calling.yaml"
            shell:
                "bedtools multiinter {params.extra} -i {input.beds} > {output}"
    
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
                paired_end=paired_end,
                mode="broad",
                fdr=fdr,
                genome=resources.genome,
                data_dir=lambda w, output: os.path.dirname(output[0]),
                extra=config["peak_calling_macs2"]["extra"]
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/macs2_broad_single/fdr{fdr}/{dir}/{bg_sample}.log"
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
                "logs/consensus_peaks_macs2_broad_fdr{fdr}/{bg_sample}.log"
            conda:
                "../envs/peak_calling.yaml"
            shell:
                "bedtools multiinter {params.extra} -i {input.peaks} > {output}"

        
        rule filter_consensus_peaks_macs2_broad:
            input:
                bed="results/macs2_broad/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.bed",
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
                "logs/filter_consensus_peaks/macs2_broad/fdr{fdr}/{{bg_sample}}.log"
            conda:
                "../envs/peak_calling.yaml"
            script:
                "../scripts/filter_consensus_peaks.py"

                    
        rule peak_annotation_plots_macs2_broad:
            input:
                gtf=resources.gtf,
                bed=expand("results/macs2_broad/fdr{fdr}/consensus_peaks/{bg_sample}.overlap.filtered.bed", fdr=fdr,bg_sample=BG_SAMPLES),
            output:
                fd=report("results/plots/macs2_broad/fdr{fdr}/feature_distributions_overlap.pdf", caption="../report/feature_distributions.rst", category="Peak annotation MACS2 broad"),
                dt=report("results/plots/macs2_broad/fdr{fdr}/distance_to_tss.pdf_overlap", caption="../report/distance_to_tss.rst", category="Peak annotation MACS2 broad"),
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