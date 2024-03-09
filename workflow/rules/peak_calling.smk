if config["peak_calling_perl"]["run"]:
    rule peak_calling:
        input:
            fp="resources/find_peaks",
            bg="results/bedgraph/{dir}/{bg_sample}-vs-Dam.kde-norm.gatc.bedgraph",
        output:
            gff="results/peaks/{dir}/{bg_sample}.peaks.gff",
            data="results/peaks/{dir}/{bg_sample}.data",
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
            "logs/find_peaks/{dir}/{bg_sample}.log"
        script:
            "../scripts/run_find_peaks.py"


if config["peak_calling_macs2"]["run"]:
    if config["peak_calling_macs2"]["mode"] == "narrow":
        rule peak_calling_MACS2_narrow:
            input:
                treatment=expand("results/bam/{dir}/{bg_sample}{bamext}.bam", bg_sample=BG_SAMPLES, bamext=BAM_EXT, dir=DIRS),
                control=expand("results/bam/{dir}/{dam_sample}{bamext}.bam", dam_sample=DAM_SAMPLES, bamext=BAM_EXT, dir=DIRS),
            output:
                multiext("results/macs2_narrow/{bg_sample}",
                        "_peaks.xls",
                        "_peaks.narrowPeak",
                        "_summits.bed"
                        )
            params:
                macs2_params(),
            log:
                "logs/macs2_narrow/{bg_sample}.log"
            wrapper:
                "v3.4.1/bio/macs2/callpeak"
        

        rule peak_annotation_plots_macs2_narrow:
            input:
                gtf=resources.gtf,
                bed=expand("results/macs2_narrow/{bg_sample}_peaks.narrowPeak", bg_sample=BG_SAMPLES),
            output:
                fd=report("results/plots/macs2_narrow/feature_distributions.pdf", caption="../report/feature_distributions.rst", category="Peak annotation"),
                dt=report("results/plots/macs2_narrow/distance_to_tss.pdf", caption="../report/distance_to_tss.rst", category="Peak annotation MACS2 narrow"),

        
        rule annotate_peaks_macs_narrow:
            input:
                bed="results/macs2_narrow/{bg_sample}_peaks.narrowPeak",
                adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
                gtf=resources.gtf,
            output:
                txt=report("results/macs2_narrow/{bg_sample}.annotated.txt", caption="../report/annotated_peaks.rst", category="Annotated peaks MACS2 narrow"),
            params:
                extra=""
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/annotate_peaks_macs2_narrow/{bg_sample}.log"
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/annotate_peaks.R"

    elif config["peak_calling_macs2"]["mode"] == "broad":
        rule peak_calling_MACS2_broad:
            input:
                treatment=expand("results/bam/{dir}/{bg_sample}{bamext}.bam", bg_sample=BG_SAMPLES, bamext=BAM_EXT, dir=DIRS),
                control=expand("results/bam/{dir}/{dam_sample}{bamext}.bam", dam_sample=DAM_SAMPLES, bamext=BAM_EXT, dir=DIRS),
            output:
                multiext("results/macs2_broad/{bg_sample}",
                        "_peaks.xls",
                        "_peaks.broadPeak",
                        "_peaks.gappedPeak"
                        )
            params:
                macs2_params(),
            log:
                "logs/macs2_broad/{bg_sample}.log"
            wrapper:
                "v3.4.1/bio/macs2/callpeak"
        

        rule peak_annotation_plots_macs2_broad:
            input:
                gtf=resources.gtf,
                bed=expand("results/macs2_broad/{bg_sample}_peaks.broadPeak", bg_sample=BG_SAMPLES),
            output:
                fd=report("results/plots/macs2_broad/feature_distributions.pdf", caption="../report/feature_distributions.rst", category="Peak annotation"),
                dt=report("results/plots/macs2_broad/distance_to_tss.pdf", caption="../report/distance_to_tss.rst", category="Peak annotation MACS2 broad"),

        
        rule annotate_peaks_macs_broad:
            input:
                bed="results/macs2_broad/{bg_sample}_peaks.broadPeak",
                adb=f"resources/{resources.genome}_{resources.build}_annotation.Rdata",
                gtf=resources.gtf,
            output:
                txt=report("results/macs2_broad/{bg_sample}.annotated.txt", caption="../report/annotated_peaks.rst", category="Annotated peaks MACS2 broad"),
            params:
                extra=""
            threads: config["resources"]["deeptools"]["cpu"]
            resources:
                runtime=config["resources"]["deeptools"]["time"]
            log:
                "logs/annotate_peaks_macs2_broad/{bg_sample}.log"
            conda:
                "../envs/R.yaml"
            script:
                "../scripts/annotate_peaks.R"