rule peak_calling:
    input:
        fp="resources/find_peaks",
        bg=expand("results/bedgraph/{dir}/{sample, ^((?!Dam).)*$}-vs-Dam.gatc.bedgraph", dir=DIRS, sample=SAMPLES), # exclude Dam sample from sample wildcard
    output:
        ""

rules peaks2genes:
    input:
        fp="resources/find_peaks",
        gtf=resources.gtf,
        peaks="",
    output:
        "",