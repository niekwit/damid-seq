def targets():
    """Returns file targets for rule all
    """
    TARGETS = [
        "results/qc/multiqc/multiqc.html",
        #"results/plots/PCA.pdf",
        #"results/plots/scree.pdf",
        "results/plots/sample_correlation.pdf",
        "results/plots/heatmap.pdf",
        "results/plots/profile_plot.pdf",
        "results/plots/peaks/feature_distributions.pdf",
        "results/plots/peaks/distance_to_tss.pdf",
        "results/plots/mapping_rates.pdf",
        expand("results/peaks/overlapping_peaks/{bg_sample}.annotated.txt", bg_sample=BG_SAMPLES),
        ]

    return TARGETS


def dirs():
    """Each dir contains one replicate sets of fastq files
    """
    DIRS = glob.glob("reads/*")
    DIRS = [os.path.basename(d) for d in DIRS]

    if len(DIRS) == 0:
        raise ValueError("No replicate directories found in reads directory...")
        
    return DIRS
    

def samples(bedgraph=False):
    """Checks sample names/files and returns sample wildcard values for Snakemake
    """
    SAMPLES = csv["sample"]
    
    # Check if sample names contain any characters that are not alphanumeric or underscore
    illegal = []
    for sample in SAMPLES:
        if not re.match("^[a-zA-Z0-9_]*$", sample):
            illegal.append(sample)
    if len(illegal) != 0:
        illegal = "\n".join(illegal)
        raise ValueError(f"Following samples contain illegal characters:\n{illegal}")

    # Check if sample names match file names
    not_found = []
    for sample in SAMPLES:
        for dir in DIRS:
            if paired_end:
                r1= f"reads/{dir}/{sample}_R1_001.fastq.gz"
                r2= f"reads/{dir}/{sample}_R2_001.fastq.gz"
                if not os.path.isfile(r1):
                    not_found.append(r1)
                if not os.path.isfile(r2):
                    not_found.append(r2)
            else:
                r1= f"reads/{dir}/{sample}.fastq.gz"
                if not os.path.isfile(r1):
                    not_found.append(r1)
    if len(not_found) != 0:
        not_found = "\n".join(not_found)
        raise ValueError(f"Following files not found:\n{not_found}")
    
    # Only return non-Dam samples if bedgraph is True
    if bedgraph:
        SAMPLES = [s for s in SAMPLES if not "Dam" in s]

    return SAMPLES


def paired_end():
    """Checks if paired-end reads are used
    """
    # Get one fastq file
    reads = glob.glob("reads/*/*.gz")
    assert len(reads) != 0, "No fastq files found..."
    fastq = reads[0]

    # Check file extension to see if paired-end reads are used
    if fastq.endswith("_R1_001.fastq.gz"):
        return True
    elif fastq.endswith("_R2_001.fastq.gz"):
        return True
    else:
        return False
    

def computematrix_args(region_labels=False):
    """Returns computeMatrix arguments as string based on config file.
    """
    # Add mode argument
    mode = config["deeptools"]["matrix"]["mode"]
    if mode == "scale-regions".lower():
        rbl = config["deeptools"]["matrix"]["regionBodyLength"]
        args = f"scale-regions --regionBodyLength {rbl} "
    elif mode == "reference-point".lower():
        rp = config["deeptools"]["matrix"]["referencePoint"]
        args = f"reference-point --referencePoint {rp} "
    else:
        raise ValueError(f"Deeptools matrix mode {mode} not supported")
    
    
    # Add common arguments
    b = config["deeptools"]["matrix"]["upstream"]
    a = config["deeptools"]["matrix"]["downstream"]
    bs = config["deeptools"]["matrix"]["binSize"]
    atb = config["deeptools"]["matrix"]["averageTypeBins"]
        
    args = f"{args} --upstream {b} --downstream {a} --binSize {bs} --averageTypeBins {atb} "

    # Add region argument
    r = config["deeptools"]["matrix"]["regionsFileName"]
    no_whole_genome = config["deeptools"]["matrix"]["no_whole_genome"]

    # Check if files in r exist
    if r != "":
        _r = r.split(",")
        for region in _r:
            if not os.path.isfile(region):
                raise ValueError(f"File {region} not found...")

    if region_labels: # for plotProfile/Heatmap
        # Multiple files may be parsed in config file
            labels = r.split(",")
            labels = [os.path.basename(l).split(".")[0] for l in labels]

    if no_whole_genome and r:
        if region_labels: 
            return " ".join(labels)
        args = f"{args} --regionsFileName {r} "
    elif not no_whole_genome and r:
        if region_labels: # for plotProfile/Heatmap
            # Multiple files may be parsed in config file
            labels.insert(0, f'"Genome-wide ({resources.genome})"')
            return " ".join(labels)
        args = f"{args} --regionsFileName {resources.gtf} {r} "
    else: 
        args = f"{args} --regionsFileName {resources.gtf} "
        if region_labels: # for plotProfile/Heatmap
            return f'"Genome-wide ({resources.genome})"'
    return args

