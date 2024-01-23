def dirs():
    """Each dir contains one replicate sets of fastq files
    """
    DIRS = glob.glob("reads/*")
    DIRS = [os.path.basename(d) for d in DIRS]
        
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
        raise ValueError(f"ERROR: following samples contain illegal characters:\n{illegal}")

    # Check if sample names match file names
    not_found = []
    for sample in SAMPLES:
        for dir in DIRS:
            if config["paired_end"]:
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
        raise ValueError(f"ERROR: following files not found:\n{not_found}")
    
    # Only return non-Dam samples if bedgraph is True
    if bedgraph:
        SAMPLES = [s for s in SAMPLES if not "Dam" in s]

    return SAMPLES


def targets():
    """Returns file targets for rule all
    """
    TARGETS = [
        #expand("results/bedgraph/{dir}/{sample, ^((?!Dam).)*$}-vs-Dam.gatc.bedgraph", dir=DIRS, sample=SAMPLES),
        "results/qc/multiqc/multiqc.html",
        #expand("results/bedgraph/{dir}", dir=DIRS),
        #expand("results/bigwig/{dir}/{bg_sample}.bw", dir=DIRS, bg_sample=BG_SAMPLES),
        expand("results/bigwig/average_bw/{bg_sample}.bw", bg_sample=BG_SAMPLES),
        "results/deeptools/PCA.pdf",
        "results/deeptools/scree.pdf",

    ]

    return TARGETS


def paired_end():
    """Checks if paired-end reads are used
    """
    # Get one fastq file
    reads = glob.glob("reads/*/*.gz")
    assert len(reads) != 0, "No fastq files found..."
    fastq = reads[0]

    # Check file extension to see if paired-end reads are used
    if fastq.endswith(".fastq.gz"):
        return False
    elif fastq.endswith("_R1_001.fastq.gz"):
        return True
    elif fastq.endswith("_R2_001.fastq.gz"):
        return True
    else:
        raise ValueError("ERROR: Could not determine if paired-end reads are used\n"
                         "Please check if fastq files end with either .fastq.gz (SE) or _R1_001.fastq.gz/_R2_001.fastq.gz (PE)")
    
    