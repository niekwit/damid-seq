#from itertools import product
import re
import numpy as np

def targets():
    """
    Returns file targets for rule all
    """
    # Base target files
    TARGETS = [
        "results/qc/multiqc/multiqc.html",
        #"results/plots/PCA.pdf",
        #"results/plots/scree.pdf",
        "results/plots/sample_correlation.pdf",
        "results/plots/heatmap.pdf",
        "results/plots/profile_plot.pdf",
        "results/plots/mapping_rates.pdf",
        expand("results/bigwig_rev_log2/average_bw/{bg_sample}.bw", bg_sample=BG_SAMPLES),
        ]
    
    if config["peak_calling_perl"]["run"]:
        TARGETS.extend([
            "results/plots/peaks/feature_distributions.pdf",
            "results/plots/peaks/distance_to_tss.pdf",
            expand("results/peaks/overlapping_peaks/{bg_sample}.annotated.txt", bg_sample=BG_SAMPLES),
            ])
    
    if config["peak_calling_macs2"]["run"]:
        if config["peak_calling_macs2"]["mode"] == "narrow":
            TARGETS.extend([
                "results/plots/macs2_narrow/feature_distributions.pdf",
                "results/plots/macs2_narrow/distance_to_tss.pdf",
                expand("results/macs2_narrow/{bg_sample}_peaks.xls", bg_sample=BG_SAMPLES),
                expand("results/macs2_narrow/{bg_sample}_peaks.narrowPeak", bg_sample=BG_SAMPLES),
                expand("results/macs2_narrow/{bg_sample}_summits.bed", bg_sample=BG_SAMPLES),
                ])
        elif config["peak_calling_macs2"]["mode"] == "broad":
            TARGETS.extend([
                "results/plots/macs2_broad/broad/feature_distributions.pdf",
                "results/plots/macs2_broad/broad/distance_to_tss.pdf",
                expand("results/macs2_broad/{bg_sample}_peaks.xls", bg_sample=BG_SAMPLES),
                expand("results/macs2_broad/{bg_sample}_peaks.broadPeak", bg_sample=BG_SAMPLES),
                expand("results/macs2_broad/{bg_sample}_peaks.gappedPeak", bg_sample=BG_SAMPLES),
                ])
    
    return TARGETS


def data_type():
    """
    Detect whether Dam only or Dam-POI sample numbers match or not.
    If they do not match, all samples will be in reads/ directory, wihtout any subdirectories. 
    If they do match, each biological replicate will be in a separate subdirectory in reads/.
    """
    # Check if any subdirectories exist in reads/
    subdirs = glob.glob("reads/*/*")
    if len(subdirs) == 0:
        return "matrix"
    else:
        return "paired"


def matrix_samples():
    """
    Match each Dam only sample to all Dam-POI single replicate samples

    Note: only works for paired-end reads at the moment.
    """
    print("Matrix samples detected...")
    print("Preparing directory structure for all combinations of Dam controls vs Dam-fusion(s)...")
    # Get all R1 files in reads/
    r1 = glob.glob("reads/*_R1_001.fastq.gz")

    # Get all R1 Dam only files in reads/
    dam = glob.glob("reads/*Dam*_R1_001.fastq.gz")

    # Get all R1 Dam-POI files in reads/
    fusion = [f for f in r1 if f not in dam]

    # Get base names of all fusion R1 files
    fusion_base_names = list(set([os.path.basename(f).replace("_R1_001.fastq.gz","") for f in fusion]))

    # For each base name match each replicate R1 file into nested list
    base_replicates = []
    for base in fusion_base_names:
        tmp = [f for f in fusion if base in f].sorted()
        base_replicates.append(tmp)

    # Check if each nested list has the same length
    if len(set([len(l) for l in base_replicates])) != 1:
        raise ValueError("Number of Dam-POI replicates does not match...")
    
    # Create nested list with the first sample of each sample from base_replicates (transpose list)
    all_base_replicates = np.array(base_replicates).T.tolist()

    # For each Dam sample, match it to all Dam-POI samples (read 1)
    matched_samples_r1 = []
    for d in dam:
        for f in all_base_replicates:
            f.insert(0, d)
            matched_samples_r1.append(f)
    
    # Get read 2 files
    matched_samples_r2 = []
    for l in matched_samples_r1:
        tmp = [f.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz") for f in l]
        matched_samples_r2.append(tmp)
    
    # Create directory names for each matched sample list
    # Symlink files to these directories
    def symlink_files(list_):
        """
        Symlink each file in list to unique dir_
        Remove _[0-9] before paired-end end extension
        All replicates will have same file name but will be in different directories
        """
        log = []
        for i, l in enumerate(list_):
            dir_ = f"reads/repl_{i}"
            os.makedirs(dir_, exist_ok=True)
            
            # Keep log of what goes where
            tmp = [dir_]
            tmp.extend(l)
            log.append(tmp)
            for f in l:
                dest = re.sub(r"_\d{1,2}_R", "_R", f)  
                os.symlink(f, f"{dir_}/{os.path.basename(dest)}")
        return log

    log_r1 = symlink_files(matched_samples_r1)
    log_r2 = symlink_files(matched_samples_r2)
                
    # Add log data to data frame
    dirs_ = [os.path.basename(d[0]) for d in log_r1]
    df = pd.DataFrame([log_r1, log_r2], columns=dirs_)

    # Save log data to csv
    df.to_csv("reads/sample_matrix.csv", index=False)


def dirs():
    """
    Each dir contains one replicate sets of fastq files
    """
    DIRS = glob.glob("reads/*")
    DIRS = [os.path.basename(d) for d in DIRS]

    if len(DIRS) == 0:
        raise ValueError("No replicate directories found in reads directory...")
        
    return DIRS
    

def samples(bedgraph=False, dam=False):
    """
    Checks sample names/files and returns sample wildcard values for Snakemake
    """
    SAMPLES = csv["sample"].tolist()
    
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
    if bedgraph and not dam:
        SAMPLES = [s for s in SAMPLES if not "Dam" in s]
    elif bedgraph and dam:
        SAMPLES = [s for s in SAMPLES if "Dam" in s]

    return SAMPLES


def paired_end():
    """
    Checks if paired-end reads are used
    """
    # Get one fastq file
    reads = glob.glob("reads/*/*.gz")
    assert len(reads) != 0, "No fastq files found..."
    fastq = reads[0]

    # Check file extension to see if paired-end reads are used
    if fastq.endswith("_R1_001.fastq.gz"):
        return True, ""
    elif fastq.endswith("_R2_001.fastq.gz"):
        return True, ""
    else:
        return False, "-ext300"
    

def computematrix_args(region_labels=False):
    """
    Returns computeMatrix arguments as string based on config file.
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


def masked_genes():
    """
    Returns string with genes that were masked in fasta file
    """
    genes = config["fusion_genes"]

    # If no genes are given, return no_genes
    if genes == "":
        genes = "no_genes"
    else:
        # Make sure there are no spaces
        genes = genes.replace(" ", "")

        # Check if gene names are Ensemble IDs
        for gene in genes.split(","):
            if "hg" in resources.genome:
                prefix="ENSG"
                count="11"
            if "mm" in resources.genome:
                prefix="ENSMUSG"
                count="11"
            if "dm" in resources.genome:
                prefix="FBgn"
                count="7"
            if not re.match("".join([prefix,"[0-9]","{",count,"}"]), gene):
                raise ValueError(f"Gene {gene} is not an Ensemble ID")

        # Replace comma with underscore
        genes = genes.replace(",", "_")
        
        return genes


def macs2_params():
    if paired_end:
        format_ = "BAMPE"
    else:
        format_ = "BAM"
    
    if "hg" in resources.genome:
        genome = "hs"
    elif "mm" in resources.genome:
        genome = "mm"
    elif "dm" in resources.genome:
        genome = "dm"
    
    if config["peak_calling_macs2"]["mode"] == "broad":
        cutoff = config["peak_calling_macs2"]["broad_cutoff"]
        broad = "--broad --broad-cutoff {cutoff} "
    else:
        broad = ""
    
    qvalue = config["peak_calling_macs2"]["qvalue"]
    extra = config["peak_calling_macs2"]["extra"]

    return f"-f {format_} -g {genome} -q {qvalue} {broad} {extra}"