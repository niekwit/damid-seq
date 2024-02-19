import random
import convert_bed2fasta as utils
import pandas as pd

bed = utils.load_bed(snakemake.input["bed"])
fasta = utils.load_fasta(snakemake.input["fasta"]) # whole genome fasta
chrom_sizes = snakemake.input["cs"]
chrom_sizes = pd.read_csv(chrom_sizes, 
                          sep="\t", 
                          header=None,
                          names=["chr", "size"])

# Remove non-canonical chromosomes (beginning with GL, KI) from chrom_sizes
chrom_sizes = chrom_sizes[~chrom_sizes["chr"].str.contains("GL|KI")]

# Add weight to each chromosome based on size
chrom_sizes["weight"] = chrom_sizes["size"] / chrom_sizes["size"].sum()

# Get chromosome names from fasta dict
chr_names = chrom_sizes["chr"].tolist()

# Get value lengths of each key in fasta dict
chr_lengths = chrom_sizes["size"].tolist()

# Create dict with chromosome names as keys and chromosome lengths as values
chr_seq = dict(zip(chr_names, chr_lengths))

def random_seq(length):
    """
    Returns random genomic sequence with no Ns in sequence
    Length of sequence is the same as the input peaks
    """
    while True:
        # Randomly select chromosome
        chr = random.choices(population=chrom_sizes["chr"],
                             weights=chrom_sizes["weight"],
                             k=1)[0]
        
        # Get chromosome length
        chr_length = chr_seq[chr]
        
        # Randomly select end and start position
        while True: # Make sure that no negative start positions are generated
            end = random.randint(1, chr_length)
            if end - length > 0:
                break
        
        start = end - length
        
        # Extract sequence
        seq = fasta[chr][start:end]
        
        # Check if sequence contains any Ns, if so pick another sequence
        if "N" not in seq:
            return seq, chr, start, end

# Dict to store random regions
random_regions = {}

# For each peak in bed file, create a random sequence with the same length
for row in bed.itertuples(index=False):
    # Get peak length
    start = int(row.start)
    end = int(row.end)
    name = row.name
    peak_number = name.split("_")[1]
    length = end - start

    # Get random sequence, chromosome, start position and end position
    _seq, _chr, _start, _end = random_seq(length)
        
    # Store sequence in dict
    random_regions[f"random_{peak_number}_{_chr}_{_start}_{_end}"] = _seq

# Write random regions to fasta file
utils.write_dict2fasta(random_regions, snakemake.output["out"])

