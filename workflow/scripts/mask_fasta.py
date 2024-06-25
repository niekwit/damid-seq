import sys
import subprocess
from Bio import SeqIO

"""
Replaces gene sequences set in config:fusion_genes:genes 
in fasta file with Ns. Either whole gene sequence or
just exons can be masked (config:fusion_genes:feature_to_mask).

Reason: plasmid expressing Dam fusion 
genes can be methylated at very high levels
"""

# Load Snakemake variables
gtf = snakemake.input["gtf"]
genes2mask = snakemake.params["g2m"]
feature2mask = snakemake.params["f2m"]
fasta = snakemake.input["fa"]
masked_fasta = snakemake.output["out"]

def write_dict2fasta(d, out):
    """
    Write dict to fasta file
    """
    with open(out, "w") as f:
        for name, seq in d.items():
            f.write(f">{name}\n{seq}\n")

# Load fasta file as dict
chr_seq = {}
for chr_ in SeqIO.parse(fasta, "fasta"):
    chr_seq[chr_.id] = chr_.seq

# Mask gene sequences with Ns
if genes2mask == "no_genes":
    print(f"No genes to mask from {fasta}...")
    
    # Write unmasked fasta as masked fasta file
    write_dict2fasta(chr_seq, masked_fasta)
else:
    for gene in genes2mask.split("_"):
        print(f"Masking {gene} sequence from {fasta} (feature {feature2mask})...")
                
        # Get genomic coordinates of selected feature of gene to mask from GTF file
        cmd = f"""sed '1,4d' {gtf} | awk '{{if ($3 == "{feature2mask}") {{print $0}} }}' | grep {gene}"""
        try:
            lines = subprocess.check_output(cmd, shell=True).decode().split("\n")
        except subprocess.CalledProcessError:
            print(f"Gene {gene} not found in {gtf}...")
            sys.exit(1)
        
        # Mask each feature of gene with Ns
        for line in lines:
            try:
                chr, db, t, start, end, *args = line.split("\t")
            except ValueError:
                continue # Skip empty line (last one)
            
            # Load chromosome sequence where gene feuture is located
            seq = chr_seq[chr]
            
            # Correct start and end positions for 0-based indexing
            start = int(start) - 1
            end = int(end) - 1
            
            # Mask gene feature sequence with Ns
            seq_masked = seq[:start] + "N" * (end - start) + seq[end:]
            
            # Replace sequence in dict
            chr_seq[chr] = seq_masked

    # Write masked fasta to file
    print(f"Writing masked sequence(s) to {masked_fasta}...")
    write_dict2fasta(chr_seq, masked_fasta)