import convert_bed2fasta as utils
import subprocess

"""
Replaces gene sequences set in config:fusion_genes 
in fasta file with Ns.

Reason: plasmid expressing Dam fusion 
genes can be methylated at very high levels
"""

# Load Snakemake variables
gtf = snakemake.input["gtf"]
genes2mask = snakemake.params["g2m"]
fasta = snakemake.input["fa"]
masked_fasta = snakemake.output["out"]

# Load fasta file as dict
chr_seq = utils.load_fasta(fasta)

if genes2mask == "no_genes":
    print(f"No genes to mask from {fasta}...")
    
    # Write unmasked fasta as masked fasta file
    utils.write_dict2fasta(chr_seq, masked_fasta)
else:
    for gene in genes2mask.split("_"):
        print(f"Masking {gene} sequence from {fasta}...")
        
        # Get genomic coordinates of genes to mask from GTF file
        cmd = f"""sed '1,4d' {gtf} | awk '{{if ($3 == "gene") {{print $0}} }}' | grep {gene}"""
        try:
            line = subprocess.check_output(cmd, shell=True).decode()
        except subprocess.CalledProcessError:
            print(f"ERROR: gene {gene} not found in {gtf}...")
            sys.exit(1)
        chr, db, t, start, end, *args = line.split("\t")
        
        # Load chromosome sequence where gene is located
        seq = chr_seq[chr]
        
        # Correct start and end positions for 0-based indexing
        start = int(start) - 1
        end = int(end) - 1
        
        # Mask gene sequence with Ns
        seq_masked = seq[:start] + "N" * (end - start) + seq[end:]
        
        # Replace sequence in dict
        chr_seq[chr] = seq_masked

    # Write masked fasta to file
    print(f"Writing masked sequence(s) to {masked_fasta}...")
    utils.write_dict2fasta(chr_seq, masked_fasta)