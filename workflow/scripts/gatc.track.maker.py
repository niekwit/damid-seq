#!/usr/bin/env python3

import argparse
import gzip
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO

# Command line argument processing
def process_cli():
    parser = argparse.ArgumentParser(description='Generates a GFF file containing the locations of all GATC sites in the genome sequence')
    parser.add_argument("--fasta", "-f",
                        type=str, 
                        help='Input genome FASTA file')
    parser.add_argument("--outfile", "-o", 
                        type=str, 
                        required=True,
                        help="Output file name")
    parser.add_argument("--genome", "-g",
                        type=str,
                        default="hg38",
                        choices=["hg38", "mm38", "dm6", 
                                 "dm3", "hg19", "mm39"],
                        required=True,
                        help='Genome build (default: hg38)')
    parser.add_argument("--threads", "-t", 
                        type=int, 
                        default=1, 
                        help="Number of parallel processes")
    parser.add_argument("--motif", "-m", 
                        type=str, 
                        default="GATC", 
                        help="Motif to search for (default: GATC)")
    parser.add_argument("--scaffolds", 
                        action="store_true", 
                        help="Process scaffold assemblies (not recommended)")
    parser.add_argument("--mito", 
                        action="store_true", 
                        help="Process mitochondrial chromosome (not recommended)")
    return parser.parse_args()


def regex_patterns(genome):
    # Dictionary to store scaffold and mitochondrial regex patterns
    patterns = {
        "hg38": [r"",""],
        "hg19": [r"",""],
        "mm38": [r"",""],
        "mm39": [r"",""],
        "dm6": [r"Scaffold|\d{15}","mitochondrion_genome"],
        "dm3": [r"",""], 
    }
    return patterns[genome]


def load_fasta(fasta, args):
    """
    Load fasta file and return dict with chr as key and sequence as value.
    Omit scaffold and mitochondrial sequences if specified
    """
    chr_seq = {}
    for chr_ in SeqIO.parse(fasta, "fasta"):
        chr_seq[chr_.id] = chr_.seq
    
    scaffold_pattern, mito_pattern = regex_patterns(args.genome)
    
    # Remove scaffold chromosomes and mitochondrial chromosome if specified
    if not args.scaffolds:
        chr_seq = {k: v for k, v in chr_seq.items() if not re.search(scaffold_pattern, k, re.IGNORECASE)}
    if not args.mito:
        chr_seq = {k: v for k, v in chr_seq.items() if not re.search(mito_pattern, k, re.IGNORECASE)}
    
    return chr_seq
   
    
# Generate track function
def generate_track(args):
    motif = "GATC"
    motif_len = len(motif)
    
    fasta = args.fasta
    print(f"Opening {fasta}...")
    chr_seq = load_fasta(fasta, args)
    
    gff = args.outfile
    print(f"Writing data to {gff}...")

    with open(gff, "w") as track:
              
        chromosomes = chr_seq.keys()
         
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            futures = []
            for chromosome in chromosomes:
                print(f"Processing chromosome: {chromosome}")
                seq_str = str(chr_seq[chromosome])
                futures.append(executor.submit(process, chromosome, seq_str, motif, motif_len))

            for future in as_completed(futures):
                chromosome, results = future.result()
                for result in results:
                    track.write(result)

    print("All done!")


# Process function
def process(chr_name, seq, motif, motif_len):
    #print(f"Processing {chr_name} ...", file=sys.stderr)
    results = motif_hash(seq, chr_name, motif, motif_len)
    return chr_name, results


# Motif hash function
def motif_hash(seq, chr_name, motif, motif_len):
    results = []
    seq_len = len(seq)
    for base in range(seq_len - motif_len + 1):
        window = seq[base:base + motif_len]
        if re.search(motif, window, re.IGNORECASE):
            results.append(f"{chr_name}\t.\t.\t{base}\t{base + motif_len}\t1\t+\t.\t.\n")
    print(f"{chr_name} done!")
    return results

# Main function
def main():
    args = process_cli()
    generate_track(args)

if __name__ == "__main__":
    main()