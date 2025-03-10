"""
Converts a bed file to a fasta file. 
As DamID is not strand specific, the sequence 
will be extracted from the positive strand only.
"""
import pandas as pd
from Bio import SeqIO


def load_bed(bed):
    """
    Load bed file as pandas df
    """
    # Load bed file as pandas df and add header names
    # This is not the most efficient way to do this (optimise later)
    bed = pd.read_csv(bed, sep="\t", header=None)
    bed.columns = ["chr", "start", "end", "name", "strand", "score"]

    return bed


def load_fasta(fasta):
    """
    Load fasta file and return dict with chr as key and sequence as value
    """
    chr_seq = {}
    for chr_ in SeqIO.parse(fasta, "fasta"):
        chr_seq[chr_.id] = chr_.seq

    return chr_seq


def write_dict2fasta(d, out):
    """
    Write dict to fasta file
    """
    with open(out, "w") as f:
        for name, seq in d.items():
            f.write(f">{name}\n{seq}\n")


def main(bed, chr_seq):
    # Dict to store sequence per peak
    peak_seq = {}

    # For each region in bed file, extract sequence and write to fasta file
    for row in bed.itertuples(index=False):
        chr = row.chr
        start = row.start
        end = row.end
        name = row.name

        # Extract sequence
        seq = chr_seq[chr][start:end]

        # Store sequence in dict
        peak_seq[f"{name}_{chr}_{start}_{end}"] = seq

    # Write to fasta file
    write_dict2fasta(peak_seq, snakemake.output["out"])


if __name__ == "__main__":
    bed = load_bed(snakemake.input["bed"])
    chr_seq = load_fasta(snakemake.input["fasta"])

    main(bed, chr_seq)
