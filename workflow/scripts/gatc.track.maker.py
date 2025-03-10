#!/usr/bin/env python3

import re
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from Bio import SeqIO

genome = snakemake.params["genome"]
motif = snakemake.params["motif"]
fasta = snakemake.input["fa"]
outfile = snakemake.output["out"]
threads = snakemake.threads

# Setup logging
log = snakemake.log[0]
logging.basicConfig(
    format="%(levelname)s:%(message)s",
    level=logging.DEBUG,
    handlers=[logging.FileHandler(log)],
)


def load_fasta(fasta):
    """
    Load fasta file and return dict with chr as key and sequence as value.
    Omit scaffold and mitochondrial sequences if specified
    """
    chr_seq = {}
    for chr_ in SeqIO.parse(fasta, "fasta"):
        chr_seq[chr_.id] = chr_.seq

    return chr_seq


# Generate track function
def generate_track(motif, fasta, gff, threads):
    motif = motif
    motif_len = len(motif)

    fasta = fasta
    logging.info(f"Opening {fasta}...")
    chr_seq = load_fasta(fasta)

    gff = outfile
    logging.info(f"Writing data to {gff}...")

    with open(gff, "w") as track:

        chromosomes = chr_seq.keys()

        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = []
            for chromosome in chromosomes:
                logging.info(f"Processing chromosome: {chromosome}")
                seq_str = str(chr_seq[chromosome])
                futures.append(
                    executor.submit(process, chromosome, seq_str, motif, motif_len)
                )
                logging.info(f"{chromosome} done!")

            for future in as_completed(futures):
                chromosome, results = future.result()
                for result in results:
                    track.write(result)


# Process function
def process(chr_name, seq, motif, motif_len):
    results = motif_hash(seq, chr_name, motif, motif_len)
    return chr_name, results


# Motif hash function
def motif_hash(seq, chr_name, motif, motif_len):
    results = []
    seq_len = len(seq)
    for base in range(seq_len - motif_len + 1):
        window = seq[base : base + motif_len]
        if re.search(motif, window, re.IGNORECASE):
            results.append(
                f"{chr_name}\t.\t.\t{base}\t{base + motif_len}\t1\t+\t.\t.\n"
            )
    return results


logging.info(f"Mapping {motif} sites in {fasta}")
generate_track(motif, fasta, outfile, threads)
logging.info("All done!")
