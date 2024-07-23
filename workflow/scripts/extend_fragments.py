#!/usr/bin/env python3

"""
Extends single-end reads to a specified length or the nearest GATC site, 
whichever is closer.
Also filters out reads with low mapping quality.

Approach:
1. Read in GATC sites from a GFF file.
    a. For each chromosome, get the start position of each GATC site.
2. For each read, check if the read is below a certain mapping quality threshold.
3. If the read is below the threshold, skip to the next read.
4. Get the chromosome, start and end position of the read.
5. Get the closest GATC site to the end of the read (omit the ones that came before the end of the read).
6. If the distance from the end of the read to the closest GATC site is less than n,
    set the length to add to the read to the distance.
7. If the distance is greater than n, set the length to add to n.
8. Update the CIGAR string of the read to reflect the new length.
9. Create a string for the new fragment in SAM format.
10. Write the new fragment to a new BAM file.
"""

import logging
import re
import subprocess
import pysam
import numpy as np
import gffpandas.gffpandas as gffpd

input_bam = snakemake.input[0]
output_bam = snakemake.output[0]
#output_sam = output_bam.replace(".bam", ".sam")
gatc_gff = snakemake.input["gatc_gff"]
n = snakemake.params["n"]
q = snakemake.params["q"]
log = snakemake.log[0]
threads = snakemake.threads # for later

# Setup log file
logging.basicConfig(format='%(levelname)s:%(message)s', 
                    level=logging.DEBUG,
                    handlers=[logging.FileHandler(log),
                              logging.StreamHandler()])

# For each chromosome obtain start sites of GATC sites in np array
logging.info("Reading GATC sites from GFF file...")
gatc_gff = gffpd.read_gff3(gatc_gff).df

chr_gatc = {}
for chr in gatc_gff["seq_id"].unique():
    chr_gatc[chr] = gatc_gff[gatc_gff["seq_id"] == chr]
    
    # Only keep the start position
    chr_gatc[chr] = chr_gatc[chr]["start"].values + 2 # Correct for DpnI cut site
    
    # Convert to Numpy array for speed
    chr_gatc[chr] = np.sort(chr_gatc[chr])

# Counter for number of sequences MQ filtered/extended 
# and total number of sequences
extended_seqs = 0
all_seqs = 0

# idea for parallisation:
# each worker can do a subset of the chromosomes

logging.info(f"Extending reads up to {n} bases or nearest GATC site...")
samfile_in = pysam.AlignmentFile(input_bam, "rb")
#samfile_out = pysam.AlignmentFile(output_sam, "wb", template=samfile)

for read in samfile_in.fetch():
    print(read)
    # First filter out reads with low mapping quality
    mq = read.mapping_quality
    if mq < q:
        all_seqs += 1
        continue
    
    # Get read information
    chr_ = read.reference_name
    start = read.reference_start
    length = read.query_length
    print(length)
    end = start + length
    
    # Get the GATC sites for the current chromosome
    tmp_gff = chr_gatc[chr_]
    
    # Remove GATC sites that are before the read end
    tmp_gff = tmp_gff[tmp_gff >= end]
    
    # Find index of GATC site closest to the read end
    closest = np.argmin(abs(end - tmp_gff))
    
    # Check if distance from read start to closest GATC site is less than n
    if tmp_gff[closest] - start < n:
        # Get distance from fragment end to closest GATC site
        # This is the length that needs to be added to the read
        length_to_add = tmp_gff[closest] - end
    else:
        length_to_add = n - length
        
    new_sequence_length = length + length_to_add
        
    # Update CIGAR string for the new read length
    cigar = read.cigarstring
    if re.search(r"^\d+M$", cigar):
        cigar = f"{length + length_to_add}M"
    elif re.search(r"\d+M$", cigar):
        # Get the length of the \d+M part at the end
        l = int(re.search(r"\d+M$", cigar).group().replace("M", ""))
        
        # Add the length to add to CIGAR string
        cigar = re.sub(r"\d+M$", f"{l + length_to_add}M", cigar)
    else:
        cigar = f"{cigar}{length_to_add}M"
    
    """
    # Update quality scores
    old_qual = read.query_alignment_qualities
    
    # Create updated read and write to file
    new_read = pysam.AlignedSegment()
    new_read.query_name = read.query_name
    new_read.query_sequence = "*"
    new_read.flag = read.flag
    new_read.reference_id = read.reference_id
    new_read.reference_start = read.reference_start
    new_read.mapping_quality = read.mapping_quality
    new_read.cigarstring = cigar
    new_read.next_reference_id = read.next_reference_id
    new_read.next_reference_start = read.next_reference_start
    new_read.template_length = length #new_sequence_length
    new_read.query_qualities = pysam.qualitystring_to_array("*")
    new_read.tags = read.tags
    #print(new_read)
    """
    
    read_out = "\t".join([
        read.query_name,
        str(read.flag),
        str(read.reference_id),
        str(start),
        str(mq),
        cigar,
        "*",
        "0",
        "0",
        "*",
        "*",
        read.tags
        ])
    print(read_out)
    # Write to BAM file
    command = f"samtools view -Shb - > {output_bam}"
    process = subprocess.Popen(command, stdin=subprocess.PIPE, shell=True)
    process.stdin.write(bf"{read_out}\n")
    
    all_seqs += 1
    extended_seqs += 1
    break

process.stdin.close()
process.wait()
#samfile_in.close()
#samfile_out.close()

# Convert to BAM
#pysam.view("-@", threads, "-o", output_bam, samfile_out)

logging.info(f"Total sequences: {all_seqs}")
logging.info(f"Sequences with MQ < {q}: {all_seqs - extended_seqs}")
logging.info(f"Sequences extended: {extended_seqs}")