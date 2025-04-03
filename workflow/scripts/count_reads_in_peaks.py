import pybedtools as pbt
import pysam
import pandas as pd
import os

# Get Snakemake variables
bams = snakemake.input["bams"]
peaks = snakemake.input["peaks"]

# Create data frame to store name, total counts,
# counts in peaks, and fraction of reads in peaks
df = pd.DataFrame(
    columns=["name", "total_reads", "reads_in_peaks", "fraction_in_peaks"]
)

# Load peaks into bedtools object
for i, bam in enumerate(bams):
    # Get total reads in bam file
    with pysam.AlignmentFile(bam, "rb") as bamfile:
        num_alignments = bamfile.count()

    # Convert bam to bedtools object
    bed = pbt.BedTool(bam).bamtobed()

    # Intersect with peaks
    intersect = bed.intersect(peaks[i], u=True)

    # Get number of intersecting reads
    num_intersecting = intersect.count()

    # Calculate fraction of reads in peaks
    fraction_in_peaks = num_intersecting / num_alignments

    name = os.path.join(os.path.basename(os.path.dirname(bam)), os.path.basename(bam))
    name = name.replace(".sorted.bam", "")

    # Append to data frame
    df = pd.concat(
        [
            df,
            pd.DataFrame(
                [
                    {
                        "name": name,
                        "total_reads": num_alignments,
                        "reads_in_peaks": num_intersecting,
                        "fraction_in_peaks": fraction_in_peaks,
                    }
                ]
            ),
        ],
        ignore_index=True,
    )

# Save to file
df.to_csv(snakemake.output["csv"], index=False)
