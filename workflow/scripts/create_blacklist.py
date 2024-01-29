"""
Based on parsed genes names, create a bed file of gene locations of these genes.
The plasmid used to express the Dam-Protein fusion can give a lot of signal, so
we want to remove these loci from the analysis.
"""

from snakemake.shell import shell

# Load Snakemake variables
genes = snakemake.params["genes"]
genes = genes.split(",")
gtf = snakemake.input["gtf"]
txt = snakemake.output["txt"]
bed = snakemake.output["bed"]

print(f"Creating blacklist file (.bed) for loci of {','.join(genes)}")

# Create file with genes names to search for (input for grep)
with open(txt, "a") as f:
  for gene in genes:
    print(gene, file=f)

# Create bed file with gene locations
# Filter for protein_coding, otherwise potential AS, pseudogenes, etc. will be included
shell(
    """
    awk '{{if ($3 == "gene") print $0}}' {gtf} | grep protein_coding | grep -wf genes.txt | gff2bed > {bed}
    """
)

# Check whether bed file contains as many lines as genes requested
# If not, print warning/error
num_genes = len(genes)
with open(bed) as f:
   num_lines = sum(1 for _ in f)

if num_lines < num_genes: # Warning
    print("Warning: number of genes in bed file is less than number of genes requested.\nCheck if genes are present in {gtf}.")
elif num_lines > num_genes: # Should not happen
    raise ValueError("Number of genes in bed file is greater than number of genes requested.")
