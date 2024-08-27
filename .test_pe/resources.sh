# Create fasta file with just chromosome 11
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa -r regions.txt | pigz > Homo_sapiens.GRCh38.dna.primary_assembly_subset_chr11.fa.gz