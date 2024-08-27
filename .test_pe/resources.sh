# Create fasta file with regions were just peaks are located
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa -r regions.txt | sed 's/:.*$//' |  pigz > Homo_sapiens.GRCh38.dna.primary_assembly_subset.fa.gz

