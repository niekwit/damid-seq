import os

class Resources:
    """Gets URLs and file names of fasta and GTF files for a given genome and build
    """
    # create genome directory
    os.makedirs("resources/", exist_ok=True)
    
    def __init__(self, genome, build):
        self.genome = genome
        self.build = build
                
        # base URLs
        base_url_ens = f"https://ftp.ensembl.org/pub/release-{build}/"
                
        if "hg" in genome:
            if genome == "hg19":
                name = "GRCh37"
            elif genome == "hg38":
                name = "GRCh38"
                
            # create URLs for genome files
            self.fasta_url = f"{base_url_ens}fasta/homo_sapiens/dna/Homo_sapiens.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = f"{base_url_ens}gtf/homo_sapiens/Homo_sapiens.{name}.{build}.gtf.gz"
                      
        elif "mm" in genome:
            if genome == "mm38":
                name = "GRCm38"
            elif genome == "mm39":
                name = "GRCm39"
                
            # create URLs for genome files
            self.fasta_url = f"{base_url_ens}fasta/mus_musculus/dna/Mus_musculus.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = f"{base_url_ens}gtf/mus_musculus/Mus_musculus.{name}.{build}.gtf.gz"        
        
        elif "dm" in genome:
            if genome == "dm6":
                name = "BDGP6.46"

            self.fasta_url = f"{base_url_ens}fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.{name}.dna.toplevel.fa.gz"
            self.gtf_url = f"{base_url_ens}gtf/drosophila_melanogaster/Drosophila_melanogaster.{name}.{build}.gtf.gz" 
        
        elif "test" in genome:
            if genome == "dm6":
                name = "BDGP6.46"

            self.fasta_url = "https://github.com/niekwit/damid-seq/raw/main/.test_pe/Homo_sapiens.GRCh38.dna.primary_assembly_subset.fa.gz"
            self.gtf_url = "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
            
        else:
            raise ValueError("Genome {genome} not supported")
        
        # downloaded unzipped file names
        self.fasta = self._file_from_url(self.fasta_url)
        self.gtf = self._file_from_url(self.gtf_url)
        
        
    def _file_from_url(self, url):
        """Returns file path for unzipped downloaded file
        """
        return f"resources/{os.path.basename(url).replace('.gz','')}"
    
    