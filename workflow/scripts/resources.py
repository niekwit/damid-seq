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
                        
            # set sha256sums for unzipped genome files
            self.fasta_sha256 = "1e74081a49ceb9739cc14c812fbb8b3db978eb80ba8e5350beb80d8ad8dfef3b"
            self.gtf_sha256 = "12582b0db02ebe19c29c5733c6edaa62599fe934af593cb7f24423a14db3186c"
                      
        elif "mm" in genome:
            if genome == "mm9":
                name = "GRCm38"
            elif genome == "mm10":
                name = "GCRm39"
                
            # create URLs for genome files
            self.fasta_url = f"{base_url_ens}fasta/mus_musculus/dna/Mus_musculus.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = f"{base_url_ens}gtf/mus_musculus/Mus_musculus.{name}.{build}.gtf.gz"
            
            # set sha256sums for unzipped genome files
            self.fasta_sha256 = "14571f7559e292baf0a40f9d155c41ede19a04d80fdeb59a0c2dfe566db90552"
            self.gtf_sha256 = "6efbe1fdbd41d4321daf6d550db240656473b41a107648d6faaf9d61cfdb6c4d"
            
        # downloaded unzipped file names
        self.fasta = self._file_from_url(self.fasta_url)
        self.gtf = self._file_from_url(self.gtf_url)
        
    def _file_from_url(self, url):
        """Returns file path for unzipped downloaded file
        """
        
        return f"resources/{os.path.basename(url).replace('.gz','')}"
    
    
  
        
        
            
            
    