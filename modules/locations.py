import os

class Locations:
    def __init__(self, workingDirectory):
        self.workingDirectory = workingDirectory
    
    @property
    def workingDirectory(self):
        return self._workingDirectory
    
    @workingDirectory.setter
    def workingDirectory(self, value):
        value = os.path.abspath(value)
        if not os.path.isdir(value):
            raise FileNotFoundError(f"-d working directory '{value}' does not exist!")
        self._workingDirectory = value
    
    # Naive directory properties
    @property
    def depthDir(self):
        return os.path.join(self.workingDirectory, "depth")
    
    @property
    def callDir(self):
        return os.path.join(self.workingDirectory, "call")
    
    # Naive file properties
    @property
    def depthSuffix(self):
        return ".depth.tsv"
    
    @property
    def finalDeletionFile(self):
        return os.path.join(self.depthDir, "psQTL_deletions.vcf.gz") # Make sure this is gzipped
    
    @property
    def bamListFile(self):
        return os.path.join(self.callDir, "bamlist.txt")
    
    @property
    def vcfFile(self):
        return os.path.join(self.callDir, "psQTL_variants.vcf.gz")
    
    @property
    def filteredVcfFile(self):
        return os.path.join(self.callDir, "psQTL_variants.filtered.vcf.gz")
    
    @property
    def variantEdFile(self):
        return os.path.join(self.workingDirectory, "psQTL_variants.ed.tsv.gz")
    
    @property
    def depthEdFile(self):
        return os.path.join(self.workingDirectory, "psQTL_depth.ed.tsv.gz")
    