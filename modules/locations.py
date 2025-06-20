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
    
    @property
    def splsdaDir(self):
        return os.path.join(self.workingDirectory, "splsda")
    
    # Naive file properties
    @property
    def depthSuffix(self):
        return ".depth.tsv"
    
    @property
    def finalDepthFile(self):
        return os.path.join(self.depthDir, "psQTL_depth.vcf.gz")
    
    @property
    def bamListFile(self):
        return os.path.join(self.callDir, "bamlist.txt")
    
    @property
    def vcfFile(self):
        return os.path.join(self.callDir, "psQTL_call.vcf.gz")
    
    @property
    def filteredVcfFile(self):
        return os.path.join(self.callDir, "psQTL_call.filtered.vcf.gz")
    
    @property
    def allelesEdFile(self):
        return os.path.join(self.workingDirectory, "psQTL_call.alleles_ed.tsv.gz")
    
    @property
    def inheritanceEdFile(self):
        return os.path.join(self.workingDirectory, "psQTL_call.inheritance_ed.tsv.gz")
    
    @property
    def genotypesEdFile(self):
        return os.path.join(self.workingDirectory, "psQTL_call.genotypes_ed.tsv.gz")
    
    @property
    def variantRecodedFile(self):
        return os.path.join(self.splsdaDir, "psQTL_call.recode.tsv.gz")
    
    @property
    def variantSplsdaSelectedFile(self):
        return os.path.join(self.splsdaDir, "psQTL_call.selected.tsv")
    
    @property
    def variantSplsdaBerFile(self):
        return os.path.join(self.splsdaDir, "psQTL_call.BER.tsv")
    
    @property
    def variantSplsdaRdataFile(self):
        return os.path.join(self.splsdaDir, "psQTL_call.Rdata")
    
    @property
    def depthEdFile(self):
        return os.path.join(self.workingDirectory, "psQTL_depth.ed.tsv.gz")
    
    @property
    def depthRecodedFile(self):
        return os.path.join(self.splsdaDir, "psQTL_depth.recode.tsv.gz")
    
    @property
    def depthSplsdaSelectedFile(self):
        return os.path.join(self.splsdaDir, "psQTL_depth.selected.tsv")
    
    @property
    def depthSplsdaBerFile(self):
        return os.path.join(self.splsdaDir, "psQTL_depth.BER.tsv")
    
    @property
    def depthSplsdaRdataFile(self):
        return os.path.join(self.splsdaDir, "psQTL_depth.Rdata")
    
    @property
    def integrativeSplsdaSelectedFile(self):
        return os.path.join(self.splsdaDir, "psQTL_integrative.selected.tsv")
    
    @property
    def windowedSplsdaRscript(self):
        return os.path.join(os.path.dirname(os.path.dirname(__file__)), "utilities", "windowed_splsda.R")
    
    @property
    def integrativeSplsdaRscript(self):
        return os.path.join(os.path.dirname(os.path.dirname(__file__)), "utilities", "integrative_splsda.R")
    
    # Attributes with value input
    def allelesEdPickleFile(self, value):
        return os.path.join(self.workingDirectory, f"psQTL_call.alleles_ed.{value}.pkl")
    
    def inheritanceEdPickleFile(self, value):
        return os.path.join(self.workingDirectory, f"psQTL_call.inheritance_ed.{value}.pkl")
    
    def genotypesEdPickleFile(self, value):
        return os.path.join(self.workingDirectory, f"psQTL_call.genotypes_ed.{value}.pkl")
    
    def depthEdPickleFile(self, value):
        return os.path.join(self.workingDirectory, f"psQTL_depth.ed.{value}.pkl")
