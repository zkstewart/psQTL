import os, json

from .parsing import parse_vcf_stats, parse_deletion_stats, parse_metadata

class ParameterCache:
    def __init__(self, workingDirectory):
        if not os.path.isdir(workingDirectory):
            raise FileNotFoundError(f"Directory '{workingDirectory}' is not a directory.")
        self.workingDirectory = os.path.abspath(workingDirectory)
        
        self._metadataFile = None
        self._vcfFile = None
        self._filteredVcfFile = None
        self._deletionFile = None
        self._bamSuffix = None
        self._bamFiles = None
        self._windowSize = None
        self._qualFilter = None
        self._missingFilter = None
    
    def merge(self, args):
        '''
        Takes the values in the parameter cache and merges them into the argparse object
        if the argparse object has empty values for the cacheable parameters. In practice,
        this allows the 'initialise' submodule of psQTL_prep to be run, with subsequent
        submodules using the cache values if the user does not specify them in the command
        line.
        
        Also merges the argparse values into the cache, so that the cache is kept up-to-date
        with the most recently provided/updated values.
        '''
        self.load()
        for _param, cacheValue in self.__dict__.items():
            param = _param.strip("_") # remove leading underscore from property names
            # Implant cache values into argparse object if key does not exist
            if not hasattr(args, param):
                args.__dict__[param] = cacheValue
            # Merge cache values into argparse object
            elif args.__dict__[param] == None or args.__dict__[param] == []: # only use existing cache values if args are empty
                if cacheValue != None or cacheValue != []: # only use cache values if they are _not_ empty
                    setattr(args, param, self.__dict__[_param])
            # Overrule cache values with argparse values
            elif args.__dict__[param] != self.__dict__[_param]:
                setattr(self, param, args.__dict__[param])
        self.save()
    
    def initialise(self, args):
        '''
        Parameters:
            args -- an argparse.Namespace object containing attributes stored in the cache
        '''
        if os.path.exists(self.cacheFile):
            raise FileExistsError(f"Directory '{self.workingDirectory}' has already been initialised.")
        
        self.metadataFile = args.metadataFile if hasattr(args, "metadataFile") else None
        self.vcfFile = args.vcfFile if hasattr(args, "vcfFile") else None
        self.filteredVcfFile = args.filteredVcfFile if hasattr(args, "filteredVcfFile") else None
        self.deletionFile = args.deletionFile if hasattr(args, "deletionFile") else None
        self.bamSuffix = args.bamSuffix if hasattr(args, "bamSuffix") else None
        self.bamFiles = args.bamFiles if hasattr(args, "bamFiles") else None
        self.windowSize = args.windowSize if hasattr(args, "windowSize") else None
        self.qualFilter = args.qualFilter if hasattr(args, "qualFilter") else None
        self.missingFilter = args.missingFilter if hasattr(args, "missingFilter") else None
        self.save()
    
    def load(self):
        if os.path.exists(self.cacheFile):
            with open(self.cacheFile, "r") as fileIn:
                data = json.load(fileIn)
            # Load data while circumventing the setters
            "Circumvent the setters to avoid re-validation"
            self._metadataFile = data["metadataFile"]
            self._vcfFile = data["vcfFile"]
            self._filteredVcfFile = data["filteredVcfFile"]
            self._deletionFile = data["deletionFile"]
            self._bamSuffix = data["bamSuffix"]
            self._bamFiles = data["bamFiles"]
            self._windowSize = data["windowSize"]
            self._qualFilter = data["qualFilter"]
            self._missingFilter = data["missingFilter"]
            # Now set values using the setter
            "This is necessary to ensure that the values are validated without the comparison being to None"
            self.metadataFile = data["metadataFile"]
            self.vcfFile = data["vcfFile"]
            self.filteredVcfFile = data["filteredVcfFile"]
            self.deletionFile = data["deletionFile"]
            self.bamSuffix = data["bamSuffix"]
            self.bamFiles = data["bamFiles"]
            self.windowSize = data["windowSize"]
            self.qualFilter = data["qualFilter"]
            self.missingFilter = data["missingFilter"]
        else:
            raise FileNotFoundError(f"Working directory '{self.workingDirectory}' has not been initialised.")
    
    def save(self):
        data = {
            "metadataFile": self.metadataFile,
            "vcfFile": self.vcfFile,
            "filteredVcfFile": self.filteredVcfFile,
            "deletionFile": self.deletionFile,
            "bamSuffix": self.bamSuffix,
            "bamFiles": self.bamFiles,
            "windowSize": self.windowSize,
            "qualFilter": self.qualFilter,
            "missingFilter": self.missingFilter
        }
        with open(self.cacheFile, "w") as fileOut:
            json.dump(data, fileOut)
    
    @property
    def cacheFile(self):
        return os.path.join(self.workingDirectory, "parameters.json")
    
    @property
    def metadataFile(self):
        return self._metadataFile
    
    @metadataFile.setter
    def metadataFile(self, value):
        # Ignore unchanged values
        if value == self.metadataFile:
            return
        # Validate file existence
        if value != None:
            if not os.path.isfile(value):
                raise FileNotFoundError(f"Metadata file '{value}' is not a file.")
            value = os.path.abspath(value) # store absolute path
            # Ignore again if value is unchanged after getting absolute path
            if value == self.metadataFile:
                return
        # Store updated value
        updateMsg = f"# Parameter cache: 'metadataFile' changed from '{self._metadataFile}' to '{value}'"
        self._metadataFile = value
        # Propagate change to metadata cache
        metadataCache = MetadataCache(self.workingDirectory)
        metadataCache.establish()
        metadataCache.metadataFile = value
        # Save after metadata cache is updated
        self.save() # try to keep the files in sync by allowing errors to be raised there first
        print(updateMsg)
    
    @property
    def vcfFile(self):
        return self._vcfFile
    
    @vcfFile.setter
    def vcfFile(self, value):
        # Ignore unchanged values
        if value == self.vcfFile:
            return
        # Validate file existence
        if value != None:
            if not os.path.isfile(value):
                raise FileNotFoundError(f"VCF file '{value}' is not a file.")
            value = os.path.abspath(value) # store absolute path
            # Ignore again if value is unchanged after getting absolute path
            if value == self.vcfFile:
                return
        # Store updated value
        updateMsg = f"# Parameter cache: 'vcfFile' changed from '{self._vcfFile}' to '{value}'"
        self._vcfFile = value
        # Propagate change to VCF cache
        vcfCache = VcfCache(self.workingDirectory)
        vcfCache.establish()
        vcfCache.vcfFile = value
        # Save after VCF cache is updated
        self.save() # try to keep the files in sync by allowing errors to be raised there first
        print(updateMsg)
    
    @property
    def filteredVcfFile(self):
        return self._filteredVcfFile
    
    @filteredVcfFile.setter
    def filteredVcfFile(self, value):
        # Ignore unchanged values
        if value == self.filteredVcfFile:
            return
        # Validate file existence
        if value != None:
            if not os.path.isfile(value):
                raise FileNotFoundError(f"Filtered VCF file '{value}' is not a file.")
            value = os.path.abspath(value) # store absolute path
            # Ignore again if value is unchanged after getting absolute path
            if value == self.filteredVcfFile:
                return
        # Store updated value
        updateMsg = f"# Parameter cache: 'filteredVcfFile' changed from '{self._filteredVcfFile}' to '{value}'"
        self._filteredVcfFile = value
        # Propagate change to VCF cache
        vcfCache = VcfCache(self.workingDirectory)
        vcfCache.establish()
        vcfCache.filteredVcfFile = value
        # Save after VCF cache is updated
        self.save() # try to keep the files in sync by allowing errors to be raised there first
        print(updateMsg)
    
    @property
    def deletionFile(self):
        return self._deletionFile
    
    @deletionFile.setter
    def deletionFile(self, value):
        # Ignore unchanged values
        if value == self.deletionFile:
            return
        # Validate file existence
        if value != None:
            if not os.path.isfile(value):
                raise FileNotFoundError(f"Deletion file '{value}' is not a file.")
            value = os.path.abspath(value) # store absolute path
            # Ignore again if value is unchanged after getting absolute path
            if value == self.deletionFile:
                return
        # Store updated value
        updateMsg = f"# Parameter cache: 'deletionFile' changed from '{self._deletionFile}' to '{value}'"
        self._deletionFile = value
        # Propagate change to deletion cache
        deletionCache = DeletionCache(self.workingDirectory)
        deletionCache.establish()
        deletionCache.deletionFile = value
        # Save after VCF cache is updated
        self.save() # try to keep the files in sync by allowing errors to be raised there first
        print(updateMsg)
    
    @property
    def bamSuffix(self):
        '''
        bamSuffix must be set BEFORE bamfiles, as the latter is validated using the former.
        '''
        return self._bamSuffix
    
    @bamSuffix.setter
    def bamSuffix(self, value):
        # Ignore unchanged values
        if value == self.bamSuffix:
            return
        # Validate that the suffix is not empty
        if value != None:
            if value == "":
                raise ValueError(f"--bamSuffix cannot be an empty string.")
        # Store and save
        updateMsg = f"# Parameter cache: 'bamSuffix' changed from '{self._bamSuffix}' to '{value}'"
        self._bamSuffix = value
        self.save()
        print(updateMsg)
    
    @property
    def bamFiles(self):
        return self._bamFiles
    
    @bamFiles.setter
    def bamFiles(self, value):
        # Ignore unchanged values
        if value == self.bamFiles:
            return
        # Locate and validate files
        foundBAMs = value
        if value != None:
            # Error out if bamSuffix is not set
            if self.bamSuffix == None:
                raise ValueError("--bamSuffix must be set for --bam files to be validated.")
            
            # Iterate through provided locations
            bamPrefixes = set()
            foundBAMs = []
            for location in value:
                location = os.path.abspath(location)
                
                # Handle an existing file
                if os.path.isfile(location):
                    if not location.endswith(self.bamSuffix):
                        raise ValueError(f"BAM file '{location}' does not end with the specified suffix '{self.bamSuffix}'")
                    else:
                        bamPrefix = os.path.basename(location).rsplit(self.bamSuffix, maxsplit=1)[0]
                        if bamPrefix in bamPrefixes:
                            raise ValueError(f"Duplicate BAM prefix found: '{bamPrefix}'")
                        
                        foundBAMs.append(location)
                # Handle an existing directory
                elif os.path.isdir(location):
                    foundAny = False
                    for f in os.listdir(location):
                        if f.endswith(self.bamSuffix):
                            bamPrefix = os.path.basename(location).rsplit(self.bamSuffix, maxsplit=1)[0]
                            if bamPrefix in bamPrefixes:
                                raise ValueError(f"Duplicate BAM prefix found: '{bamPrefix}'")
                            
                            foundBAMs.append(os.path.join(location, f))
                            foundAny = True
                    if not foundAny:
                        raise FileNotFoundError(f"No BAM files found in directory '{location}' ending with '{self.bamSuffix}'")
                # Error out if location does not exist
                else:
                    raise FileNotFoundError(f"Input BAM file or directory '{location}' not found!")
        # Ignore again if value is unchanged
        """The value from the cmdline may differ from the cached value if user input includes directories
        which get expanded out into full constituent file paths. This isn't exactly planned behaviour but
        it does allow for re-validation of the files in the cache which is nice."""
        if foundBAMs == self.bamFiles:
            return
        # Store and save
        updateMsg = f"# Parameter cache: 'bamFiles' changed from '{self._bamFiles}' to '{foundBAMs}'"
        self._bamFiles = foundBAMs
        self.save()
        print(updateMsg)
    
    @property
    def windowSize(self):
        return self._windowSize
    
    @windowSize.setter
    def windowSize(self, value):
        # Ignore unchanged values
        if value == self.windowSize:
            return
        # Validate value format and sensibility
        if value != None:
            if not type(value) == int:
                raise ValueError("--windowSize must be an integer.")
            if value < 1:
                raise ValueError("--windowSize must be a positive integer.")
        # Store and save
        updateMsg = f"# Parameter cache: 'windowSize' changed from '{self._windowSize}' to '{value}'"
        self._windowSize = value
        self.save()
        print(updateMsg)
    
    @property
    def qualFilter(self):
        return self._qualFilter
    
    @qualFilter.setter
    def qualFilter(self, value):
        # Ignore unchanged values
        if value == self.qualFilter:
            return
        # Validate value format and sensibility
        if value != None:
            try:
                value = float(value)
            except ValueError:
                raise ValueError("--qual must be a float or integer number.")
            if value < 0:
                raise ValueError("--qual must be at least 0.")
        # Store and save
        updateMsg = f"# Parameter cache: 'qualFilter' changed from '{self._qualFilter}' to '{value}'"
        self._qualFilter = value
        self.save()
        print(updateMsg)
    
    @property
    def missingFilter(self):
        return self._missingFilter
    
    @missingFilter.setter
    def missingFilter(self, value):
        # Ignore unchanged values
        if value == self.missingFilter:
            return
        # Validate value format and sensibility
        if value != None:
            try:
                value = float(value)
            except ValueError:
                raise ValueError("--missing must be a float or integer number.")
            if value < 0:
                raise ValueError("--missing must be at least 0.")
            if value > 1:
                raise ValueError("--missing must be at most 1.")
        # Store and save
        updateMsg = f"# Parameter cache: 'missingFilter' changed from '{self._missingFilter}' to '{value}'"
        self._missingFilter = value
        self.save()
        print(updateMsg)

class VcfCache:
    def __init__(self, workingDirectory):
        if not os.path.isdir(workingDirectory):
            raise FileNotFoundError(f"Directory '{workingDirectory}' is not a directory.")
        self.workingDirectory = os.path.abspath(workingDirectory)
        
        self._vcfFile = None
        self._variants = None
        self._samples = None
        self._contigs = None
        self._filteredVcfFile = None
        self._filteredVariants = None
        self._filteredSamples = None
        self._filteredContigs = None
    
    def establish(self):
        if not os.path.exists(self.cacheFile):
            self.initialise()
        else:
            self.load()
    
    def initialise(self):
        if os.path.exists(self.cacheFile):
            raise FileExistsError(f"Directory '{self.workingDirectory}' has already been initialised.")
        
        self.vcfFile = None
        self.filteredVcfFile = None
        self._variants = None # no setter for this and below
        self._samples = None
        self._contigs = None
        self._filteredVariants = None
        self._filteredSamples = None
        self._filteredContigs = None
        self.save()
    
    def load(self):
        if os.path.exists(self.cacheFile):
            with open(self.cacheFile, "r") as fileIn:
                data = json.load(fileIn)
            # Load data while circumventing the setters
            "Circumvent the setters to avoid re-parsing the VCF file(s)"
            self._vcfFile = data["vcfFile"]
            self._variants = data["variants"]
            self._samples = data["samples"]
            self._contigs = data["contigs"]
            self._filteredVcfFile = data["filteredVcfFile"]
            self._filteredVariants = data["filteredVariants"]
            self._filteredSamples = data["filteredSamples"]
            self._filteredContigs = data["filteredContigs"]
        else:
            raise FileNotFoundError(f"Working directory '{self.workingDirectory}' has not had its VCF cache initialised.")
    
    def save(self):
        data = {
            "vcfFile": self.vcfFile,
            "variants": self.variants,
            "samples": self.samples,
            "contigs": self.contigs,
            "filteredVcfFile": self.filteredVcfFile,
            "filteredVariants": self.filteredVariants,
            "filteredSamples": self.filteredSamples,
            "filteredContigs": self.filteredContigs
        }
        with open(self.cacheFile, "w") as fileOut:
            json.dump(data, fileOut)
    
    def _parse_vcf(self):
        if self._vcfFile == None or not os.path.isfile(self._vcfFile):
            raise FileNotFoundError(f"VCF file '{self._vcfFile}' is not a file.")
        self._variants, self._samples, self._contigs = parse_vcf_stats(self._vcfFile)
    
    def _parse_filtered_vcf(self):
        if self._filteredVcfFile == None or not os.path.isfile(self._filteredVcfFile):
            raise FileNotFoundError(f"Filtered VCF file '{self._filteredVcfFile}' is not a file.")
        self._filteredVariants, self._filteredSamples, self._filteredContigs = parse_vcf_stats(self._filteredVcfFile)
    
    @property
    def cacheFile(self):
        return os.path.join(self.workingDirectory, "vcf.json")
    
    @property
    def vcfFile(self):
        return self._vcfFile
    
    @vcfFile.setter
    def vcfFile(self, value):
        # Ignore unchanged values
        if value == self.vcfFile:
            return
        # Handle None values
        if value == None:
            print(f"# VCF cache: 'vcfFile' being blanked to '{value}'")
            self._vcfFile = None
            self._variants = None
            self._samples = None
            self._contigs = None
            self.save()
            return
        # Validate file existence
        if not os.path.isfile(value):
            raise FileNotFoundError(f"VCF file '{value}' is not a file.")
        value = os.path.abspath(value) # store absolute path
        # Ignore again if value is unchanged after getting absolute path
        if value == self.vcfFile:
            return
        # Store and parse
        updateMsg = f"# VCF cache: 'vcfFile' changed from '{self._vcfFile}' to '{value}'"
        self._vcfFile = value
        self._parse_vcf()
        self.save() # save after all changes are made
        print(updateMsg)
    
    @property
    def variants(self):
        return self._variants
    
    @property
    def samples(self):
        return self._samples
    
    @property
    def contigs(self):
        return self._contigs
    
    @property
    def filteredVcfFile(self):
        return self._filteredVcfFile
    
    @filteredVcfFile.setter
    def filteredVcfFile(self, value):
        # Ignore unchanged values
        if value == self.filteredVcfFile:
            return
        # Handle None values
        if value == None:
            print(f"# VCF cache: 'filteredVcfFile' being blanked to '{value}'")
            self._filteredVcfFile = None
            self._filteredVariants = None
            self._filteredSamples = None
            self._filteredContigs = None
            self.save()
            return
        # Validate file existence
        if not os.path.isfile(value):
            raise FileNotFoundError(f"Filtered VCF file '{value}' is not a file.")
        value = os.path.abspath(value) # store absolute path
        # Ignore again if value is unchanged after getting absolute path
        if value == self.filteredVcfFile:
            return
        # Store and parse
        updateMsg = f"# VCF cache: 'filteredVcfFile' changed from '{self._filteredVcfFile}' to '{value}'"
        self._filteredVcfFile = value
        self._parse_filtered_vcf()
        self.save() # save after all changes are made
        print(updateMsg)
    
    @property
    def filteredVariants(self):
        return self._filteredVariants
    
    @property
    def filteredSamples(self):
        return self._filteredSamples
    
    @property
    def filteredContigs(self):
        return self._filteredContigs

class DeletionCache:
    def __init__(self, workingDirectory):
        if not os.path.isdir(workingDirectory):
            raise FileNotFoundError(f"Directory '{workingDirectory}' is not a directory.")
        self.workingDirectory = os.path.abspath(workingDirectory)
        
        self._deletionFile = None
        self._bins = None
        self._deletionBins = None
        self._samples = None
        self._contigs = None
    
    def establish(self):
        if not os.path.exists(self.cacheFile):
            self.initialise()
        else:
            self.load()
    
    def initialise(self):
        if os.path.exists(self.cacheFile):
            raise FileExistsError(f"Directory '{self.workingDirectory}' has already been initialised.")
        
        self.deletionFile = None
        self._bins = None # no setter for this and below
        self._deletionBins = None
        self._samples = None
        self._contigs = None
        self.save()
    
    def load(self):
        if os.path.exists(self.cacheFile):
            with open(self.cacheFile, "r") as fileIn:
                data = json.load(fileIn)
            # Load data while circumventing the setters
            "Circumvent the setters to avoid re-parsing the VCF-like file"
            self._deletionFile = data["deletionFile"]
            self._bins = data["bins"]
            self._deletionBins = data["deletionBins"]
            self._samples = data["samples"]
            self._contigs = data["contigs"]
        else:
            raise FileNotFoundError(f"Working directory '{self.workingDirectory}' has not had its deletion cache initialised.")
    
    def save(self):
        data = {
            "deletionFile": self.deletionFile,
            "bins": self.bins,
            "deletionBins": self.deletionBins,
            "samples": self.samples,
            "contigs": self.contigs
        }
        with open(self.cacheFile, "w") as fileOut:
            json.dump(data, fileOut)
    
    def _parse_deletion_vcf(self):
        if self.deletionFile == None or not os.path.isfile(self.deletionFile):
            raise FileNotFoundError(f"Deletion file '{self.deletionFile}' is not a file.")
        self._bins, self._deletionBins, self._samples, self._contigs = parse_deletion_stats(self.deletionFile)
    
    @property
    def cacheFile(self):
        return os.path.join(self.workingDirectory, "deletion.json")
    
    @property
    def deletionFile(self):
        return self._deletionFile
    
    @deletionFile.setter
    def deletionFile(self, value):
        # Ignore unchanged values
        if value == self.deletionFile:
            return
        # Handle None values
        if value == None:
            print(f"# Deletion cache: 'deletionFile' being blanked to '{value}'")
            self._deletionFile = None
            self._bins = None
            self._deletionBins = None
            self._samples = None
            self._contigs = None
            self.save()
            return
        # Validate file existence
        if not os.path.isfile(value):
            raise FileNotFoundError(f"Deletion file '{value}' is not a file.")
        value = os.path.abspath(value) # store absolute path
        # Ignore again if value is unchanged after getting absolute path
        if value == self.deletionFile:
            return
        # Store and parse
        updateMsg = f"# Deletion cache: 'deletionFile' changed from '{self._deletionFile}' to '{value}'"
        self._deletionFile = value
        self._parse_deletion_vcf()
        self.save() # save after all changes are made
        print(updateMsg)
    
    @property
    def bins(self):
        return self._bins
    
    @property
    def deletionBins(self):
        return self._deletionBins
    
    @property
    def samples(self):
        return self._samples
    
    @property
    def contigs(self):
        return self._contigs

class MetadataCache:
    def __init__(self, workingDirectory):
        if not os.path.isdir(workingDirectory):
            raise FileNotFoundError(f"Directory '{workingDirectory}' is not a directory.")
        self.workingDirectory = os.path.abspath(workingDirectory)
        
        self._metadataFile = None
        self._bulk1 = None
        self._bulk2 = None
    
    def establish(self):
        if not os.path.exists(self.cacheFile):
            self.initialise()
        else:
            self.load()
    
    def initialise(self):
        if os.path.exists(self.cacheFile):
            raise FileExistsError(f"Directory '{self.workingDirectory}' has already been initialised.")
        
        self.metadataFile = None
        self._bulk1 = None # no setter for this and below
        self._bulk2 = None
        self.save()
    
    def load(self):
        if os.path.exists(self.cacheFile):
            with open(self.cacheFile, "r") as fileIn:
                data = json.load(fileIn)
            # Load data while circumventing the setters
            "Circumventing is only necessary here because of how None values are handled later"
            self._metadataFile = data["metadataFile"]
            self._bulk1 = data["bulk1"]
            self._bulk2 = data["bulk2"]
        else:
            raise FileNotFoundError(f"Working directory '{self.workingDirectory}' has not had its metadata cache initialised.")
    
    def save(self):
        data = {
            "metadataFile": self.metadataFile,
            "bulk1": self.bulk1,
            "bulk2": self.bulk2
        }
        with open(self.cacheFile, "w") as fileOut:
            json.dump(data, fileOut)
    
    def _parse_metadata(self):
        if self.metadataFile == None or not os.path.isfile(self.metadataFile):
            raise FileNotFoundError(f"Metadata file '{self.metadataFile}' is not a file.")
        metadataDict = parse_metadata(self.metadataFile)
        self._bulk1, self._bulk2 = metadataDict["bulk1"], metadataDict["bulk2"]
    
    @property
    def cacheFile(self):
        return os.path.join(self.workingDirectory, "metadata.json")
    
    @property
    def metadataFile(self):
        return self._metadataFile
    
    @metadataFile.setter
    def metadataFile(self, value):
        # Ignore unchanged values
        if value == self.metadataFile:
            return
        # Handle None values
        if value == None:
            print(f"# Metadata cache: 'metadataFile' being blanked to '{value}'")
            self._metadataFile = None
            self._bulk1 = None
            self._bulk2 = None
            self.save()
            return
        # Validate file existence
        if not os.path.isfile(value):
            raise FileNotFoundError(f"Metadata file '{value}' is not a file.")
        value = os.path.abspath(value) # store absolute path
        # Ignore again if value is unchanged after getting absolute path
        if value == self.metadataFile:
            return
        # Store and parse
        updateMsg = f"# Metadata cache: 'metadataFile' changed from '{self._metadataFile}' to '{value}'"
        self._metadataFile = value
        self._parse_metadata()
        self.save() # save after all changes are made
        print(updateMsg)
    
    @property
    def bulk1(self):
        return self._bulk1
    
    @property
    def bulk2(self):
        return self._bulk2
