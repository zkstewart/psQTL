import os, gzip, codecs
from contextlib import contextmanager

def get_codec(fileName):
    try:
        f = codecs.open(fileName, encoding='utf-8', errors='strict')
        for line in f:
            break
        return "utf-8"
    except:
        try:
            f = codecs.open(fileName, encoding='utf-16', errors='strict')
            for line in f:
                break
            return "utf-16"
        except UnicodeDecodeError:
            print(f"Can't tell what codec '{fileName}' is!!")

@contextmanager
def read_gz_file(filename):
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            yield f
    else:
        with open(filename, "r", encoding=get_codec(filename)) as f:
            yield f

def parse_metadata(metadataFile):
    '''
    Parameters:
        metadataFile -- a string indicating the path to a metadata file
    Returns:
        metadataDict -- a dictionary mapping population names to one of two bulks
    '''
    ACCEPTED_BULK1 = ['bulk1', '1', 'bulk 1', 'b1']
    ACCEPTED_BULK2 = ['bulk2', '2', 'bulk 2', 'b2']
    
    metadataDict = {}
    foundSamples = set()
    with read_gz_file(metadataFile) as fileIn:
        for line in fileIn:
            # Split line based on delimiter
            if "\t" in line:
                sl = line.rstrip("\r\n ").split("\t")
            elif "," in line:
                sl = line.rstrip("\r\n ").split(",")
            else:
                raise ValueError(f"Metadata file is not tab or comma-delimited; offending line is '{line}'")
            
            # Skip comment lines
            if line.startswith("#"):
                continue
            
            # Skip blank lines
            if sl == []:
                continue
            
            # Parse out relevant information
            try:
                sample, pop = sl
                sample, pop = sample.strip(), pop.strip().lower() # make lowercase for easier comparison
            except ValueError:
                raise ValueError(f"Metadata file does not have two columns; offending line is '{line}'")
            
            # Validate that the population is one of the expected values
            if not pop in ACCEPTED_BULK1 + ACCEPTED_BULK2:
                raise ValueError(f"'{pop}' is not in the expected format for denoting bulks; offending line is '{line}'")
            
            # Unify the population names
            if pop in ACCEPTED_BULK1:
                pop = 'bulk1'
            elif pop in ACCEPTED_BULK2:
                pop = 'bulk2'
            
            # Validate that the sample is non-redundant
            if sample in foundSamples:
                raise ValueError(f"Sample '{sample}' is listed more than once in the metadata file")
            foundSamples.add(sample)
            
            # Store the sample in the appropriate population
            metadataDict.setdefault(pop, set())
            metadataDict[pop].add(sample)
    
    # Make sure that the metadata file has the expected number of populations
    if len(metadataDict) == 0:
        raise ValueError("Metadata file is empty; please provide at least two samples from different bulks")
    if len(metadataDict) == 1:
        foundBulk = list(metadataDict.keys())[0]
        raise ValueError(f"Metadata file should have 2 populations ('bulk1' and 'bulk2'); I only found '{foundBulk}'")
    
    # Reformat sets and lists then return
    for pop in metadataDict:
        metadataDict[pop] = list(metadataDict[pop])
    return metadataDict

def parse_samtools_depth_tsv(depthFile):
    '''
    Parameters:
        depthFile -- a string indicating the path to a samtools depth file
    Yields:
        contigID -- a string indicating the contig ID
        pos -- an integer indicating the position on the contig
        depth -- an integer indicating the depth at this position
    '''
    with read_gz_file(depthFile) as fileIn:
        for line in fileIn:
            # Parse out relevant details from this line
            try:
                contigID, pos, depth = line.rstrip("\r\n ").split("\t")
                pos, depth = int(pos), int(depth)
            except ValueError:
                raise(ValueError(f"samtools depth file is improperly formatted; offending line is '{line}'"))
            yield contigID, pos, depth

def parse_binned_tsv(binFile):
    '''
    Receives a file created by depth binning of the 'samtools depth' output, and
    returns a dictionary mapping contig IDs to depth bins.
    
    Parameters:
        binFile -- a strings indicating the path to a depth bin file
    Returns:
        histoDict -- a dictionary mapping contig IDs to depth bins
    '''
    histoDict = {}
    with read_gz_file(binFile) as fileIn:
        for line in fileIn:
            # Parse out relevant details from this line
            try:
                contigID, pos, depth = line.rstrip("\r\n ").split("\t")
                pos, depth = int(pos), int(depth)
            except ValueError:
                raise(ValueError(f"Depth bin file is improperly formatted; offending line is '{line}'"))
            
            # Store the depth for this contig
            histoDict.setdefault(contigID, {})
            histoDict[contigID][pos] = depth
    return histoDict

class SimpleGenotypeIterator:
    '''
    A class to iterate through a VCF file and yield the genotype calls.
    The first value yielded is the sample IDs, with subsequent yields
    providing the following parameters:
        - chrom -- a string indicating the chromosome
        - pos -- an int indicating the position
        - ref -- a string indicating the reference allele
        - alt -- a list of strings indicating the alternate allele(s)
        - posGenotypeDict -- a dictionary with sample IDs as keys and
                             genotype calls as lists of integers (0 for ref, 1 for alt,
                             and so on...)
    '''
    def __init__(self, file_location):
        assert os.path.exists(file_location), \
            f"SimpleGenotypeIterator can't find file existing at '{file_location}'"
        
        self.fileLocation = file_location
    
    def parse(self):
        with read_gz_file(self.fileLocation) as fileIn:
            for line in fileIn:
                sl = line.rstrip("\r\n").replace('"', '').split("\t") # remove quotations to help with files opened by Excel
                
                # Handle header lines
                if line.startswith("#CHROM"):
                    samples = sl[9:] # This gives us the ordered sample IDs
                    yield samples
                    continue
                if line.startswith("#"):
                    continue
                
                # Extract relevant details
                chrom = sl[0]
                pos = int(sl[1])
                ref = sl[3]
                alt = sl[4].split(",")
                
                # Determine which field position we're extracting to get our GT value
                fieldsDescription = sl[8]
                if ":" not in fieldsDescription:
                    gtIndex = 0
                else:
                    gtIndex = fieldsDescription.split(":").index("GT")
                
                # Format a dictionary to store sample genotypes for this position
                posGenotypeDict = {}
                ongoingCount = 0 # This gives us the index for our samples header list 
                for sampleResult in sl[9:]: # This gives us the results for each sample as per fieldsDescription
                    # Grab our genotype
                    if gtIndex != -1:
                        genotype = sampleResult.split(":")[gtIndex]
                    else:
                        genotype = sampleResult
                    
                    # Edit genotype to have a consistently predictable separator
                    "We don't care if the VCF is phased or not for this function"
                    genotype = genotype.replace("/", "|")
                    
                    # Skip uncalled genotypes
                    if "." in genotype:
                        ongoingCount += 1
                        continue
                    
                    # Parse and store genotype
                    samplePopulation = samples[ongoingCount]
                    posGenotypeDict[samplePopulation] = list(map(int, genotype.split("|")))
                    
                    ongoingCount += 1
                
                # Yield result
                yield chrom, pos, ref, alt, posGenotypeDict
    
    def __iter__(self):
        for x in self.parse():
            yield x
