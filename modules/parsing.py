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
        metadataDict -- a dictionary with structure like:
                        {
                            "bulk1": [ "sample1", "sample2", ... ],
                            "bulk2": [ "sample3", "sample4", ... ]
                        }
    '''
    ACCEPTED_BULK1 = ['bulk1', '1', 'bulk 1', 'b1']
    ACCEPTED_BULK2 = ['bulk2', '2', 'bulk 2', 'b2']
    
    metadataDict = {}
    foundSamples = set()
    with read_gz_file(metadataFile) as fileIn:
        for line in fileIn:
            l = line.rstrip("\r\n ")
            # Skip blank lines
            if l == "":
                continue
            
            # Skip comment lines
            if line.startswith("#"):
                continue
            
            # Split line based on delimiter
            if "\t" in line:
                sl = l.split("\t")
            elif "," in line:
                sl = l.split(",")
            else:
                raise ValueError(f"Metadata file is not tab or comma-delimited; offending line is '{l}'")
            
            # Strip any junk and Excel quotes
            sl = [s.strip("\"' ") for s in sl]
            
            # Parse out relevant information
            try:
                sample, pop = sl
                sample, pop = sample.strip(), pop.strip().lower() # make lowercase for easier comparison
            except ValueError:
                raise ValueError(f"Metadata file does not have two columns; offending line is '{l}'")
            
            # Validate that the population is one of the expected values
            if not pop in ACCEPTED_BULK1 + ACCEPTED_BULK2:
                raise ValueError(f"'{pop}' is not in the expected format for denoting bulks; offending line is '{l}'")
            
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

def parse_vcf_stats(vcfFile, refAllele="0"):
    '''
    Parses a VCF file to determine the number of variants, the sample IDs, and
    contig IDs.
    
    Parameters:
        vcfFile -- a string indicating the path to a VCF file
        refAllele -- (OPTIONAL) a string indicating the reference allele to use for
                     determining whether a position has a variant; default is "0"
    Returns:
        positions -- an integer indicating the number of positions listed in the VCF file
        variants -- an integer indicating the number of positions with a variant in the VCF file
        samples -- a list of strings indicating the sample IDs in the VCF file
        contigs -- a list of strings indicating the contig IDs in the VCF file
    '''
    positions = 0
    variants = 0
    samples = None
    contigs = {} # using as an ordered set
    with read_gz_file(vcfFile) as fileIn:
        for line in fileIn:
            l = line.strip('\r\n\t "') # remove quotations to help with files opened by Excel
            
            # Skip blank lines
            if l == "":
                continue
            
            # Handle header lines
            if l.startswith("#CHROM"):
                try:
                    samples = l.split("\t")[9:]
                except:
                    raise ValueError(f"#CHROM header line is malformed; offending line is '{l}'")
            elif l.startswith("#"):
                continue
            
            # Handle variant lines
            else:
                # Split line based on expected delimiter
                try:
                    sl = l.split("\t")
                    contig = sl[0]
                except:
                    raise ValueError(f"VCF file lacks tab-delimited formatting; offending line is '{l}'")
                
                # Validate line length
                if len(sl) < 11: # 9 fixed columns + 2 genotype columns
                    raise ValueError(f"VCF file has too few columns; offending line is '{l}'")
                formatField = sl[8]
                
                # Determine the GT field index
                if not ":" in formatField:
                    gtIndex = 0
                else:
                    gtIndex = formatField.split(":").index("GT")
                
                # Store contig ID if not already present
                if contig not in contigs:
                    contigs[contig] = None
                
                # Tally positions and variants
                positions += 1
                variants += 1 if any([
                        allele != refAllele
                        for sampleResult in sl[9:]
                        for allele in sampleResult.split(":")[gtIndex].replace("|", "/").split("/")
                    ]) else 0
    
    # Raise errors for faulty VCF files
    if positions == 0 or contigs == {}:
        raise ValueError("VCF file is empty or has no variants!")
    if samples is None:
        raise ValueError("VCF file has no #CHROM header line!")
    
    # Format and return results
    contigs = list(contigs.keys()) # convert to list
    return positions, variants, samples, contigs

def vcf_header_to_metadata_validation(vcfSamples, metadataDict, strict=True, quiet=True):
    '''
    Validates that the VCF file and metadata file are compatible with each other,
    printing warnings and producing errors as necessary.
    
    Parameters:
        vcfSamples -- a list or set of strings indicating the sample IDs in the VCF file.
        metadataDict -- a dictionary with structure like:
                        {
                            "bulk1": set([ "sample1", "sample2", ... ]),
                            "bulk2": set([ "sample3", "sample4", ... ])
                        }
        strict -- (OPTIONAL) a boolean indicating whether to error out if the files
                  show any discrepancies. Default is True.
        quiet -- (OPTIONAL) a boolean indicating whether to hide the final sample
                 lists to the console. Default is True.
    '''
    vcfSamples = set(vcfSamples) # convert to set if not already one
    
    # Extract metadata sample IDs from dict
    metadataSamples = set([
        sample
        for samples in metadataDict.values()
        for sample in samples
    ])
    
    # Check for sample ID discrepancies between the two files
    if vcfSamples != metadataSamples:
        # Find differences between the two sample ID sets
        vcfDiff = vcfSamples.difference(metadataSamples)
        metadataDiff = metadataSamples.difference(vcfSamples)
        
        # Error out if files are incompatible
        if len(metadataDiff) == len(metadataSamples):
            raise ValueError("Metadata file has no samples in common with VCF file; " +
                             f"Metadata samples: {metadataSamples}\nVCF samples: {vcfSamples}")
        
        # Error out or warn if some samples are missing (depending on 'strict' parameter)
        if len(vcfDiff) > 0:
            if strict:
                raise ValueError(f"VCF file has samples not present in metadata file: {', '.join(vcfDiff)}")
            else:
                print("# WARNING: In your VCF, the following samples exist which are " +
                      "absent from the metadata: ", ", ".join(vcfDiff))
                print("# These samples will not be considered when generating the resulting file.")
        if len(metadataDiff) > 0:
            if strict:
                raise ValueError(f"Metadata file has samples not present in VCF file: {', '.join(metadataDiff)}")
            else:
                print("# WARNING: In your metadata, the following samples exist which are " + 
                      "absent from the VCF: ", ", ".join(metadataDiff))
                print("# These samples will not be considered when generating the resulting file.")
    
    # Error out if we don't have samples from both bulks
    b1Samples = [ sample for sample in vcfSamples if sample in metadataDict["bulk1"] ]
    if len(b1Samples) == 0:
        raise ValueError("No samples from bulk 1 are present in the VCF file.")
    
    b2Samples = [ sample for sample in vcfSamples if sample in metadataDict["bulk2"] ] 
    if len(b2Samples) == 0:
        raise ValueError("No samples from bulk 2 are present in the VCF file.")
    
    # Notify user of samples that will be used
    if not quiet:
        print(f"# Samples used as part of bulk 1 (n={len(b1Samples)}) include: " + ", ".join(b1Samples))
        print(f"# Samples used as part of bulk 2 (n={len(b2Samples)}) include: " + ", ".join(b2Samples))

def parse_vcf_genotypes(formatField, sampleFields, samples):
    '''
    Parameters:
        formatField -- a string from a split VCF line; if the split line is 'sl',
                       this should be 'sl[8]'
        sampleFields -- a list of string from a VCF line; if the split line is 'sl',
                        this should be 'sl[9:]'
        samples -- a list of strings indicating the sample IDs as obtained from the 
                   header line of the VCF file
    Returns:
        posGenotypeDict -- a dictionary with sample IDs as keys and genotype calls as
                           lists of integers; for example, '0/0' in the VCF file would
                           be stored as [0, 0] in the dictionary.
    '''
    # Determine which field position we're extracting to get our GT value
    if not ":" in formatField:
        gtIndex = 0
    else:
        gtIndex = formatField.split(":").index("GT")
    
    # Format a dictionary to store sample genotypes for this position
    posGenotypeDict = {}
    ongoingCount = 0 # This gives us the index for our samples header list 
    for sampleResult in sampleFields:
        # Grab our genotype
        genotype = sampleResult.split(":")[gtIndex]
        
        # Edit genotype to have a consistently predictable separator
        "Phasing information is irrelevant for psQTL analysis"
        genotype = genotype.replace("|", "/")
        
        # Skip uncalled genotypes
        if "." in genotype:
            ongoingCount += 1
            continue
        
        # Parse and store genotype
        sampleID = samples[ongoingCount]
        posGenotypeDict[sampleID] = list(map(int, genotype.split("/")))
        
        ongoingCount += 1
    return posGenotypeDict
