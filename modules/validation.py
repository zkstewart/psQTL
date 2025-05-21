import os, subprocess

from .locations import Locations
from .parameters import ParameterCache

def validate_prep_args(args):
    # Validate working directory
    args.workingDirectory = os.path.abspath(args.workingDirectory)
    if not os.path.exists(args.workingDirectory):
        # If we're initialising, we can create the directory
        if args.mode in ["initialise", "init"]:
            if not os.path.exists(os.path.dirname(args.workingDirectory)):
                raise FileNotFoundError(f"Parent directory for proposed -w location " + 
                                        f"'{os.path.dirname(args.workingDirectory)}' does not exist!")
            os.mkdir(args.workingDirectory)
            print(f"# Created working directory '{args.workingDirectory}'")
        # Otherwise, we need to error out
        else:
            raise FileNotFoundError(f"-d working directory '{args.workingDirectory}' does not exist!")
    
    # Validate cache existence
    if not args.mode in ["initialise", "init"]:
        paramsCache = ParameterCache(args.workingDirectory)
        paramsCache.load() # raises FileNotFoundError if cache does not exist
    
    # Establish locations object
    locations = Locations(args.workingDirectory)
    return locations

def validate_uncached(args):
    '''
    This function just needs to validate any arguments that are not found within the
    parameter cache.
    '''
    # Validate threads
    if args.threads < 1:
        raise ValueError("Number of threads must be at least 1!")
    
    # Validate genome FASTA file
    args.genomeFasta = os.path.abspath(args.genomeFasta)
    if not os.path.isfile(args.genomeFasta):
        raise FileNotFoundError(f"Genome FASTA file '{args.genomeFasta}' does not exist!")

def validate_proc_args(args):
    locations = Locations(args.workingDirectory) # performs validation implicitly
    
    # Validate cache existence & merge into args
    paramsCache = ParameterCache(args.workingDirectory)
    paramsCache.merge(args) # raises FileNotFoundError if cache does not exist
    
    # Validate metadata file
    if args.metadataFile == None:
        raise FileNotFoundError("Working directory has not been initialised with a metadata file!")
    elif not os.path.isfile(args.metadataFile):
            raise FileNotFoundError(f"Metadata file '{args.metadataFile}' was identified in " +
                                    "the parameters cache, but it doesn't exist or is not a file!")
    
    return locations

def validate_c(args):
    '''
    Params cache should have been merged into args before calling this function.
    '''
    # Choose which VCF file to use
    args.vcfFile = args.filteredVcfFile if args.filteredVcfFile != None else args.vcfFile
    
    # Validate VCF file
    if args.vcfFile == None:
        raise FileNotFoundError("Working directory has not been initialised with a VCF file!")
    else:
        args.vcfFile = os.path.abspath(args.vcfFile)
        if not os.path.isfile(args.vcfFile):
            raise FileNotFoundError(f"VCF file '{args.vcfFile}' was identified in " +
                                    "the parameters cache, but it doesn't exist or is not a file!")

def validate_d(args):
    '''
    Params cache should have been merged into args before calling this function.
    '''
    if args.deletionFile == None:
        raise FileNotFoundError("Working directory has not been initialised with a deletion file!")
    else:
        args.deletionFile = os.path.abspath(args.deletionFile)
        if not os.path.isfile(args.deletionFile):
            raise FileNotFoundError(f"Deletion file '{args.deletionFile}' was identified in " +
                                    "the parameters cache, but it doesn't exist or is not a file!")

def validate_s(args):
    # Validate numeric arguments
    if args.windowSize < 1:
        raise ValueError(f"--windowSize value '{args.windowSize}' must be >= 1!")
    if args.threads < 1:
        raise ValueError(f"--threads value '{args.threads}' must be >= 1!")
    if args.numRepeats < 1:
        raise ValueError(f"--nrepeat value '{args.numRepeats}' must be >= 1!")
    if args.maxIterations < 1:
        raise ValueError(f"--maxiters value '{args.maxIterations}' must be >= 1!")
    
    # Validate ratio arguments
    if args.mafFilter < 0 or args.mafFilter > 0.5:
        raise ValueError(f"--maf value '{args.mafFilter}' must be between 0 and 0.5!")
    if args.berFilter < 0 or args.berFilter > 0.5:
        raise ValueError(f"--ber value '{args.berFilter}' must be between 0 and 0.5!")

def validate_post_args(args):
    locations = Locations(args.workingDirectory) # performs validation implicitly
    
    # Validate cache existence & merge into args
    paramsCache = ParameterCache(args.workingDirectory)
    paramsCache.merge(args) # raises FileNotFoundError if cache does not exist
    
    # Validate metadata file
    if args.metadataFile == None:
        raise FileNotFoundError(f"Metadata file has not been initialised with psQTL_prep.py!")
    elif not os.path.isfile(args.metadataFile):
        raise FileNotFoundError(f"Metadata file '{args.metadataFile}' was initialised previously " + 
                                "but no longer exists!")
    else:
        args.metadataDict = parse_metadata(args.metadataFile)
    
    # Validate input file
    if args.inputType == "call":
        args.inputFile = os.path.join(args.workingDirectory, "psQTL_variants.ed.tsv.gz")
    else:
        args.inputFile = os.path.join(args.workingDirectory, "psQTL_depth.ed.tsv.gz")
    
    if not os.path.isfile(args.inputFile):
        raise FileNotFoundError(f"Euclidean distance file '{args.inputFile}' does not exist!")
    elif not os.path.isfile(args.inputFile + ".ok"):
        raise FileNotFoundError(f"Euclidean distance file '{args.inputFile}' does not have a '.ok' flag!")
    
    # Validate genome FASTA file
    if not os.path.isfile(args.genomeFasta):
        raise FileNotFoundError(f"-f '{args.genomeFasta}' is not a file!")
    
    # Validate numeric arguments
    if args.power < 1:
        raise ValueError(f"--power value '{args.power}' must be >= 1!")
    if args.missingFilter < 0 or args.missingFilter > 1:
        raise ValueError(f"--missing value '{args.missingFilter}' must be between 0 and 1!")
    
    # Validate annotation GFF3 file
    if args.annotationGFF3 != None:
        if not os.path.isfile(args.annotationGFF3):
            raise FileNotFoundError(f"-a/--annotation file '{args.annotationGFF3}' is not a file!")
        else:
            args.gff3Obj = GFF3(args.annotationGFF3) # parsing now to raise errors early
            args.gff3Obj.create_ncls_index("gene")
    
    # Validate output file
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o output file '{args.outputFileName}' already exists!")
    args.outputFileName = os.path.abspath(args.outputFileName)
    
    if not os.path.isdir(os.path.dirname(args.outputFileName)):
        raise FileNotFoundError(f"-o parent directory '{os.path.dirname(args.outputFileName)}' does not exist!")
    
    return locations

def validate_regions(args, lengthsDict):
    # Validate regions
    regions = []
    regionsRegex = re.compile(r"^([^:]+):(\d+)-(\d+)$")
    for region in args.regions:
        reMatch = regionsRegex.match(region)
        
        # Handle chr:start-end format
        if reMatch != None:
            contigID, start, end = reMatch.groups()
            start = int(start)
            end = int(end)
            
            # Validate contig ID
            if not contigID in lengthsDict:
                raise ValueError(f"--region contig ID '{contigID}' not found in the -f FASTA!")
            # Validate start and end positions
            if start < 0:
                raise ValueError(f"--region start position '{start}' is < 0!")
            reverse = False
            if start >= end:
                start, end = end, start
                reverse = True
            # Store region
            regions.append([contigID, start, end, reverse])
        # Handle invalid format
        elif ":" in region:
            raise ValueError(f"Invalid region input '{region}'; you included a ':' but did " + 
                             "not format the region as 'chr:start-end'!")
        # Handle chr format
        else:
            if not region in lengthsDict:
                raise ValueError(f"--region contig ID '{region}' not found in the -f FASTA!")
            regions.append([region, 0, lengthsDict[region], False])
    
    # Handle empty regions
    if regions == []:
        regions = [[contigID, 0, lengthsDict[contigID], False] for contigID in lengthsDict]
    
    args.regions = regions

def validate_p(args):
    # Validate numeric arguments
    if args.wmaSize < 1:
        raise ValueError(f"--wma value '{args.wmaSize}' must be >= 1!")
    if args.width != None:
        if args.width < 1:
            raise ValueError(f"--width value '{args.width}' must be >= 1!")
    if args.height != None:
        if args.height < 1:
            raise ValueError(f"--height value '{args.height}' must be >= 1!")
    if args.binSize < 2:
        raise ValueError(f"--bin value '{args.binSize}' must be >= 2!")
    if args.binThreshold < 0:
        raise ValueError(f"--threshold value '{args.binThreshold}' must be >= 0!")
    
    # Validate plot types
    if len(set(args.plotTypes)) != len(args.plotTypes):
        raise ValueError(f"-p must not contain duplicate plot types!")
    if "genes" in args.plotTypes and args.annotationGFF3 == None:
        raise ValueError(f"Cannot plot gene locations without providing an --annotation GFF3 file!")
    if "histogram" in args.plotTypes and args.inputType == "depth":
        raise ValueError(f"Cannot plot histogram for -i depth data!")
    
    # Validate samples for coverage plot
    if "coverage" in args.plotTypes:
        for sampleID in args.sampleCoverage:
            if not sampleID in args.metadataDict["bulk1"] + args.metadataDict["bulk2"]:
                raise ValueError(f"Sample '{sampleID}' specified in --sampleCoverage not found in metadata!")
        if len(args.sampleCoverage) > NUM_SAMPLE_LINES:
            raise ValueError(f"Cannot plot more than {NUM_SAMPLE_LINES} samples using --sampleCoverage for clarity")
    
    # Validate output file suffix
    if not (args.outputFileName.endswith(".pdf") or args.outputFileName.endswith(".png") or args.outputFileName.endswith(".svg")):
        raise ValueError(f"-o output file '{args.outputFileName}' must end with '.pdf', '.png', or '.svg'!")

def validate_r(args):
    # Validate numeric arguments
    if args.radiusSize < 0:
        raise ValueError(f"--radius value '{args.radiusSize}' must be >= 0!")
    
    # Validate output file suffix
    if not args.outputFileName.endswith(".tsv") and not args.outputFileName.endswith(".csv"):
        raise ValueError(f"-o output file '{args.outputFileName}' must end with '.tsv' or '.csv'!")
