import os, subprocess, re
from pycirclize.parser import Gff

from .locations import Locations
from .parameters import ParameterCache
from .parsing import parse_metadata
from .gff3 import GFF3Graph
from .plot import NUM_SAMPLE_LINES

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
    
    # Validate ploidy number (only for 'call' mode)
    if args.ploidyNum < 1:
        raise ValueError(f"--ploidy value '{args.ploidyNum}' must be >= 1!")
    
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
    if args.depthFile == None:
        raise FileNotFoundError("Working directory has not been initialised with a depth file!")
    else:
        args.depthFile = os.path.abspath(args.depthFile)
        if not os.path.isfile(args.depthFile):
            raise FileNotFoundError(f"Depth file '{args.depthFile}' was identified in " +
                                    "the parameters cache, but it doesn't exist or is not a file!")

def validate_s(args):
    # Validate numeric arguments
    if args.splsdaWindowSize < 1:
        raise ValueError(f"--windowSize value '{args.splsdaWindowSize}' must be >= 1!")
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
    
    # Format inputType and measurementType for compatibility between plot and report modes
    "Report mode only accepts a modified measurementType value, so we need to reformat it accordingly."
    if args.mode == "report":
        if args.measurementType == "ed-call":
            args.inputType = ["call"]
            args.measurementType = ["ed"]
        elif args.measurementType == "ed-depth":
            args.inputType = ["depth"]
            args.measurementType = ["ed"]
        elif args.measurementType == "splsda":
            args.inputType = []
            # Add 'call' input if the output files exist
            if os.path.isfile(locations.variantSplsdaSelectedFile) and os.path.isfile(locations.variantSplsdaSelectedFile + ".ok"):
                args.inputType.append("call")
            # Add 'depth' input if the output files exist
            if os.path.isfile(locations.depthSplsdaSelectedFile) and os.path.isfile(locations.depthSplsdaSelectedFile + ".ok"):
                args.inputType.append("depth")
            # If neither input type is available, raise an error
            if args.inputType == []:
                raise FileNotFoundError("No sPLS-DA selection files found for 'call' or 'depth'!")
            args.measurementType = ["splsda"]
    
    # Validate ED input type
    if "call" in args.inputType:
        if args.edType == "alleles":
            args.edFile = locations.allelesEdFile
            args.pickleFile = locations.allelesEdPickleFile # store the function
        elif args.edType == "inheritance":
            args.edFile = locations.inheritanceEdFile
            args.pickleFile = locations.inheritanceEdPickleFile
        elif args.edType == "genotypes":
            args.edFile = locations.genotypesEdFile
            args.pickleFile = locations.genotypesEdPickleFile
        else:
            raise ValueError(f"Invalid --ed '{args.edType}'! Must be one of 'alleles', 'inheritance', or 'genotypes'.")
    
    # Locate and validate input files
    if "call" in args.inputType:
        if "ed" in args.measurementType:
            if not os.path.isfile(args.edFile):
                raise FileNotFoundError(f"'call' ED file '{args.edFile}' does not exist!")
            if not os.path.isfile(args.edFile + ".ok"):
                raise FileNotFoundError(f"'call' ED file '{args.edFile}' does not have a '.ok' flag!")
        
        if "splsda" in args.measurementType:
            if not os.path.isfile(locations.variantSplsdaSelectedFile):
                raise FileNotFoundError(f"sPLS-DA file for 'call' selected features '{locations.variantSplsdaSelectedFile}' does not exist!")
            if not os.path.isfile(locations.variantSplsdaSelectedFile + ".ok"):
                raise FileNotFoundError(f"sPLS-DA file for 'call' selected features '{locations.variantSplsdaSelectedFile}' does not have a '.ok' flag!")
            
            if not os.path.isfile(locations.variantSplsdaBerFile):
                raise FileNotFoundError(f"sPLS-DA file for 'call' Balanced Error Rate '{locations.variantSplsdaBerFile}' does not exist!")
            if not os.path.isfile(locations.variantSplsdaBerFile + ".ok"):
                raise FileNotFoundError(f"sPLS-DA file for 'call' Balanced Error Rate '{locations.variantSplsdaBerFile}' does not have a '.ok' flag!")
    if "depth" in args.inputType:
        if "ed" in args.measurementType:
            if not os.path.isfile(locations.depthEdFile):
                raise FileNotFoundError(f"'depth' ED file '{locations.depthEdFile}' does not exist!")
            if not os.path.isfile(locations.depthEdFile + ".ok"):
                raise FileNotFoundError(f"'depth' ED file '{locations.depthEdFile}' does not have a '.ok' flag!")
        
        if "splsda" in args.measurementType:
            if not os.path.isfile(locations.depthSplsdaSelectedFile):
                raise FileNotFoundError(f"sPLS-DA file for 'depth' selected features '{locations.depthSplsdaSelectedFile}' does not exist!")
            if not os.path.isfile(locations.depthSplsdaSelectedFile + ".ok"):
                raise FileNotFoundError(f"sPLS-DA file for 'depth' selected features '{locations.depthSplsdaSelectedFile}' does not have a '.ok' flag!")
            
            if not os.path.isfile(locations.depthSplsdaBerFile):
                raise FileNotFoundError(f"sPLS-DA file for 'depth' Balanced Error Rate '{locations.depthSplsdaBerFile}' does not exist!")
            if not os.path.isfile(locations.depthSplsdaBerFile + ".ok"):
                raise FileNotFoundError(f"sPLS-DA file for 'depth' Balanced Error Rate '{locations.depthSplsdaBerFile}' does not have a '.ok' flag!")
    if "call" in args.inputType and "depth" in args.inputType:
        if "splsda" in args.measurementType:
            if not os.path.isfile(locations.integrativeSplsdaSelectedFile):
                print("# Note: '-i call depth' and '-m splsda' was set, but no integrated sPLS-DA result was found; if you have not run " +
                      "'psQTL_proc.py splsda -i call depth' yet, you may want to do that first. Otherwise, if running that command provided " +
                      "a message indicating that 'Integrative sPLS-DA is not possible or necessary', then you should ignore this message.")
    
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
            args.gff3Obj = GFF3Graph(args.annotationGFF3) # parsing now to raise errors early
            
            # Validate that GFF3 is not empty
            if args.gff3Obj.features == {}:
                raise ValueError(f"-a/--annotation file '{args.annotationGFF3}' is empty?")
            
            # Validate that GFF3 contains genes and mRNAs
            if not "gene" in args.gff3Obj.ftypes or not "mRNA" in args.gff3Obj.ftypes:
                raise ValueError(f"-a/--annotation file '{args.annotationGFF3}' does not contain 'gene' and/or 'mRNA' features; " +
                                 "psQTL requires these features to be present to make use of a GFF3 annotation file.")
            
            # Create NCLS index and perform QC
            args.gff3Obj.create_ncls_index("gene")
            args.gff3Obj.qc(typesToCheck=["gene", "mRNA"]) # prints out warnings if any issues found
    
    # Validate output file
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o output file '{args.outputFileName}' already exists!")
    args.outputFileName = os.path.abspath(args.outputFileName)
    
    if not os.path.isdir(os.path.dirname(args.outputFileName)):
        raise FileNotFoundError(f"-o parent directory '{os.path.dirname(args.outputFileName)}' does not exist!")
    
    return locations

def validate_regions(regions, mode, plotStyle, lengthsDict, argName="--region"):
    '''
    Returns:
        parsedRegions -- a list of lists with structure like:
                         [
                             [contigID, start, end, reverse],
                             ...
                         ]
    '''
    # Parse regions
    parsedRegions = []
    regionsRegex = re.compile(r"^([^:]+):(\d+)-(\d+)$")
    for region in regions:
        reMatch = regionsRegex.match(region)
        
        # Handle chr:start-end format
        if reMatch != None:
            contigID, start, end = reMatch.groups()
            start = int(start)
            end = int(end)
            
            # Validate contig ID
            if not contigID in lengthsDict:
                raise ValueError(f"{argName} contig ID '{contigID}' not found in the -f FASTA!")
            
            # Validate start position
            if start < 0:
                raise ValueError(f"{argName} start position '{start}' is < 0!")
            if start == end:
                raise ValueError(f"{argName} start position '{start}' is equal to end position '{end}'!")
            
            # Detect reverse orientation and swap start/end if necessary
            reverse = False
            if start > end:
                if mode == "plot":
                    # Prevent reverse orientation if plotStyle == "circos"
                    if plotStyle == "circos":
                        raise ValueError(f"{argName} '{contigID, start, end}' cannot be in reverse orientation " + 
                                        "with '-s circos'; only '-s horizontal' plot style can plot in reverse.")
                    # Otherwise, swap start and end
                    else:
                        start, end = end, start
                        reverse = True
                # If running 'report' mode, just swap start and end
                else:
                    start, end = end, start
                    reverse = True # has no effect in 'report' mode, but keeps logic consistent
            
            # Validate end position
            if end > lengthsDict[contigID]:
                raise ValueError(f"{argName} '{contigID, start, end}' end position is > contig length '{lengthsDict[contigID]}'!")
            
            # Store region
            parsedRegions.append([contigID, start, end, reverse])
        
        # Handle invalid format
        elif ":" in region:
            raise ValueError(f"Invalid region input '{region}'; you included a ':' but did " + 
                             "not format the region as 'chr:start-end'!")
        # Handle chr format
        else:
            if not region in lengthsDict:
                raise ValueError(f"{argName} contig ID '{region}' not found in the -f FASTA!")
            parsedRegions.append([region, 0, lengthsDict[region], False])
    
    # Handle empty regions
    if parsedRegions == []:
        parsedRegions = [[contigID, 0, lengthsDict[contigID], False] for contigID in lengthsDict]
    
    return parsedRegions

def validate_depth_files(depthDir, metadataDict, windowSize):
    '''
    Parameters:
        depthDir -- a string indicating the path to the directory containing
                    depth files
        metadataDict -- a dictionary with structure like:
                        {
                            "group1": [sampleID, ...],
                            "group2": [sampleID, ...]
                        }
        windowSize -- the window size used when generating the depth files
    Returns:
        depthFileDict --  a dictionary with structure like:
                          {
                              "group1": [[sampleID, depthFile], ...],
                              "group2": [[sampleID, depthFile], ...]
                          }
    '''
    notFound = []
    depthFileDict = {"group1": [], "group2": []}
    for group, sampleList in metadataDict.items():
        for sample in sampleList:
            depthFile = os.path.join(depthDir, f"{sample}.binned.{windowSize}.tsv")
            if not os.path.isfile(depthFile):
                notFound.append(sample)
            else:
                depthFileDict[group].append([sample, depthFile])
    if notFound != []:
        raise FileNotFoundError(f"Could not find depth files with bin size of {windowSize} " +
                                f"for samples: {', '.join(notFound)}")
    return depthFileDict

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
    if args.axisSpace < 0:
        raise ValueError(f"--space value '{args.axisSpace}' must be >= 0!")
    if args.axisSpace >= 360:
        raise ValueError(f"--space value '{args.axisSpace}' must be < 360!")
    
    # Alert user to Circos dimensions logic
    if args.plotStyle == "circos":
        if (args.width != None and args.height != None):
            if args.width != 8 or args.height != 8:
                print("# Note: psQTL is not yet equipped to handle non-default width or height while maintaining " +
                      "clarity for Circos plots, so width and height will be set to 8.")
                args.width = 8
                args.height = 8
            if args.width != args.height:
                print(f"# Note: Circos plotting enforces square dimensions, so width and height " +
                    f"will be set to the smallest value provided i.e., '{min(args.width, args.height)}'.")
    
    # Validate plot types
    if len(set(args.plotTypes)) != len(args.plotTypes):
        raise ValueError(f"-p must not contain duplicate plot types!")
    if "genes" in args.plotTypes and args.annotationGFF3 == None:
        raise ValueError(f"Cannot plot gene locations without providing an --annotation GFF3 file!")
    
    # Validate samples for coverage plot
    if "coverage" in args.plotTypes:
        for sampleID in args.coverageSamples:
            if not sampleID in args.metadataDict["group1"] + args.metadataDict["group2"]:
                raise ValueError(f"Sample '{sampleID}' specified in --coverageSamples not found in metadata!")
        if len(args.coverageSamples) > NUM_SAMPLE_LINES:
            raise ValueError(f"Cannot plot more than {NUM_SAMPLE_LINES} samples using --coverageSamples (for clarity reasons)!")
    
    # Validate argument logic
    if "coverage" in args.plotTypes and not "depth" in args.inputType:
        raise ValueError(f"Cannot specify '-p coverage' without also specifying '-i depth'!")
    
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
