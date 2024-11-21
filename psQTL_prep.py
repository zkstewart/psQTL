#! python3
# psQTL_prep.py
# Represents step 1 in the psQTL pipeline, where the user can initialise a working directory,
# view the metadata of an analysis directory, call variants, and calculate depth files for samples.
# Called variants and depth files can be used as input to the psQTL_proc.py script.

import os, argparse, sys
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.parameters import ParameterCache, VcfCache, DeletionCache, MetadataCache
from modules.samtools_handling import validate_samtools_exists, run_samtools_depth, run_samtools_faidx, \
                                      bin_samtools_depth
from modules.bcftools_handling import validate_bcftools_exists, validate_bgzip_exists, validate_vt_exists, \
                                      run_bcftools_call, run_bcftools_index, run_normalisation, \
                                      run_bcftools_concat, run_bgzip
from modules.depth import call_deletions_from_depth
from modules.parsing import parse_metadata
from modules.filter import filter_vcf

def validate_args(args):
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

def main():
    usage = """%(prog)s manages several steps in preparation for the psQTL pipeline.
    A working directory must be 'initialise'd before downstream processing and post-
    processing takes place. An analysis directory can be 'view'ed to see what samples
    are present. Otherwise, if you have not done so already, you may optionally
    use this script to 'call' variants or calculate the 'depth' of coverage for your
    samples before proceeding with the psQTL pipeline.
    """
    # Establish main parser
    p = argparse.ArgumentParser()
    
    # Set arguments shared by subparsers
    p.add_argument("-d", dest="workingDirectory",
                    required=True,
                    help="Specify the location where the analysis is being performed")
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser(description=usage)
    subparsers = subParentParser.add_subparsers(dest="mode",
                                                required=True)
    
    iparser = subparsers.add_parser("initialise",
                                    aliases=["init"],
                                    parents=[p],
                                    add_help=False,
                                    help="Initialise a working directory for the psQTL pipeline")
    iparser.set_defaults(func=imain)
    
    dparser = subparsers.add_parser("depth",
                                    parents=[p],
                                    add_help=False,
                                    help="Generate samtools depth files per sample")
    dparser.set_defaults(func=dmain)
    
    cparser = subparsers.add_parser("call",
                                    parents=[p],
                                    add_help=False,
                                    help="Call variants per sample")
    cparser.set_defaults(func=cmain)
    
    vparser = subparsers.add_parser("view",
                                    parents=[p],
                                    add_help=False,
                                    help="View the metadata of an analysis directory")
    vparser.set_defaults(func=vmain)
    
    # Init-subparser arguments
    iparser.add_argument("--meta", dest="metadataFile",
                         required=False,
                         help="""Optionally, specify the location of a metadata TSV file containing two
                         columns indicating 1) sample ID and 2) the bulk it belongs to""")
    iparser.add_argument("--vcf", dest="vcfFile",
                         required=False,
                         help="""Optionally, specify a VCF file containg per-sample variant
                         calls that you have already produced""")
    iparser.add_argument("--fvcf", dest="filteredVcfFile",
                         required=False,
                         help="""Optionally, specify a filtered VCF file containg per-sample variant
                         calls that you have already produced""")
    iparser.add_argument("--deletion", dest="deletionFile",
                         required=False,
                         help="""Optionally, specify a deletion VCF-like file containing per-sample
                         deletion calls that you have already produced""")
    iparser.add_argument("--bam", dest="bamFiles",
                         required=False,
                         nargs="+",
                         help="""Optionally, specify one or more locations of BAM files and/or
                         directories containing BAM files for variant calling and/or depth
                         calculations""")
    iparser.add_argument("--bamSuffix", dest="bamSuffix",
                         required=False,
                         help="""Optionally, specify the suffix used to denote BAM files;
                         relevant if directories are provided to --bam""")
    iparser.add_argument("--windowSize", dest="windowSize",
                         type=int,
                         required=False,
                         help="""Optionally, specify the window size that reads will be
                         binned into for deletion calling""")
    iparser.add_argument("--qual", dest="qualFilter",
                         type=float,
                         required=False,
                         help="""Optionally, specify the QUAL value that variants must equal or
                         exceed to be included in the final VCF file (recommended: 30.0)""")
    iparser.add_argument("--missing", dest="missingFilter",
                         type=float,
                         required=False,
                         help="""Optionally, specify the proportion of missing data that is
                         tolerated in both bulk populations before a variant is filtered out
                         (recommended: 0.25)""")
    
    # Depth-subparser arguments
    dparser.add_argument("-f", dest="genomeFasta",
                         required=True,
                         help="""Specify the location of the genome FASTA file that BAM files
                         are aligned to""")
    dparser.add_argument("--windowSize", dest="windowSize",
                         type=int,
                         required=False,
                         help="""Optionally, specify the window size that reads will be
                         binned into for deletion calling (recommended: 1000)""")
    dparser.add_argument("--bam", dest="bamFiles",
                         required=False,
                         nargs="+",
                         help="""Optionally, specify one or more locations of BAM files and/or
                         directories containing BAM files for depth calculations""")
    dparser.add_argument("--bamSuffix", dest="bamSuffix",
                         required=False,
                         help="""Optionally, specify the suffix used to denote BAM files;
                         relevant if directories are provided to --bam""")
    dparser.add_argument("--threads", dest="threads",
                         type=int,
                         required=False,
                         help="""Optionally, specify the number of threads to use for samtools
                         depth calculations (default: 1)""",
                         default=1)
    
    # Call-subparser arguments
    cparser.add_argument("-f", dest="genomeFasta",
                         required=True,
                         help="""Specify the location of the genome FASTA file that BAM files
                         are aligned to""")
    cparser.add_argument("--qual", dest="qualFilter",
                         type=float,
                         required=False,
                         help="""Optionally, specify the QUAL value that variants must equal or
                         exceed to be included in the final VCF file (recommended: 30.0)""")
    cparser.add_argument("--missing", dest="missingFilter",
                         type=float,
                         required=False,
                         help="""Optionally, specify the proportion of missing data that is
                         tolerated in both bulk populations before a variant is filtered out
                         (recommended: 0.25)""")
    cparser.add_argument("--bam", dest="bamFiles",
                         required=False,
                         nargs="+",
                         help="""Optionally, specify one or more locations of BAM files and/or
                         directories containing BAM files for variant calling and/or depth
                         calculations""")
    cparser.add_argument("--bamSuffix", dest="bamSuffix",
                         required=False,
                         help="""Optionally, specify the suffix used to denote BAM files;
                         relevant if directories are provided to --bam""")
    cparser.add_argument("--threads", dest="threads",
                         type=int,
                         required=False,
                         help="""Optionally, specify the number of threads to use for bcftools
                         variant calling (default: 1)""",
                         default=1)
    
    # View-subparser arguments
    # N/A
    
    args = subParentParser.parse_args()
    validate_args(args)
    
    # Split into mode-specific functions
    if args.mode in ["initialise", "init"]:
        print("## psQTL_prep.py - Initialisation ##")
        imain(args)
    elif args.mode == "depth":
        validate_uncached(args)
        print("## psQTL_prep.py - Depth Calculation ##")
        dmain(args)
    elif args.mode == "call":
        validate_uncached(args)
        print("## psQTL_prep.py - Variant Calling ##")
        cmain(args)
    elif args.mode == "view":
        print("## psQTL_prep.py - View Directory ##")
        vmain(args)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def imain(args):
    paramsCache = ParameterCache(args.workingDirectory)
    try:
        paramsCache.initialise(args)
    except FileExistsError:
        paramsCache.merge(args)
    
    print("Initialisation complete!")

def dmain(args):
    DEPTH_SUFFIX = ".depth.tsv"
    DEPTH_DIR = os.path.join(args.workingDirectory, "depth")
    os.makedirs(DEPTH_DIR, exist_ok=True)
    
    # Validate that necessary programs exist
    validate_samtools_exists()
    
    # Merge params and args
    paramsCache = ParameterCache(args.workingDirectory)
    paramsCache.merge(args) # raises FileNotFoundError if cache does not exist
    
    # Validate that necessary arguments are provided
    if args.bamFiles == None or args.bamFiles == []:
        raise ValueError("--bam files not yet provided for depth calculation!")
    if args.windowSize == None:
        raise ValueError("--windowSize not yet provided for depth calculation!")
    
    # Get the sample prefixes from the BAM files
    bamPrefixes = [ os.path.basename(f).rsplit(args.bamSuffix, maxsplit=1)[0] for f in args.bamFiles ]
    
    # Determine which depth files need to be generated
    depthIO = []
    for bamFile, bamPrefix in zip(args.bamFiles, bamPrefixes):
        depthFile = os.path.join(DEPTH_DIR, f"{bamPrefix}{DEPTH_SUFFIX}")
        
        # Skip if the file already exists
        if os.path.isfile(depthFile) and os.path.isfile(depthFile + ".ok"):
            continue
        depthIO.append([bamFile, depthFile])
    
    # Run samtools depth
    run_samtools_depth(depthIO, args.threads)
    
    # Parse the genome FASTA file to get contig lengths
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }
    
    # Determine which depth files need to be binned
    binIO = []
    for bamPrefix in bamPrefixes:
        depthFile = os.path.join(DEPTH_DIR, f"{bamPrefix}{DEPTH_SUFFIX}")
        
        # Error out if depth file is missing
        if not os.path.isfile(depthFile):
            raise FileNotFoundError(f"Depth file '{depthFile}' not found!")
        
        # Format the binned file name for this depth file
        binFile = os.path.join(DEPTH_DIR, f"{bamPrefix}.binned.{args.windowSize}.tsv")
        
        # Skip if the binned file already exists
        if os.path.isfile(binFile) and os.path.isfile(binFile + ".ok"):
            continue
        binIO.append([depthFile, binFile])
    
    # Bin the depth files
    bin_samtools_depth(binIO, lengthsDict, args.threads, args.windowSize)
    
    # Get all sample prefixes and their associated bin file
    samplePairs = []
    for bamPrefix in bamPrefixes:
        binFile = os.path.join(DEPTH_DIR, f"{bamPrefix}.binned.{args.windowSize}.tsv")
        
        # Error out if bin file is missing
        if not os.path.isfile(binFile):
            raise FileNotFoundError(f"Bin file '{binFile}' not found!")
        
        # Store file paied with sample prefix
        samplePairs.append([bamPrefix, binFile])
    
    # Collate binned files into a VCF-like format of deletion calls
    FINAL_DELETION_FILE = os.path.join(DEPTH_DIR, "psQTL_deletions.vcf")
    if (not os.path.isfile(FINAL_DELETION_FILE)) or (not os.path.isfile(FINAL_DELETION_FILE + ".ok")):
        print("# Generating deletion file...")
        call_deletions_from_depth(samplePairs, FINAL_DELETION_FILE, args.windowSize)
        open(FINAL_DELETION_FILE + ".ok", "w").close() # touch a .ok file to indicate success
    else:
        print(f"# Deletion file '{FINAL_DELETION_FILE}' exists; skipping ...")
    
    # Update param cache with (potentially) newly produced deletion file
    paramsCache = ParameterCache(args.workingDirectory)
    paramsCache.load() # reload in case we're running call simultaneously
    paramsCache.deletionFile = FINAL_DELETION_FILE
    paramsCache.windowSize = args.windowSize
    
    print("Depth file generation complete!")

def cmain(args):
    '''
    One note of caution: this script uses '.ok' files to track files that we know have
    been successfully processed. However, we do not apply the same for VCF indices and
    FASTA indices. As such, any unusual issues or bugs that occur may well be related to
    an index that was partially created. It is expected that stderr values will provide
    this information but this hasn't been validated yet.
    '''
    CALL_DIR = os.path.join(args.workingDirectory, "call")
    os.makedirs(CALL_DIR, exist_ok=True)
    
    # Validate that necessary programs exist
    validate_samtools_exists()
    validate_bcftools_exists()
    validate_bgzip_exists()
    validate_vt_exists()
    
    # Merge params and args
    paramsCache = ParameterCache(args.workingDirectory)
    paramsCache.merge(args) # raises FileNotFoundError if cache does not exist
    
    # Validate that necessary parameters arguments are provided
    if args.bamFiles == None or args.bamFiles == []:
        raise ValueError("--bam files not yet provided for variant calling!")
    if args.qualFilter == None:
        raise ValueError("--qual not yet provided for variant calling!")
    if args.missingFilter == None:
        raise ValueError("--missing not yet provided for variant calling!")
    
    # Index the reference genome (if necessary)
    if not os.path.isfile(args.genomeFasta + ".fai"):
        run_samtools_faidx(args.genomeFasta)
    
    # Create a bamlist file
    BAMLIST_FILE = os.path.join(CALL_DIR, "bamlist.txt")
    with open(BAMLIST_FILE, "w") as bamlistFile: # allowed to overwrite existing files
        for bamFile in args.bamFiles:
            bamlistFile.write(f"{bamFile}\n")
    
    # Run bcftools mpileup->call on each contig
    run_bcftools_call(BAMLIST_FILE, args.genomeFasta, CALL_DIR, args.threads) # handles skipping internally
    
    # Index each VCF file
    for vcfFile in [ os.path.join(CALL_DIR, f) for f in os.listdir(CALL_DIR) ]:
        if vcfFile.endswith(".vcf.gz") and not os.path.isfile(vcfFile + ".csi"):
            run_bcftools_index(vcfFile)
    
    # Run normalisation on each contig's VCF
    run_normalisation(args.genomeFasta, CALL_DIR, args.threads) # handles skipping internally
    
    # Index each VCF file
    for vcfFile in [ os.path.join(CALL_DIR, f) for f in os.listdir(CALL_DIR) ]:
        if vcfFile.endswith(".vcf.gz") and not os.path.isfile(vcfFile + ".csi"):
            run_bcftools_index(vcfFile)
    
    # Concatenate all VCF files
    CONCAT_VCF_FILE = os.path.join(CALL_DIR, "psQTL_variants.vcf.gz")
    if (not os.path.isfile(CONCAT_VCF_FILE)) or (not os.path.isfile(CONCAT_VCF_FILE + ".ok")):
        print("# Concatenating VCF files...")
        run_bcftools_concat(args.genomeFasta, CALL_DIR, CONCAT_VCF_FILE)
    else:
        print(f"# Concatenated VCF file '{CONCAT_VCF_FILE}' exists; skipping ...")
    
    # Index the concatenated VCF file
    if not os.path.isfile(CONCAT_VCF_FILE + ".csi"):
        run_bcftools_index(CONCAT_VCF_FILE)
    
    # Update param cache with newly produced VCF file
    paramsCache = ParameterCache(args.workingDirectory)
    paramsCache.load() # reload in case we're running depth simultaneously
    paramsCache.vcfFile = CONCAT_VCF_FILE
    
    # Filter the VCF file
    FILTERED_FILE = os.path.join(CALL_DIR, "psQTL_variants.filtered.vcf")
    FINAL_FILTERED_FILE = FILTERED_FILE + ".gz"
    if (not os.path.isfile(FINAL_FILTERED_FILE)) or (not os.path.isfile(FINAL_FILTERED_FILE + ".ok")):        
        # Skip filtering if metadata is not initialised
        if args.metadataFile == None:
            print("# Metadata file not initialised; skipping final VCF filtering...")
            print("# Run 'psQTL_prep.py initialise' to set the metadata file and you can filter the VCF later.")
        
        # Raise error if metadata is initialised, but file is now missing
        elif not os.path.isfile(args.metadataFile):
            raise FileNotFoundError(f"Metadata file '{args.metadataFile}' can no longer be found!")
        
        # Otherwise, proceed with filtering
        else:
            # Parse the metadata file
            metadataDict = parse_metadata(args.metadataFile)
            
            # Filter the VCF file
            print("# Filtering VCF file...")
            filter_vcf(CONCAT_VCF_FILE, FILTERED_FILE, metadataDict, args.missingFilter, args.qualFilter)
            
            # bgzip the filtered VCF file
            print("# bgzipping filtered VCF file...")
            run_bgzip(FILTERED_FILE)
            
            # Index the concatenated VCF file
            run_bcftools_index(FINAL_FILTERED_FILE)
            open(FINAL_FILTERED_FILE + ".ok", "w").close() # touch a .ok file to indicate success
    else:
        print(f"# Filtered VCF file '{FINAL_FILTERED_FILE}' exists; skipping ...")
    
    # Update param cache with (potentially) newly produced filtered VCF file
    paramsCache = ParameterCache(args.workingDirectory)
    paramsCache.load() # reload in case we're running depth simultaneously
    paramsCache.filteredVcfFile = FINAL_FILTERED_FILE
    
    print("Variant calling complete!")

def vmain(args):
    paramsCache = ParameterCache(args.workingDirectory)
    paramsCache.merge(args) # raises FileNotFoundError if cache does not exist
    
    # Present standard parameters
    print("# Parameters:")
    print(f"Working directory: {args.workingDirectory}")
    
    if args.bamFiles != []:
        print(f"BAM files: {args.bamFiles}")
    else:
        print("BAM files: None")
    
    # Present metadata cache
    print() # blank line for spacing
    print("# Metadata details:")
    if args.metadataFile is not None:
        metadataCache = MetadataCache(args.workingDirectory)
        metadataCache.establish()
        if metadataCache.metadataFile == None:
            print("## Metadata cache not found; re-initialising...")
            metadataCache.metadataFile = args.metadataFile
        
        print(f"Metadata file: {args.metadataFile}")
        print(f"Bulk 1 samples (n={len(metadataCache.bulk1)}): {metadataCache.bulk1}")
        print(f"Bulk 2 samples (n={len(metadataCache.bulk2)}): {metadataCache.bulk2}")
    else:
        print("Metadata file: None")
    
    # Present VCF cache
    print() # blank line for spacing
    print("# Variants VCF details:")
    if args.vcfFile is not None:
        vcfCache = VcfCache(args.workingDirectory)
        vcfCache.establish()
        if vcfCache.vcfFile == None:
            print("## VCF cache not found; re-initialising...")
            vcfCache.vcfFile = args.vcfFile
        
        print(f"VCF file: {args.vcfFile}")
        print(f"Num. variants: {vcfCache.variants}")
        print(f"Samples (n={len(vcfCache.samples)}): {vcfCache.samples}")
        print(f"Contigs (n={len(vcfCache.contigs)}): {vcfCache.contigs}")
    else:
        print("VCF file: None")
    
    print() # blank line for spacing
    print("# Filtered variants VCF details:")
    if args.filteredVcfFile is not None:
        vcfCache = VcfCache(args.workingDirectory)
        vcfCache.establish()
        if vcfCache.filteredVcfFile == None:
            print("## VCF cache not found; re-initialising...")
            vcfCache.filteredVcfFile = args.filteredVcfFile
        
        print(f"Filtered VCF file: {args.filteredVcfFile}")
        if paramsCache.qualFilter == None:
            print(f"Quality filter: unknown")
        else:
            print(f"Quality filter: {paramsCache.qualFilter}")
        if paramsCache.missingFilter == None:
            print(f"Missing filter: unknown")
        else:
            print(f"Missing filter: {paramsCache.missingFilter}")
        print(f"Num. filtered variants: {vcfCache.filteredVariants}")
        print(f"Filtered samples (n={len(vcfCache.filteredSamples)}): {vcfCache.filteredSamples}")
        print(f"Filtered contigs (n={len(vcfCache.filteredContigs)}): {vcfCache.filteredContigs}")
    else:
        print("Filtered VCF file: None")
    print() # blank line for spacing
    
    # Present deletion cache
    print("# Deletion VCF-like details:")
    if args.deletionFile is not None:
        deletionCache = DeletionCache(args.workingDirectory)
        deletionCache.establish()
        if deletionCache.deletionFile == None:
            print("## Deletion cache not found; re-initialising...")
            deletionCache.deletionFile = args.deletionFile
        
        print(f"Deletion file: {args.deletionFile}")
        if paramsCache.windowSize == None:
            print(f"Window size: unknown")
        else:
            print(f"Window size: {paramsCache.windowSize} bp")
        print(f"Total num. bins: {deletionCache.bins}")
        print(f"Num. bins with deletion: {deletionCache.deletionBins}")
        print(f"Samples (n={len(deletionCache.samples)}): {deletionCache.samples}")
        print(f"Contigs (n={len(deletionCache.contigs)}): {deletionCache.contigs}")
        
    else:
        print("Deletion file: None")
    print() # blank line for spacing
    
    # Identify potential conflicts or issues
    issues = []
    
    ## Sample issues
    keepFindingIssues = True
    try:
        metaSamples = set(metadataCache.bulk1 + metadataCache.bulk2)
        if len(metaSamples) == 2:
            issues.append("Metadata only indicates two samples; analysis interpretation may be limited")
    except:
        issues.append("Metadata file needs to be set before further analysis can proceed")
        keepFindingIssues = False
    
    if keepFindingIssues:
        try:
            vcfSamples = set(vcfCache.filteredSamples) if vcfCache.filteredSamples != None else set(vcfCache.samples)
            if metaSamples != vcfSamples:
                issues.append(f"Metadata samples (n={len(metaSamples)}) do not match VCF samples (n={len(vcfCache.samples)})")
        except:
            issues.append("VCF file may need to be set or generated before further analysis can proceed")
            keepFindingIssues = False
    
    if keepFindingIssues:
        try:
            deletionSamples = set(deletionCache.samples)
            if metaSamples != deletionSamples:
                issues.append(f"Metadata samples (n={len(metaSamples)}) do not match deletion samples (n={len(deletionCache.samples)})")
            if vcfSamples != deletionSamples:
                issues.append(f"VCF samples (n={len(vcfSamples)}) do not match deletion samples (n={len(deletionSamples)})")
        except:
            issues.append("Deletion file may need to be set or generated before further analysis can proceed")
            keepFindingIssues = False
    
    ## Contig issues
    if keepFindingIssues:
        vcfContigs = set(vcfCache.filteredContigs) if vcfCache.filteredContigs != None else set(vcfCache.contigs)
        deletionContigs = set(deletionCache.contigs)
        if vcfContigs != deletionContigs:
            issues.append(f"VCF contigs (n={len(vcfContigs)}) do not match deletion contigs (n={len(deletionContigs)})")
    
    # Present potential issues
    print("# Potential issues:")
    if issues == []:
        print("No issues found!")
    else:
        print("\n".join(issues))
    print() # blank line for spacing
    
    print("Analysis viewing complete!")

if __name__ == "__main__":
    main()
