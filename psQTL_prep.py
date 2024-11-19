#! python3
# psQTL_prep.py
# Represents step 1 in the psQTL pipeline, where the user can initialise a working directory,
# view the metadata of an analysis directory, call variants, and calculate depth files for samples.
# Called variants and depth files can be used as input to the psQTL_proc.py script.

import os, argparse, sys
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.cache_handling import initialise_param_cache, update_param_cache, load_param_cache, \
                            load_vcf_cache, update_vcf_cache, initialise_deletion_cache, \
                            load_deletion_cache, initialise_metadata_cache, load_metadata_cache, \
                            merge_cache_into_args
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
        paramsDict = load_param_cache(args.workingDirectory)
        if paramsDict == {}:
            raise FileNotFoundError("Working directory has not been initialised;" +
                                    f"parameter cache not found in '{args.workingDirectory}'!")

def validate_bam_files(args):
    '''
    Performs validations to make sure that 1) any indicated files exists, 2) any directories
    provided contain BAM files, and 3) that no duplicate prefixes are found, as that may
    cause issues when producing VCF and depth files.
    '''
    # Validate BAM files
    bamPrefixes = set()
    foundBAMs = []
    for location in args.bamFiles:
        location = os.path.abspath(location)
        
        # Handle an existing file
        if os.path.isfile(location):
            if not location.endswith(args.bamSuffix):
                raise ValueError(f"BAM file '{location}' does not end with the specified suffix '{args.bamSuffix}'")
            else:
                bamPrefix = os.path.basename(location).rsplit(args.bamSuffix, maxsplit=1)[0]
                if bamPrefix in bamPrefixes:
                    raise ValueError(f"Duplicate BAM prefix found: '{bamPrefix}'")
                
                foundBAMs.append(location)
        # Handle an existing directory
        elif os.path.isdir(location):
            foundAny = False
            for f in os.listdir(location):
                if f.endswith(args.bamSuffix):
                    bamPrefix = os.path.basename(location).rsplit(args.bamSuffix, maxsplit=1)[0]
                    if bamPrefix in bamPrefixes:
                        raise ValueError(f"Duplicate BAM prefix found: '{bamPrefix}'")
                    
                    foundBAMs.append(os.path.join(location, f))
                    foundAny = True
            if not foundAny:
                raise FileNotFoundError(f"No BAM files found in directory '{location}' ending with '{args.bamSuffix}'")
        # Error out if location does not exist
        else:
            raise FileNotFoundError(f"Input BAM file or directory '{location}' not found!")
    args.bamFiles = foundBAMs

def validate_i(args):
    '''
    The goal of this function is to make sure that at least one input was given, and that
    any provided files actually exist.
    '''
    # Check that at least one input was given
    if args.metadataFile == None and args.vcfFile == None and args.bamFiles == [] \
    and args.deletionFile == []:
        raise ValueError("At least one of --meta, --vcf, --deletion, or --bam must be specified!")
    
    # Check that we some sort of useful file input
    if args.bamFiles == None and args.vcfFile == None and args.deletionFile == None:
        raise ValueError("Downstream psQTL processing requires at least one of --vcf, --deletion, or --bam!")
    
    # Validate metadata file
    if args.metadataFile is not None:
        args.metadataFile = os.path.abspath(args.metadataFile)
        if not os.path.isfile(args.metadataFile):
            raise FileNotFoundError(f"Metadata file '{args.metadataFile}' does not exist!")
    
    # Validate VCF file
    if args.vcfFile is not None:
        args.vcfFile = os.path.abspath(args.vcfFile)
        if not os.path.isfile(args.vcfFile):
            raise FileNotFoundError(f"VCF file '{args.vcfFile}' does not exist!")
    
    # Validate deletion file
    if args.deletionFile is not None:
        args.deletionFile = os.path.abspath(args.deletionFile)
        if not os.path.isfile(args.deletionFile):
            raise FileNotFoundError(f"Deletion VCF-like file '{args.deletionFile}' does not exist!")
    
    # Validate BAM files
    if args.bamFiles != []:
        validate_bam_files(args)

def validate_d(args):
    '''
    The goal of this function is to determine whether this submodule's depth calculations
    are necessary i.e., an existing deletion file is not provided. If so, BAM file locations
    are verified to ensure that they exist in order for depth calculations to proceed.
    '''
    # Validate threads
    if args.threads < 1:
        raise ValueError("Number of threads must be at least 1!")
    
    # Validate window size
    if args.windowSize < 1:
        raise ValueError("Window size must be at least 1!")
    
    # Validate BAM files
    if args.bamFiles == []:
        raise ValueError("No BAM files have yet been provided for depth calculations!")
    else:
        validate_bam_files(args)
    
    # Validate genome FASTA file
    args.genomeFasta = os.path.abspath(args.genomeFasta)
    if not os.path.isfile(args.genomeFasta):
        raise FileNotFoundError(f"Genome FASTA file '{args.genomeFasta}' does not exist!")

def validate_c(args):
    # Validate threads
    if args.threads < 1:
        raise ValueError("Number of threads must be at least 1!")
    
    # Validate qual and missing filters
    if args.qualFilter < 0:
        raise ValueError("QUAL filter value must be at least 0!")
    if args.missingFilter < 0:
        raise ValueError("Missing filter value must be at least 0!")
    elif args.missingFilter > 1:
        raise ValueError("Missing filter value must be at most 1!")
    
    # Validate BAM files
    if args.bamFiles == []:
        raise ValueError("No BAM files have yet been provided for variant calling!")
    else:
        validate_bam_files(args)
    
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
                         columns indicating 1) sample ID and 2) the bulk it belongs to""",
                         default=None)
    iparser.add_argument("--vcf", dest="vcfFile",
                         required=False,
                         help="""Optionally, specify a VCF file containg per-sample variant
                         calls that you have already produced""",
                         default=None)
    iparser.add_argument("--deletion", dest="deletionFile",
                         required=False,
                         help="""Optionally, specify a deletion VCF-like file containing per-sample
                         deletion calls that you have already produced""",
                         default=None)
    iparser.add_argument("--bam", dest="bamFiles",
                         required=False,
                         nargs="+",
                         help="""Optionally, specify one or more locations of BAM files and/or
                         directories containing BAM files for variant calling and/or depth
                         calculations""",
                         default=[])
    iparser.add_argument("--bamSuffix", dest="bamSuffix",
                         required=False,
                         help="""Optionally, specify the suffix used to denote BAM files;
                         relevant if directories are provided to --bam (default: '.bam')""",
                         default=".bam")
    
    # Depth-subparser arguments
    dparser.add_argument("-f", dest="genomeFasta",
                         required=True,
                         help="""Specify the location of the genome FASTA file that BAM files
                         are aligned to""")
    dparser.add_argument("--windowSize", dest="windowSize",
                         type=int,
                         required=False,
                         help="""Optionally, specify the window size that reads will be
                         binned into for deletion calling""",
                         default=1000)
    dparser.add_argument("--bam", dest="bamFiles",
                         required=False,
                         nargs="+",
                         help="""Optionally, specify one or more locations of BAM files and/or
                         directories containing BAM files for depth calculations""",
                         default=[])
    dparser.add_argument("--bamSuffix", dest="bamSuffix",
                         required=False,
                         help="""Optionally, specify the suffix used to denote BAM files;
                         relevant if directories are provided to --bam (default: '.bam')""",
                         default=[])
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
                         exceed to be included in the final VCF file (default: 30.0)""",
                         default=30.0)
    cparser.add_argument("--missing", dest="missingFilter",
                         type=float,
                         required=False,
                         help="""Optionally, specify the proportion of missing data that is
                         tolerated in both bulk populations before a variant is filtered out
                         (default: 0.25)""",
                         default=0.25)
    cparser.add_argument("--bam", dest="bamFiles",
                         required=False,
                         nargs="+",
                         help="""Optionally, specify one or more locations of BAM files and/or
                         directories containing BAM files for variant calling and/or depth
                         calculations""",
                         default=[])
    cparser.add_argument("--bamSuffix", dest="bamSuffix",
                         required=False,
                         help="""Optionally, specify the suffix used to denote BAM files;
                         relevant if directories are provided to --bam (default: '.bam')""",
                         default=[])
    cparser.add_argument("--threads", dest="threads",
                         type=int,
                         required=False,
                         help="""Optionally, specify the number of threads to use for bcftools
                         variant calling (default: 1)""",
                         default=1)
    
    args = subParentParser.parse_args()
    validate_args(args)
    merge_cache_into_args(args)
    
    # Split into mode-specific functions
    if args.mode == "depth":
        print("## psQTL_prep.py - Depth Calculation ##")
        dmain(args)
    elif args.mode == "call":
        print("## psQTL_prep.py - Variant Calling ##")
        cmain(args)
    elif args.mode in ["initialise", "init"]:
        print("## psQTL_prep.py - Initialisation ##")
        imain(args)
    elif args.mode == "view":
        print("## psQTL_prep.py - View Directory ##")
        vmain(args)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def imain(args):
    validate_i(args)
    
    # Initialise or update the param cache
    initialise_param_cache(args)
    
    # Initialise or update other caches
    if args.metadataFile != None:
        initialise_metadata_cache(args.workingDirectory)
    if args.deletionFile != None:
        "If we receive a pre-computed deletion file, our windowSize value will be unknown"
        initialise_deletion_cache(args.workingDirectory, None)
    if args.vcfFile != None:
        update_vcf_cache(args.workingDirectory, None, None)
    print("Initialisation complete!")

def dmain(args):
    DEPTH_SUFFIX = ".depth.tsv"
    validate_d(args)
    
    # Validate that necessary programs exist
    validate_samtools_exists()
    
    # Create depth directory if it doesn't exist
    DEPTH_DIR = os.path.join(args.workingDirectory, "depth")
    os.makedirs(DEPTH_DIR, exist_ok=True)
    
    # Store any newly provided values in the cache
    update_param_cache(args.workingDirectory, {"bamFiles" : args.bamFiles,
                                               "bamSuffix": args.bamSuffix})
    
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
        
        # Update the param and deletion caches
        update_param_cache(args.workingDirectory, {"deletionFile" : FINAL_DELETION_FILE})
        initialise_deletion_cache(args.workingDirectory, args.windowSize)
    else:
        print(f"# Deletion file '{FINAL_DELETION_FILE}' exists; skipping ...")
    
    print("Depth file generation complete!")

def cmain(args):
    '''
    One note of caution: this script uses '.ok' files to track files that we know have
    been successfully processed. However, we do not apply the same for VCF indices and
    FASTA indices. As such, any unusual issues or bugs that occur may well be related to
    an index that was partially created. It is expected that stderr values will provide
    this information but this hasn't been validated yet.
    '''
    validate_c(args)
    
    # Validate that necessary programs exist
    validate_samtools_exists()
    validate_bcftools_exists()
    validate_bgzip_exists()
    validate_vt_exists()
    
    # Create call directory if it doesn't exist
    CALL_DIR = os.path.join(args.workingDirectory, "call")
    os.makedirs(CALL_DIR, exist_ok=True)
    
    # Store any newly provided values in the cache
    update_param_cache(args.workingDirectory, {"bamFiles" : args.bamFiles,
                                               "bamSuffix": args.bamSuffix})
    
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
    
    # Update the param and VCF caches
    update_param_cache(args.workingDirectory, {"vcfFile" : CONCAT_VCF_FILE})
    update_vcf_cache(args.workingDirectory, None, None)
    
    # Filter the VCF file
    FILTERED_FILE = os.path.join(CALL_DIR, "psQTL_variants.filtered.vcf")
    FINAL_FILTERED_FILE = FILTERED_FILE + ".gz"
    if (not os.path.isfile(FINAL_FILTERED_FILE)) or (not os.path.isfile(FINAL_FILTERED_FILE + ".ok")):
        paramsDict = load_param_cache(args.workingDirectory)
        
        # Skip filtering if metadata is not initialised
        if paramsDict['metadataFile'] == None:
            print("# Metadata file not initialised; skipping final VCF filtering...")
            print("# Run 'psQTL_prep.py initialise' to set the metadata file and you can filter the VCF later.")
        
        # Raise error if metadata is initialised, but file is now missing
        elif not os.path.isfile(paramsDict['metadataFile']):
            raise FileNotFoundError(f"Metadata file '{paramsDict['metadataFile']}' is specified in your " +
                                    "parameters cache but can no longer be found!")
        
        # Otherwise, proceed with filtering
        else:
            # Parse the metadata file
            metadataDict = parse_metadata(paramsDict['metadataFile'])
            
            # Filter the VCF file
            print("# Filtering VCF file...")
            filter_vcf(CONCAT_VCF_FILE, FILTERED_FILE, metadataDict, args.missingFilter, args.qualFilter)
            
            # bgzip the filtered VCF file
            print("# bgzipping filtered VCF file...")
            run_bgzip(FILTERED_FILE)
            
            # Index the concatenated VCF file
            run_bcftools_index(FINAL_FILTERED_FILE)
            open(FINAL_FILTERED_FILE + ".ok", "w").close() # touch a .ok file to indicate success
            
            # Update the param and VCF caches
            update_param_cache(args.workingDirectory, {"filteredVcfFile" : FINAL_FILTERED_FILE})
            update_vcf_cache(args.workingDirectory, args.qualFilter, args.missingFilter)
    else:
        print(f"# Filtered VCF file '{FINAL_FILTERED_FILE}' exists; skipping ...")
    
    print("Variant calling complete!")

def vmain(args):
    # Handle argument parsing and cache formatting
    paramsDict = load_param_cache(args.workingDirectory)
    
    # Present standard parameters
    print("# Parameters:")
    print(f"Working directory: {paramsDict['workingDirectory']}")
    
    if paramsDict['bamFiles'] != []:
        print(f"BAM files: {paramsDict['bamFiles']}")
    else:
        print("BAM files: None")
    
    # Present metadata cache
    print() # blank line for spacing
    print("# Metadata details:")
    if paramsDict['metadataFile'] is not None:
        metadataDict = load_metadata_cache(args.workingDirectory)
        if metadataDict == {}:
            initialise_metadata_cache(args.workingDirectory)
            metadataDict = load_metadata_cache(args.workingDirectory)
        
        print(f"Metadata file: {paramsDict['metadataFile']}")
        print(f"Bulk 1 samples (n={len(metadataDict['bulk1'])}): {metadataDict['bulk1']}")
        print(f"Bulk 2 samples (n={len(metadataDict['bulk2'])}): {metadataDict['bulk2']}")
    else:
        print("Metadata file: None")
    
    # Present VCF cache
    print() # blank line for spacing
    print("# Variants VCF details:")
    if paramsDict['vcfFile'] is not None:
        vcfDict = load_vcf_cache(args.workingDirectory)
        if vcfDict == {}:
            print("## VCF cache not found; re-initialising...")
            update_vcf_cache(args.workingDirectory, None, None)
            vcfDict = load_vcf_cache(args.workingDirectory)
        
        print(f"VCF file: {paramsDict['vcfFile']}")
        print(f"Num. variants: {vcfDict['variants']}")
        print(f"Samples (n={len(vcfDict['samples'])}): {vcfDict['samples']}")
        print(f"Contigs (n={len(vcfDict['contigs'])}): {vcfDict['contigs']}")
    else:
        print("VCF file: None")
    
    if paramsDict['filteredVcfFile'] is not None:
        vcfDict = load_vcf_cache(args.workingDirectory)
        if vcfDict == {}: # this shouldn't happen unless the user is messing with the cache
            print("## VCF cache not found; re-initialising...")
            update_vcf_cache(args.workingDirectory, None, None)
            vcfDict = load_vcf_cache(args.workingDirectory)
        
        print(f"Filtered VCF file: {paramsDict['filteredVcfFile']}")
        print(f"Num. filtered variants: {vcfDict['filteredVariants']}")
        print(f"Filtered samples (n={len(vcfDict['filteredSamples'])}): {vcfDict['filteredSamples']}")
        print(f"Filtered contigs (n={len(vcfDict['filteredContigs'])}): {vcfDict['filteredContigs']}")
    else:
        print("Filtered VCF file: None")
    print() # blank line for spacing
    
    # Present deletion cache
    print("# Deletion VCF-like details:")
    if paramsDict['deletionFile'] is not None:
        deletionDict = load_deletion_cache(args.workingDirectory)
        if deletionDict == {}:
            print("## Deletion cache not found; re-initialising...")
            initialise_deletion_cache(args.workingDirectory, None)
            deletionDict = load_deletion_cache(args.workingDirectory)
        
        print(f"Deletion file: {paramsDict['deletionFile']}")
        if deletionDict['windowSize'] == None:
            print(f"Window size: unknown")
        else:
            print(f"Window size: {deletionDict['windowSize']} bp")
        print(f"Total num. bins: {deletionDict['totalBins']}")
        print(f"Num. bins with deletion: {deletionDict['deletionBins']}")
        print(f"Samples (n={len(deletionDict['samples'])}): {deletionDict['samples']}")
        print(f"Contigs (n={len(deletionDict['contigs'])}): {deletionDict['contigs']}")
        
    else:
        print("Deletion file: None")
    print() # blank line for spacing
    
    # Identify potential conflicts or issues
    issues = []
    
    ## Sample issues
    keepFindingIssues = True
    try:
        metaSamples = set(metadataDict['bulk1'] + metadataDict['bulk2'])
        if len(metaSamples) == 2:
            issues.append("Metadata only indicates two samples; analysis interpretation may be limited")
    except:
        issues.append("Metadata file needs to be set before further analysis can proceed")
        keepFindingIssues = False
    
    if keepFindingIssues:
        try:
            vcfSamples = set(vcfDict['filteredSamples']) if vcfDict['filteredSamples'] != None else set(vcfDict['samples'])
            if metaSamples != vcfSamples:
                issues.append(f"Metadata samples (n={len(metaSamples)}) do not match VCF samples (n={len(vcfDict['samples'])})")
        except:
            issues.append("VCF file may need to be set or generated before further analysis can proceed")
            keepFindingIssues = False
    
    if keepFindingIssues:
        try:
            deletionSamples = set(deletionDict['samples'])
            if metaSamples != deletionSamples:
                issues.append(f"Metadata samples (n={len(metaSamples)}) do not match deletion samples (n={len(deletionDict['samples'])})")
            if vcfSamples != deletionSamples:
                issues.append(f"VCF samples (n={len(vcfSamples)}) do not match deletion samples (n={len(deletionSamples)})")
        except:
            issues.append("Deletion file may need to be set or generated before further analysis can proceed")
            keepFindingIssues = False
    
    ## Contig issues
    if keepFindingIssues:
        vcfContigs = set(vcfDict['filteredContigs']) if vcfDict['filteredContigs'] != None else set(vcfDict['contigs'])
        deletionContigs = set(deletionDict['contigs'])
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
