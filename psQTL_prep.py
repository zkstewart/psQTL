#! python3
# psQTL_prep.py
# Represents step 1 in the psQTL pipeline, where the user can initialise a working directory,
# view the metadata of an analysis directory, call variants, and calculate depth files for samples.
# Called variants and depth files can be used as input to the psQTL_proc.py script.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.cache_handling import initialise_param_cache, update_param_cache, load_param_cache, \
                            load_vcf_cache, initialise_vcf_cache
from modules.samtools_handling import validate_samtools_exists, run_samtools_depth, run_samtools_index
from modules.bcftools_handling import validate_bcftools_exists, validate_bgzip_exists, validate_vt_exists, \
                                      run_bcftools_call, run_bcftools_index, run_normalisation, run_bcftools_concat

def validate_args(args):
    # Validate working directory
    args.workingDirectory = os.path.abspath(args.workingDirectory)
    if not os.path.exists(args.workingDirectory):
        if not os.path.exists(os.path.dirname(args.workingDirectory)):
            raise FileNotFoundError(f"Parent directory for proposed -w location " + 
                                    f"'{os.path.dirname(args.workingDirectory)}' does not exist!")
        os.mkdir(args.workingDirectory)
        print(f"# Created working directory '{args.workingDirectory}'")

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
        
        if os.path.isfile(location) and location.endswith(args.bamSuffix):
            bamPrefix = os.path.basename(location).rsplit(args.bamSuffix, maxsplit=1)[0]
            if bamPrefix in bamPrefixes:
                raise ValueError(f"Duplicate BAM prefix found: '{bamPrefix}'")
            
            foundBAMs.append(location)
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
        else:
            raise FileNotFoundError(f"Input BAM file or directory '{location}' not found!")
    args.bamFiles = foundBAMs

def validate_depth_files(args):
    '''
    Performs validations to make sure that 1) any indicated files exists, 2) any directories
    provided contain depth files, and 3) that no duplicate prefixes are found, as that may
    cause issues when analysing depth files.
    '''
    # Validate depth files
    depthPrefixes = set()
    foundDepths = []
    for location in args.depthFiles:
        location = os.path.abspath(location)
        
        if os.path.isfile(location) and location.endswith(args.depthSuffix):
            depthPrefix = os.path.basename(location).rsplit(args.depthSuffix, maxsplit=1)[0]
            if depthPrefix in depthPrefixes:
                raise ValueError(f"Duplicate depth prefix found: '{depthPrefix}'")
            
            foundDepths.append(location)
        elif os.path.isdir(location):
            foundAny = False
            for f in os.listdir(location):
                if f.endswith(args.depthSuffix):
                    depthPrefix = os.path.basename(location).rsplit(args.depthSuffix, maxsplit=1)[0]
                    if depthPrefix in depthPrefixes:
                        raise ValueError(f"Duplicate depth prefix found: '{depthPrefix}'")
                    
                    foundDepths.append(os.path.join(location, f))
                    foundAny = True
            if not foundAny:
                raise FileNotFoundError(f"No depth files found in directory '{location}' ending with '{args.depthSuffix}'")
        else:
            raise FileNotFoundError(f"Input depth file or directory '{location}' not found!")
    args.depthFiles = foundDepths

def validate_i(args):
    '''
    The goal of this function is to make sure that at least one input was given, and that
    any provided files actually exist.
    '''
    # Check that at least one input was given
    if args.metadataFile is None and args.vcfFile is None and args.bamFiles == [] and args.depthFiles == []:
        raise ValueError("At least one of --meta, --vcf, --bam, or --depth must be specified!")
    
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
    
    # Validate BAM files
    validate_bam_files(args)
    
    # Validate depth files
    validate_depth_files(args)

def validate_d(args, paramsDict):
    '''
    The goal of this function is to make sure that either .bamFiles or .depthFiles are set,
    but not both. That gives us an easy way to determine whether we're producing depth files
    with samtools depth or using existing depth files, for the subsequent step of predicting
    deletion homozygosity/heterozygosity/non-occurrence.
    '''
    # Validate threads
    if args.threads < 1:
        raise ValueError("Number of threads must be at least 1!")
    
    # Validate logic of providing --bam and --depth files
    if args.bamFiles != [] and args.depthFiles != []:
        raise ValueError("Both --bam and --depth files have been provided; please provide only one!")
    
    # Validate logic of provided neither --bam nor --depth files
    elif args.bamFiles == [] and args.depthFiles == []:
        # Check if both exist in the cache
        if paramsDict['bamFiles'] != [] and paramsDict['depthFiles'] != []:
            print("# Found existing BAM and depth files; will use the depth files and ignore BAMs")
            args.depthFiles = paramsDict['depthFiles']
            validate_depth_files(args)
        # Check if only depth files exist in the cache
        elif paramsDict['depthFiles'] != []:
            print("# Found existing depth files; will use these for calculations")
            args.depthFiles = paramsDict['depthFiles']
            validate_depth_files(args)
        # Check if only BAM files exist in the cache
        elif paramsDict['bamFiles'] != []:
            print("# Found existing BAM files; will use these for depth calculations")
            args.bamFiles = paramsDict['bamFiles']
            validate_bam_files(args)
        else:
            raise ValueError("No BAM or depth files have yet been provided for calculation!")
    
    # Validate logic of providing --depth files
    elif args.depthFiles != []:
        if paramsDict['depthFiles'] != []:
            print("# Found existing depth files; will use the newly provided --depth files instead")
        validate_depth_files(args)
    
    # Validate logic of providing --bam files
    else:
        if paramsDict['bamFiles'] != []:
            print("# Found existing BAM files; will use the newly provided --bam files instead")
        validate_bam_files(args)

def validate_c(args, paramsDict):
    # Validate threads
    if args.threads < 1:
        raise ValueError("Number of threads must be at least 1!")
    
    # Validate BAM files
    if paramsDict['bamFiles'] == [] and args.bamFiles == []:
        raise ValueError("No BAM files have yet been provided for variant calling!")
    elif args.bamFiles != []:
        if paramsDict['bamFiles'] != []:
            print(f"# Found existing BAM files: {paramsDict['bamFiles']}")
            print("# Will use new BAM files provided to --bam...")
        validate_bam_files(args)
    else:
        args.bamFiles = paramsDict['bamFiles']
        for bamFile in args.bamFiles:
            if not os.path.isfile(bamFile):
                raise FileNotFoundError(f"BAM file '{bamFile}' in param cache no longer exists!")
            elif not bamFile.endswith(args.bamSuffix):
                raise ValueError(f"BAM file '{bamFile}' does not end with the specified suffix '{args.bamSuffix}'")
    
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
    iparser.add_argument("--depth", dest="depthFiles",
                        required=False,
                        nargs="+",
                        help="""Optionally, specify one or more locations of samtools depth
                        files and/or directories containing depth files""",
                        default=[])
    iparser.add_argument("--depthSuffix", dest="depthSuffix",
                        required=False,
                        help="""Optionally, specify the suffix used to denote samtools depth files;
                        relevant if directories are provided to --depth (default: '.tsv')""",
                        default=".tsv")
    
    # Depth-subparser arguments
    dparser.add_argument("--depth", dest="depthFiles",
                        required=False,
                        nargs="+",
                        help="""Optionally, specify one or more locations of samtools depth
                        files and/or directories containing depth files""",
                        default=[])
    dparser.add_argument("--depthSuffix", dest="depthSuffix",
                         required=False,
                         help="""Optionally, specify the suffix used to denote depth files;
                         (default: '.depth.tsv')""",
                         default=".depth.tsv")
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
                         default=".bam")
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
                         default=".bam")
    cparser.add_argument("--threads", dest="threads",
                         type=int,
                         required=False,
                         help="""Optionally, specify the number of threads to use for bcftools
                         variant calling (default: 1)""",
                         default=1)
    
    args = subParentParser.parse_args()
    validate_args(args)
    
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
    # Handle argument parsing and cache formatting
    validate_i(args)
    initialise_param_cache(args)
    print("Initialisation complete!")

def dmain(args):
    # Create depth directory if it doesn't exist
    DEPTH_DIR = os.path.join(args.workingDirectory, "depth")
    os.makedirs(DEPTH_DIR, exist_ok=True)
    
    # Validate that necessary programs exist
    validate_samtools_exists()
    
    # Handle argument parsing and cache formatting
    paramsDict = load_param_cache(args.workingDirectory)
    validate_d(args, paramsDict)
    
    if args.bamFiles != []:
        print("# Updated BAM files listed in directory cache")
        update_param_cache(args.workingDirectory, {"bamFiles" : args.bamFiles})
    elif args.depthFiles != []:
        print("# Updated depth files listed in directory cache")
        update_param_cache(args.workingDirectory, {"depthFiles" : args.depthFiles})
    
    # Generate depth files
    if args.bamFiles != []: # if this list is not empty, we don't have existing depth files
        # Figure out depth file names
        depthIO = []
        for bamFile in args.bamFiles:
            # Derive the depth file name from the BAM file
            bamPrefix = os.path.basename(bamFile).rsplit(args.bamSuffix, maxsplit=1)[0]
            depthFile = os.path.join(DEPTH_DIR, f"{bamPrefix}{args.depthSuffix}")
            
            # Skip if the file already exists
            if os.path.isfile(depthFile) and os.path.isfile(depthFile + ".ok"):
                continue
            depthIO.append([bamFile, depthFile])
        
        # Run samtools depth
        run_samtools_depth(depthIO, args.threads)
        
        # Update the param cache
        print("# Updated depth files listed in directory cache")
        update_param_cache(args.workingDirectory, {"depthFiles" : [ f[1] for f in depthIO ]})
    else:
        print(f"# {'All' if len(args.depthFiles) > 1 else 'The'} {len(args.depthFiles)} " +
              f"depth file{'s' if len(args.depthFiles) > 1 else ''} " + 
              f"already exist{'s' if len(args.depthFiles) == 0 else ''}; skipping...")
    
    # Bin depth values according to window size
    ## TBD
    
    print("Depth file generation complete!")

def cmain(args):
    '''
    One note of caution: this script uses '.ok' files to track files that we know have
    been successfully processed. However, we do not apply the same for VCF indices and
    FASTA indices. As such, any unusual issues or bugs that occur may well be related to
    an index that was partially created. It is expected that stderr values will provide
    this information but this hasn't been validated yet.
    '''
    # Create call directory if it doesn't exist
    CALL_DIR = os.path.join(args.workingDirectory, "call")
    os.makedirs(CALL_DIR, exist_ok=True)
    
    # Validate that necessary programs exist
    validate_samtools_exists()
    validate_bcftools_exists()
    validate_bgzip_exists()
    validate_vt_exists()
    
    # Handle argument parsing and cache formatting
    paramsDict = load_param_cache(args.workingDirectory)
    validate_c(args, paramsDict)
    print("# Updated BAM files listed in directory cache")
    update_param_cache(args.workingDirectory, {"bamFiles" : args.bamFiles})
    
    # Index the reference genome (if necessary)
    if not os.path.isfile(args.genomeFasta + ".fai"):
        run_samtools_index(args.genomeFasta)
    
    # Create a bamlist file
    BAMLIST_FILE = os.path.join(CALL_DIR, "bamlist.txt")
    with open(BAMLIST_FILE, "w") as bamlistFile:
        for bamFile in args.bamFiles:
            bamlistFile.write(f"{bamFile}\n")
    
    # Run bcftools mpileup->call on each contig
    run_bcftools_call(BAMLIST_FILE, args.genomeFasta, CALL_DIR, args.threads)
    
    # Index each VCF file
    for vcfFile in [ os.path.join(CALL_DIR, f) for f in os.listdir(CALL_DIR) ]:
        if vcfFile.endswith(".vcf.gz") and not os.path.isfile(vcfFile + ".csi"):
            run_bcftools_index(vcfFile)
    
    # Run normalisation on each contig's VCF
    run_normalisation(args.genomeFasta, CALL_DIR, args.threads)
    
    # Index each VCF file
    for vcfFile in [ os.path.join(CALL_DIR, f) for f in os.listdir(CALL_DIR) ]:
        if vcfFile.endswith(".vcf.gz") and not os.path.isfile(vcfFile + ".csi"):
            run_bcftools_index(vcfFile)
    
    # Concatenate all VCF files
    FINAL_VCF_FILE = os.path.join(CALL_DIR, "psQTL_prep.vcf.gz")
    if not os.path.isfile(FINAL_VCF_FILE + ".ok"):
        run_bcftools_concat(args.genomeFasta, CALL_DIR, FINAL_VCF_FILE)
    
    # Index the final VCF file
    if not os.path.isfile(FINAL_VCF_FILE + ".csi"):
        run_bcftools_index(vcfFile)
    
    # Update the param and VCF caches
    update_param_cache(args.workingDirectory, {"vcfFile" : FINAL_VCF_FILE})
    initialise_vcf_cache(args.workingDirectory)
    
    print("Variant calling complete!")

def vmain(args):
    # Handle argument parsing and cache formatting
    paramsDict = load_param_cache(args.workingDirectory)
    
    # Present standard parameters
    print("# Parameters:")
    print(f"Working directory: {paramsDict['workingDirectory']}")
    print(f"Metadata file: {paramsDict['metadataFile']}")
    
    if paramsDict['bamFiles'] != []:
        print(f"BAM files: {paramsDict['bamFiles']}")
    else:
        print("BAM files: None")
    
    if paramsDict['depthFiles'] != []:
        print(f"Depth files: {paramsDict['depthFiles']}")
    else:
        print("Depth files: None")
    
    # Present VCF cache
    print() # blank line for spacing
    print("# VCF details:")
    if paramsDict['vcfFile'] is not None:
        vcfDict = load_vcf_cache(args.workingDirectory)
        if vcfDict == {}:
            initialise_vcf_cache(args.workingDirectory)
        
        print(f"VCF file: {paramsDict['vcfFile']}")
        print(f"Num. variants: {vcfDict['variants']}")
        print(f"{len(vcfDict['samples'])} samples: {vcfDict['samples']}")
        print(f"{len(vcfDict['contigs'])} contigs: {vcfDict['contigs']}")
    else:
        print("VCF files: None")
    print() # blank line for spacing
    
    print("Metadata viewing complete!")

if __name__ == "__main__":
    main()
