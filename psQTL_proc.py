#! python3
# psQTL_proc.py
# Represents step 2 in the psQTL pipeline, where the user can take a VCF or VCF-like file and
# 'proc'ess it to calculate Euclidean distance values of variants or deletions among bulks.
# These results can then be used as input to psQTL_post.py to plot and tabulate relevant details.

import os, argparse, sys, gzip

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.validation import validate_proc_args, validate_c, validate_d, validate_s
from modules.parsing import parse_metadata
from modules.ed import parse_vcf_for_ed
from modules.splsda import validate_r_exists, validate_r_packages_installation, \
    recode_vcf, run_windowed_splsda, run_integrative_splsda
from _version import __version__

def generate_ed_file(vcfFile, metadataDict, outputFileName, parentSamples=[], isCNV=False, ignoreIdentical=True):
    '''
    Parameters:
        vcfFile -- a string indicating the path to the VCF file to be processed
        metadataDict -- a dictionary of metadata information parsed from the metadata file
        outputFileName -- a string indicating the path to the output file to be written
        isCNV -- (OPTIONAL) a boolean indicating whether the VCF file contains CNV data
                 (default: False, meaning the VCF file contains variant calls)
        ignoreIdentical -- (OPTIONAL) a boolean indicating whether to ignore variants
                           where both bulks are identical
    '''
    with gzip.open(outputFileName, "wt") as fileOut:
        # Write header line
        fileOut.write("{0}\n".format("\t".join([
            "CHROM", "POSI", "variant",
            "bulk1_alleles", "bulk2_alleles", 
            "euclideanDist"
        ])))
        
        # Iterate through Euclidean distance calculations for VCF file
        for contig, pos, variant, numAllelesB1, numAllelesB2, \
        euclideanDist in parse_vcf_for_ed(vcfFile, metadataDict, isCNV,
                                          parents=parentSamples,
                                          ignoreIdentical=ignoreIdentical):
            # Write content line
            fileOut.write(f"{contig}\t{pos}\t{variant}\t{numAllelesB1}\t" + \
                            f"{numAllelesB2}\t{euclideanDist}\n")

def main():
    usage = """%(prog)s processes VCF or VCF-like files to calculate Euclidean distance values
    of variants or deletions among bulks. The input directory is expected to have been
    'initialise'd by psQTL_prep.py such that all necessary files have been validated and
    prepared for analysis.
    """
    # Establish parser
    p = argparse.ArgumentParser()
    
    # Set arguments shared by subparsers
    p.add_argument("-d", dest="workingDirectory",
                    required=True,
                    help="Specify the location where the analysis is being performed")
    p.add_argument("-i", dest="inputType",
                   required=True,
                   nargs="+",
                   choices=["call", "depth"],
                   help="""Specify one or both of 'call' and 'depth' to indicate which
                   types of variants to process.""")
    p.add_argument("-v", "--version",
                   action="version",
                   version="psQTL_proc.py {version}".format(version=__version__))
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser(description=usage)
    subParentParser.add_argument("-v", "--version",
                                 action="version",
                                 version="psQTL_proc.py {version}".format(version=__version__))
    
    subparsers = subParentParser.add_subparsers(dest="mode",
                                                required=True)
    
    eparser = subparsers.add_parser("ed",
                                    parents=[p],
                                    add_help=False,
                                    help="Compute Euclidean distance of call and/or depth variants")
    eparser.set_defaults(func=emain)
    
    sparser = subparsers.add_parser("splsda",
                                    parents=[p],
                                    add_help=False,
                                    help="Compute local sPLS-DA of call and/or depth variants")
    sparser.set_defaults(func=smain)
    
    # ED-subparser arguments
    eparser.add_argument("--parents", dest="parentSamples",
                         required=False,
                         nargs=2,
                         help="""Optionally, provide the names of the two parents used to
                         generate the bulks; this is used to apply an alternative form of ED
                         which leverages the parents' genotypes to extract more signal out of
                         your data. If not provided, the standard ED will be used.""",
                         default=[])
    eparser.add_argument("--considerIdentical", dest="considerIdentical",
                         required=False,
                         action="store_true",
                         help="""Optionally, provide this flag to prevent filtration of
                         variants where both bulks' genotypes are identical; this can
                         occur when both bulks have the same variant with respect to
                         the reference genome. Not recommended unless you have a
                         specific reason to do so.""",
                         default=False)
    
    # sPLS-DA-subparser arguments
    sparser.add_argument("--threads", dest="threads",
                         type=int,
                         required=False,
                         help="""Optionally, specify the number of threads to use when running
                         sPLS-DA analyses (default: 1)""",
                         default=1)
    sparser.add_argument("--windowSize", dest="splsdaWindowSize",
                         type=int,
                         required=False,
                         help="""Optionally, specify the window size that sPLS-DA will
                         be run within (recommended: 100000)""",
                         default=100000)
    sparser.add_argument("--maf", dest="mafFilter",
                         type=float,
                         required=False,
                         help="""Optionally, specify the Minor Allele Frequency to filter
                         variants by prior to sPLS-DA analysis (recommended: 0.05;
                         equivalent to 5 percent MAF)""",
                         default=0.05)
    sparser.add_argument("--ber", dest="berFilter",
                         type=float,
                         required=False,
                         help="""Optionally, specify the Balanced Error Rate to select
                         variants during sPLS-DA analysis (recommended: 0.4)""",
                         default=0.4)
    sparser.add_argument("--maxiters", dest="maxIterations",
                         type=int,
                         required=False,
                         help="""Optionally, specify the maximum number of iterations
                         allowed for sPLS-DA convergence (recommended: 10000)""",
                         default=10000)
    sparser.add_argument("--nrepeat", dest="numRepeats",
                         type=int,
                         required=False,
                         help="""Optionally, specify the number of times to repeat
                         M-fold validation during sPLS-DA optimisation (recommended: 20
                         (default) or more if you have the time)""",
                         default=20)
    
    args = subParentParser.parse_args()
    locations = validate_proc_args(args)
    
    # Parse metadata file
    metadataDict = parse_metadata(args.metadataFile)
    
    # Split into mode-specific functions
    if args.mode == "ed":
        print("## psQTL_proc.py - Euclidean Distances ##")
        emain(args, metadataDict, locations)
    elif args.mode == "splsda":
        print("## psQTL_proc.py - sPLS-DA ##")
        smain(args, metadataDict, locations)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def emain(args, metadataDict, locations):
    # Validate input types before proceeding
    if "call" in args.inputType:
        validate_c(args)
    if "depth" in args.inputType:
        validate_d(args)
    
    # Run the main function for each input type
    if "call" in args.inputType:
        call_ed(args, metadataDict, locations)
    if "depth" in args.inputType:
        depth_ed(args, metadataDict, locations)

def call_ed(args, metadataDict, locations):
    if not os.path.isfile(locations.variantEdFile + ".ok"):
        generate_ed_file(args.vcfFile, metadataDict, locations.variantEdFile,
                         parentSamples=args.parentSamples,
                         isCNV=False,
                         ignoreIdentical=not args.considerIdentical) # negate the flag to ignore identical
        open(locations.variantEdFile + ".ok", "w").close() # touch a .ok file to indicate success
    else:
        raise FileExistsError(f"Euclidean distance file '{locations.variantEdFile}' already has a .ok file; " +
                              "move, rename, or delete it before re-running psQTL_proc.py!")
    print("Variant call ED file generation complete!")

def depth_ed(args, metadataDict, locations):
    if not os.path.isfile(locations.deletionEdFile + ".ok"):
        generate_ed_file(args.deletionFile, metadataDict, locations.deletionEdFile,
                         parents=None, # no parents used for CNVs
                         isCNV=True,
                         ignoreIdentical=False) # don't ignore identical
        open(locations.deletionEdFile + ".ok", "w").close() # touch a .ok file to indicate success
    else:
        raise FileExistsError(f"Euclidean distance file '{locations.deletionEdFile}' already has a .ok file; " +
                              "move, rename, or delete it before re-running psQTL_proc.py!")
    print("Deletion variant ED file generation complete!")

def smain(args, metadataDict, locations):
    # Validate sPLS-DA arguments
    validate_s(args)
    
    # Validate input types before proceeding
    if "call" in args.inputType:
        validate_c(args)
    if "depth" in args.inputType:
        validate_d(args)
    os.makedirs(locations.splsdaDir, exist_ok=True)
    
    # Validate that R and necessary packages are available
    validate_r_exists()
    validate_r_packages_installation()
    
    # Run the main function for each input type
    if "call" in args.inputType:
        call_splsda(args, metadataDict, locations)
    if "depth" in args.inputType:
        depth_splsda(args, metadataDict, locations)
    if "call" in args.inputType and "depth" in args.inputType:
        integrative_splsda(args, metadataDict, locations)

def call_splsda(args, metadataDict, locations):
    # Encode variant calls for sPLS-DA analysis
    if (not os.path.isfile(locations.variantRecodedFile)) or (not os.path.isfile(locations.variantRecodedFile + ".ok")):
        print("# Encoding variant calls for sPLS-DA analysis ...")
        recode_vcf(args.vcfFile, locations.variantRecodedFile)
        open(locations.variantRecodedFile + ".ok", "w").close() # touch a .ok file to indicate success
    else:
        print("# Variant calls already encoded for sPLS-DA analysis; skipping ...")
    
    # Run windowed sPLS-DA for variant calls
    if (not os.path.isfile(locations.variantSplsdaSelectedFile) or \
        not os.path.isfile(locations.variantSplsdaSelectedFile + ".ok")) or \
        (not os.path.isfile(locations.variantSplsdaBerFile) or \
        not os.path.isfile(locations.variantSplsdaBerFile + ".ok")) or \
        (not os.path.isfile(locations.variantSplsdaRdataFile) or \
        not os.path.isfile(locations.variantSplsdaRdataFile + ".ok")):
            print("# Running windowed sPLS-DA for variant calls ...")
            run_windowed_splsda(args.metadataFile, locations.variantRecodedFile,
                                locations.variantSplsdaSelectedFile,
                                locations.variantSplsdaBerFile,
                                locations.variantSplsdaRdataFile,
                                locations.windowedSplsdaRscript,
                                args.threads, args.splsdaWindowSize, args.berFilter,
                                args.mafFilter, args.numRepeats, args.maxIterations)
            open(locations.variantSplsdaSelectedFile + ".ok", "w").close()
            open(locations.variantSplsdaBerFile + ".ok", "w").close()
            open(locations.variantSplsdaRdataFile + ".ok", "w").close()
            print("Variant call sPLS-DA analysis complete!")
    else:
        print("# Variant calls already processed for sPLS-DA analysis; skipping ...")

def depth_splsda(args, metadataDict, locations):
    # Encode deletion variants for sPLS-DA analysis
    if (not os.path.isfile(locations.deletionRecodedFile)) or (not os.path.isfile(locations.deletionRecodedFile + ".ok")):
        print("# Encoding deletion variants for sPLS-DA analysis ...")
        recode_vcf(args.deletionFile, locations.deletionRecodedFile)
        open(locations.deletionRecodedFile + ".ok", "w").close() # touch a .ok file to indicate success
    else:
        print("# Deletion variants already encoded for sPLS-DA analysis; skipping ...")
    
    # Run windowed sPLS-DA for deletion variants
    if (not os.path.isfile(locations.deletionSplsdaSelectedFile) or \
        not os.path.isfile(locations.deletionSplsdaSelectedFile + ".ok")) or \
        (not os.path.isfile(locations.deletionSplsdaBerFile) or \
        not os.path.isfile(locations.deletionSplsdaBerFile + ".ok")) or \
        (not os.path.isfile(locations.deletionSplsdaRdataFile) or \
        not os.path.isfile(locations.deletionSplsdaRdataFile + ".ok")):
            print("# Running windowed sPLS-DA for deletion variants ...")
            run_windowed_splsda(args.metadataFile, locations.deletionRecodedFile,
                                locations.deletionSplsdaSelectedFile,
                                locations.deletionSplsdaBerFile,
                                locations.deletionSplsdaRdataFile,
                                locations.windowedSplsdaRscript,
                                args.threads, args.splsdaWindowSize, args.berFilter,
                                args.mafFilter, args.numRepeats, args.maxIterations)
            open(locations.deletionSplsdaSelectedFile + ".ok", "w").close()
            open(locations.deletionSplsdaBerFile + ".ok", "w").close()
            open(locations.deletionSplsdaRdataFile + ".ok", "w").close()
            print("Deletion variant sPLS-DA analysis complete!")
    else:
        print("# Deletion variants already processed for sPLS-DA analysis; skipping ...")

def integrative_splsda(args, metadataDict, locations):
    # Run integrative sPLS-DA for variant calls
    if (not os.path.isfile(locations.integrativeSplsdaSelectedFile) or \
        not os.path.isfile(locations.integrativeSplsdaSelectedFile + ".ok")):
            print("# Running integration of sPLS-DA for variants and deletions ...")
            run_integrative_splsda(locations.variantSplsdaRdataFile, locations.deletionSplsdaRdataFile,
                                   locations.integrativeSplsdaSelectedFile,
                                   locations.integrativeSplsdaRscript,
                                   args.threads, args.numRepeats, args.maxIterations)
            open(locations.integrativeSplsdaSelectedFile + ".ok", "w").close()
            print("Integrative sPLS-DA analysis complete!")
    else:
        print("# Integrative sPLS-DA of variants and deletions already processed; skipping ...")

if __name__ == "__main__":
    main()
