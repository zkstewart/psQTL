#! python3
# psQTL_proc.py
# Represents step 2 in the psQTL pipeline, where the user can take a VCF or VCF-like file and
# 'proc'ess it to calculate Euclidean distance values of variants or deletions among bulks.
# These results can then be used as input to psQTL_post.py to plot and tabulate relevant details.

import os, argparse, sys, gzip

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.validation import validate_proc_args, validate_c, validate_d
from modules.parsing import parse_metadata
from modules.ed import parse_vcf_for_ed
from modules.splsda import recode_vcf
from _version import __version__

def generate_ed_file(vcfFile, metadataDict, outputFileName, ignoreIdentical):
    '''
    Parameters:
        vcfFile -- a string indicating the path to the VCF file to be processed
        metadataDict -- a dictionary of metadata information parsed from the metadata file
        outputFileName -- a string indicating the path to the output file to be written
        ignoreIdentical -- a boolean indicating whether to ignore variants where both bulks
                           are identical
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
        euclideanDist in parse_vcf_for_ed(vcfFile, metadataDict, ignoreIdentical):
            # Write content line
            fileOut.write(f"{contig}\t{pos}\t{variant}\t{numAllelesB1}\t" + \
                            f"{numAllelesB2}\t{euclideanDist}\n")
    open(outputFileName + ".ok", "w").close() # touch a .ok file to indicate success

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
    eparser.add_argument("--ignoreIdentical", dest="ignoreIdentical",
                         required=False,
                         action="store_true",
                         help="""Optionally, provide this flag if you'd like variants where
                         both bulks are identical to be ignored; this can occur when both bulks
                         have the same variant with respect to the reference genome""",
                         default=False)
    
    # sPLS-DA-subparser arguments
    ## N/A
    
    args = subParentParser.parse_args()
    locations = validate_proc_args(args)
    
    # Parse metadata file
    metadataDict = parse_metadata(args.metadataFile)
    
    # Split into mode-specific functions
    if args.mode == "ed":
        print("## psQTL_proc.py - Euclidean Distances ##")
        emain(args, metadataDict, locations)
    elif args.mode == "depth":
        print("## psQTL_proc.py - sPLS-DA ##")
        dmain(args, metadataDict, locations)
    
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
        generate_ed_file(args.vcfFile, metadataDict, locations.variantEdFile, args.ignoreIdentical)
    else:
        raise FileExistsError(f"Euclidean distance file '{locations.variantEdFile}' already has a .ok file; " +
                              "move, rename, or delete it before re-running psQTL_proc.py!")

def depth_ed(args, metadataDict, locations):
    if not os.path.isfile(locations.depthEdFile + ".ok"):
        generate_ed_file(args.deletionFile, metadataDict, locations.depthEdFile, False) # don't ignore identical
    else:
        raise FileExistsError(f"Euclidean distance file '{locations.depthEdFile}' already has a .ok file; " +
                              "move, rename, or delete it before re-running psQTL_proc.py!")

def smain(args, metadataDict, locations):
    # Validate input types before proceeding
    if "call" in args.inputType:
        validate_c(args)
    if "depth" in args.inputType:
        validate_d(args)
    
    # Run the main function for each input type
    if "call" in args.inputType:
        call_splsda(args, metadataDict, locations)
    if "depth" in args.inputType:
        depth_splsda(args, metadataDict, locations)
    if "call" in args.inputType and "depth" in args.inputType:
        integrative_splsda(args, metadataDict, locations)

def call_splsda(args, metadataDict, locations):
    raise NotImplementedError("sPLS-DA for call variants is not yet implemented.")

def depth_splsda(args, metadataDict, locations):
    raise NotImplementedError("sPLS-DA for depth variants is not yet implemented.")

def integrative_splsda(args, metadataDict, locations):
    raise NotImplementedError("Integrative sPLS-DA for both call and depth variants is not yet implemented.")

if __name__ == "__main__":
    main()
