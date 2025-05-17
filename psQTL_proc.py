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
from modules.splsda import validate_r_exists, validate_r_packages_installation, \
    recode_vcf, run_windowed_splsda
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
        open(locations.variantEdFile + ".ok", "w").close() # touch a .ok file to indicate success
    else:
        raise FileExistsError(f"Euclidean distance file '{locations.variantEdFile}' already has a .ok file; " +
                              "move, rename, or delete it before re-running psQTL_proc.py!")
    print("Variant call ED file generation complete!")

def depth_ed(args, metadataDict, locations):
    if not os.path.isfile(locations.deletionEdFile + ".ok"):
        generate_ed_file(args.deletionFile, metadataDict, locations.deletionEdFile, False) # don't ignore identical
        open(locations.deletionEdFile + ".ok", "w").close() # touch a .ok file to indicate success
    else:
        raise FileExistsError(f"Euclidean distance file '{locations.deletionEdFile}' already has a .ok file; " +
                              "move, rename, or delete it before re-running psQTL_proc.py!")
    print("Deletion variant ED file generation complete!")

def smain(args, metadataDict, locations):
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
    
    # Run windowed sPLS-DA for variant calls
    if (not os.path.isfile(locations.variantSplsdaSelectedFile) or \
        not os.path.isfile(locations.variantSplsdaSelectedFile + ".ok")) or \
        (not os.path.isfile(locations.variantSplsdaBerFile) or \
        not os.path.isfile(locations.variantSplsdaBerFile + ".ok")):
            print("# Running windowed sPLS-DA for variant calls ...")
            run_windowed_splsda(args.metadataFile, locations.variantRecodedFile,
                                locations.variantSplsdaSelectedFile,
                                locations.variantSplsdaBerFile,
                                locations.windowedSplsdaRscript)
    print("Variant call sPLS-DA analysis complete!")

def depth_splsda(args, metadataDict, locations):
    # Encode deletion variants for sPLS-DA analysis
    if (not os.path.isfile(locations.variantRecodedFile)) or (not os.path.isfile(locations.variantRecodedFile + ".ok")):
        print("# Encoding deletion variants for sPLS-DA analysis ...")
        recode_vcf(args.deletionFile, locations.deletionRecodedFile)
    
    # Run windowed sPLS-DA for deletion variants
    if (not os.path.isfile(locations.deletionSplsdaSelectedFile) or \
        not os.path.isfile(locations.deletionSplsdaSelectedFile + ".ok")) or \
        (not os.path.isfile(locations.deletionSplsdaBerFile) or \
        not os.path.isfile(locations.deletionSplsdaBerFile + ".ok")):
            print("# Running windowed sPLS-DA for deletion variants ...")
            run_windowed_splsda(args.metadataFile, locations.deletionRecodedFile,
                                locations.deletionSplsdaSelectedFile,
                                locations.deletionSplsdaBerFile,
                                locations.windowedSplsdaRscript)
    print("Deletion variant sPLS-DA analysis complete!")

def integrative_splsda(args, metadataDict, locations):
    raise NotImplementedError("Integrative sPLS-DA for both call and depth variants is not yet implemented.")

if __name__ == "__main__":
    main()
