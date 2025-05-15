#! python3
# psQTL_proc.py
# Represents step 2 in the psQTL pipeline, where the user can take a VCF or VCF-like file and
# 'proc'ess it to calculate Euclidean distance values of variants or deletions among bulks.
# These results can then be used as input to psQTL_post.py to plot and tabulate relevant details.

import os, argparse, sys, gzip

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.parameters import ParameterCache
from modules.parsing import parse_metadata
from modules.ed import parse_vcf_for_ed
from _version import __version__

def validate_args(args):
    # Validate working directory
    args.workingDirectory = os.path.abspath(args.workingDirectory)
    if not os.path.isdir(args.workingDirectory):
        raise FileNotFoundError(f"-d working directory '{args.workingDirectory}' is not a directory!")
    
    # Validate cache existence & merge into args
    paramsCache = ParameterCache(args.workingDirectory)
    paramsCache.merge(args) # raises FileNotFoundError if cache does not exist
    
    # Validate metadata file
    if args.metadataFile == None:
        raise FileNotFoundError("Working directory has not been initialised with a metadata file!")
    elif not os.path.isfile(args.metadataFile):
            raise FileNotFoundError(f"Metadata file '{args.metadataFile}' was identified in " +
                                    "the parameters cache, but it doesn't exist or is not a file!")

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

def generate_ed_file(vcfFile, metadataDict, outputFileName, ignoreIdentical, parents):
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
    
    cparser = subparsers.add_parser("call",
                                    parents=[p],
                                    add_help=False,
                                    help="Analyse variant calls in the VCF file")
    cparser.set_defaults(func=cmain)
    
    dparser = subparsers.add_parser("depth",
                                    parents=[p],
                                    add_help=False,
                                    help="Analyse depth-predicted deletions in the VCF-like file")
    dparser.set_defaults(func=dmain)
    
    # Call-subparser arguments
    cparser.add_argument("--parents", dest="parents",
                         required=False,
                         nargs=2,
                         help="""Optionally, specify the names of the two parents (one from each phenotype
                         group) from which the samples used here are derived. This is used to determine
                         the bulk names in the metadata file. If not provided, psQTL will run
                         agnostically and just look for segregation
                         """,
                         default=[])
    cparser.add_argument("--ignoreIdentical", dest="ignoreIdentical",
                         required=False,
                         action="store_true",
                         help="""Optionally, provide this flag if you'd like variants where
                         both bulks are identical to be ignored; this can occur when both bulks
                         have the same variant with respect to the reference genome""",
                         default=False)
    
    # Depth-subparser arguments
    # N/A
    
    args = subParentParser.parse_args()
    validate_args(args)
    
    # Parse metadata file
    metadataDict = parse_metadata(args.metadataFile)
    
    # Validate that parent samples are in the metadata file
    if args.parents != []:
        b1Found = False
        b2Found = False
        for parent in args.parents:
            if parent in metadataDict["bulk1"]:
                b1Found = True
            elif parent in metadataDict["bulk2"]:
                b2Found = True
            else:
                raise ValueError(f"Parent '{parent}' not found in metadata file!")
        if not b1Found:
            raise ValueError("Neither parent was found in bulk1!")
        if not b2Found:
            raise ValueError("Neither parent was found in bulk2!")
    
    # Split into mode-specific functions
    if args.mode == "call":
        print("## psQTL_proc.py - Variant Distances ##")
        cmain(args, metadataDict)
    elif args.mode == "depth":
        print("## psQTL_proc.py - Depth Distances ##")
        dmain(args, metadataDict)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def cmain(args, metadataDict):
    validate_c(args)
    FINAL_ED_FILE = os.path.join(args.workingDirectory, "psQTL_variants.ed.tsv.gz")
    
    # Generate Euclidean distance file if it doesn't exist
    if not os.path.isfile(FINAL_ED_FILE + ".ok"):
        generate_ed_file(args.vcfFile, metadataDict, FINAL_ED_FILE, args.ignoreIdentical, args.parents)
    # If the .ok file exists, raise an error
    else:
        raise FileExistsError(f"Euclidean distance file '{FINAL_ED_FILE}' already has a .ok file; " +
                              "move, rename, or delete it before re-running psQTL_proc.py!")

def dmain(args, metadataDict):
    validate_d(args)
    FINAL_ED_FILE = os.path.join(args.workingDirectory, "psQTL_depth.ed.tsv.gz")
    
    # Generate Euclidean distance file if it doesn't exist
    if not os.path.isfile(FINAL_ED_FILE + ".ok"):
        generate_ed_file(args.deletionFile, metadataDict, FINAL_ED_FILE, False, args.parents) # don't ignore identical
    # If the .ok file exists, raise an error
    else:
        raise FileExistsError(f"Euclidean distance file '{FINAL_ED_FILE}' already has a .ok file; " +
                              "move, rename, or delete it before re-running psQTL_proc.py!")

if __name__ == "__main__":
    main()
