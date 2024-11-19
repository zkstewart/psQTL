#! python3
# psQTL_proc.py
# Represents step 2 in the psQTL pipeline, where the user can take a VCF or VCF-like file and
# 'proc'ess it to calculate Euclidean distance values of variants or deletions among bulks.
# These results can then be used as input to psQTL_post.py to plot and tabulate relevant details.

import os, argparse, sys, gzip

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.cache_handling import load_param_cache, load_vcf_cache, load_deletion_cache, \
                                   merge_cache_into_args, initialise_deletion_cache, initialise_vcf_cache
from modules.parsing import parse_metadata
from modules.ed import parse_vcf_for_ed

def validate_args(args):
    # Validate working directory
    args.workingDirectory = os.path.abspath(args.workingDirectory)
    if not os.path.exists(args.workingDirectory):
        raise FileNotFoundError(f"-d working directory '{args.workingDirectory}' does not exist!")
    
    # Validate cache existence
    paramsDict = load_param_cache(args.workingDirectory)
    if paramsDict == {}:
        raise FileNotFoundError("Working directory has not been initialised;" +
                                f"parameter cache not found in '{args.workingDirectory}'!")
    
    # Validate metadata file
    if paramsDict["metadataFile"] == None:
        raise FileNotFoundError("Working directory has not been initialised with a metadata file!")
    elif not os.path.isfile(paramsDict["metadataFile"]):
            raise FileNotFoundError(f"Metadata file '{paramsDict['metadataFile']}' was identified in " +
                                    "the parameters cache, but it's now absent!")

def validate_c(args):
    '''
    Params cache should have been merged into args before calling this function.
    '''
    if args.vcfFile == None:
        raise FileNotFoundError("Working directory has not been initialised with a VCF file!")
    else:
        args.vcfFile = os.path.abspath(args.vcfFile)
        if not os.path.isfile(args.vcfFile):
            raise FileNotFoundError(f"VCF file '{args.vcfFile}' was identified in " +
                                    "the parameters cache, but it's now absent!")
    
    # Validate VCF cache
    vcfDict = load_vcf_cache(args.workingDirectory)
    if vcfDict == {}:
        print("## VCF cache not found; re-initialising...")
        initialise_vcf_cache(args.workingDirectory)

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
                                    "the parameters cache, but it's now absent!")
    
    # Validate deletion cache
    deletionDict = load_deletion_cache(args.workingDirectory)
    if deletionDict == {}:
        print("## Deletion cache not found; re-initialising...")
        initialise_deletion_cache(args.workingDirectory, None)

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
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser(description=usage)
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
    merge_cache_into_args(args)
    
    # Parse metadata file
    metadataDict = parse_metadata(args.metadataFile)
    
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
        generate_ed_file(args.vcfFile, metadataDict, FINAL_ED_FILE, args.ignoreIdentical)
    # If the .ok file exists, raise an error
    else:
        raise FileExistsError(f"Euclidean distance file '{FINAL_ED_FILE}' already has a .ok file; " +
                              "move, rename, or delete it before re-running psQTL_proc.py!")

def dmain(args, metadataDict):
    validate_d(args)
    FINAL_ED_FILE = os.path.join(args.workingDirectory, "psQTL_depth.ed.tsv.gz")
    
    # Generate Euclidean distance file if it doesn't exist
    if not os.path.isfile(FINAL_ED_FILE + ".ok"):
        generate_ed_file(args.deletionFile, metadataDict, FINAL_ED_FILE, False) # don't ignore identical
    # If the .ok file exists, raise an error
    else:
        raise FileExistsError(f"Euclidean distance file '{FINAL_ED_FILE}' already has a .ok file; " +
                              "move, rename, or delete it before re-running psQTL_proc.py!")

if __name__ == "__main__":
    main()
