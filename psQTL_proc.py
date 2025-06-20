#! python3
# psQTL_proc.py
# Represents step 2 in the psQTL pipeline, where the user can take a VCF or VCF-like file and
# 'proc'ess it to calculate Euclidean distance values of variants or CNVs among groups.
# These results can then be used as input to psQTL_post.py to plot and tabulate relevant details.

import os, argparse, sys, gzip

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.validation import validate_proc_args, validate_c, validate_d, validate_s
from modules.parsing import parse_metadata, WriteGzFile
from modules.ed import parse_vcf_for_ed
from modules.splsda import validate_r_exists, validate_r_packages_installation, \
    recode_vcf, run_windowed_splsda, run_integrative_splsda
from _version import __version__

def generate_ed_files(vcfFile, metadataDict, edOutputFiles,
                      parentSamples=[], isCNV=False, ignoreIdentical=True):
    '''
    Note that parse_vcf_for_ed yields a dictionary with the following keys
    contig, pos, variant, numAllelesG1, numAllelesG2, numFilteredG1, numFilteredG2,
    alleleED, genotypeED, inheritanceED, ploidy, possibleAllelesG1, possibleAllelesG2
    
    Parameters:
        vcfFile -- a string indicating the path to the VCF file to be processed
        metadataDict -- a dictionary of metadata information parsed from the metadata file
        edOutputFiles -- a list of strings indicating the paths to the output files
                         with None used to indicate no ED file for that type;
                         order as follows: allele ED, genotype ED, inheritance ED
        parentSamples -- (OPTIONAL) a list of two strings indicating the names of the
                         two parents used to generate the groups for use in calculating
                         inheritance ED; OR an empty list if no parents are used
        isCNV -- (OPTIONAL) a boolean indicating whether the VCF file contains CNV data
                 (default: False, meaning the VCF file contains variant calls)
        ignoreIdentical -- (OPTIONAL) a boolean indicating whether to ignore variants
                           where both groups are identical
    '''
    with WriteGzFile(edOutputFiles[0]) as allelesOut, WriteGzFile(edOutputFiles[1]) as genotypesOut, \
    WriteGzFile(edOutputFiles[2]) as inheritanceOut:
        contexts = [ allelesOut, genotypesOut, inheritanceOut ] # some contexts may be None
        
        # Write header lines for each output file
        for context in contexts:
            if context != None:
                context.write("{0}\n".format("\t".join([
                    "CHROM", "POSI", "variant",
                    "group1_alleles", "group1_alleles_possible",
                    "group2_alleles", "group2_alleles_possible",
                    "euclideanDist"
                ])))
        
        # Iterate through Euclidean distance calculations for VCF file
        for resultsDict in parse_vcf_for_ed(vcfFile, metadataDict, isCNV,
                                            parents=parentSamples,
                                            ignoreIdentical=ignoreIdentical):
            # Write allele frequency ED line
            if allelesOut != None:
                allelesOut.write(f"{resultsDict['contig']}\t{resultsDict['pos']}\t{resultsDict['variant']}\t" + \
                                 f"{resultsDict['numAllelesG1']}\t{resultsDict['possibleAllelesG1']}\t" + \
                                 f"{resultsDict['numAllelesG2']}\t{resultsDict['possibleAllelesG2']}\t" + \
                                 f"{resultsDict['alleleED']}\n")
            # Write genotype frequency ED line
            if genotypesOut != None:
                genotypesOut.write(f"{resultsDict['contig']}\t{resultsDict['pos']}\t{resultsDict['variant']}\t" + \
                                   f"{resultsDict['numAllelesG1']}\t{resultsDict['possibleAllelesG1']}\t" + \
                                   f"{resultsDict['numAllelesG2']}\t{resultsDict['possibleAllelesG2']}\t" + \
                                   f"{resultsDict['genotypeED']}\n")
            # Write inheritance ED line
            if inheritanceOut != None:
                inheritanceOut.write(f"{resultsDict['contig']}\t{resultsDict['pos']}\t{resultsDict['variant']}\t" + \
                                     f"{resultsDict['numFilteredG1']}\t{resultsDict['possibleAllelesG1']}\t" + \
                                     f"{resultsDict['numFilteredG2']}\t{resultsDict['possibleAllelesG2']}\t" + \
                                     f"{resultsDict['inheritanceED']}\n")

def main():
    usage = """%(prog)s processes VCF or VCF-like files containing variant or CNV
    predictions for samples belonging to two groups. The working directory is expected
    to have been 'initialise'd by psQTL_prep.py such that all necessary files have been
    validated and prepared for analysis.
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
                                    help="Compute Euclidean distance of 'call' and/or 'depth' variants")
    eparser.set_defaults(func=emain)
    
    sparser = subparsers.add_parser("splsda",
                                    parents=[p],
                                    add_help=False,
                                    help="Compute local sPLS-DA of 'call' and/or 'depth' variants")
    sparser.set_defaults(func=smain)
    
    # ED-subparser arguments
    eparser.add_argument("--parents", dest="parentSamples",
                         required=False,
                         nargs=2,
                         help="""Optionally, provide the names of the two parents used to
                         generate the groups; this is used to apply an alternative form of ED
                         which leverages the parents' genotypes to extract more signal out of
                         your data. If not provided, the standard ED will be used.""",
                         default=[])
    eparser.add_argument("--considerIdentical", dest="considerIdentical",
                         required=False,
                         action="store_true",
                         help="""Optionally, provide this flag to prevent filtration of
                         variants where both groups' genotypes are identical; this can
                         occur when both groups have the same variant with respect to
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
    sparser.add_argument("--splsdaWindowSize", dest="splsdaWindowSize",
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
    
    # Validate based on input types
    if "call" in args.inputType:
        validate_c(args)
    if "depth" in args.inputType:
        validate_d(args)
    
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
    # Run the main function for each input type
    if "call" in args.inputType:
        call_ed(args, metadataDict, locations)
    if "depth" in args.inputType:
        depth_ed(args, metadataDict, locations)

def call_ed(args, metadataDict, locations):
    # Identify the names for variant ED files
    edOutputFiles = [locations.allelesEdFile, locations.genotypesEdFile,
                     locations.inheritanceEdFile if args.parentSamples != [] else None]
    
    # Generate ED files for variant calls
    if not all([ os.path.isfile(edFile + ".ok") for edFile in edOutputFiles if edFile != None ]) or \
    not all([ os.path.isfile(edFile) for edFile in edOutputFiles if edFile != None ]):
        generate_ed_files(args.vcfFile, metadataDict, edOutputFiles,
                          parentSamples=args.parentSamples,
                          isCNV=False,
                          ignoreIdentical=not args.considerIdentical) # negate the flag to ignore identical
        for edFile in edOutputFiles:
            if edFile != None:
                open(edFile + ".ok", "w").close() # touch a .ok file to indicate success
    else:
        formattedFiles = " and ".join([f"'{edFile}'" for edFile in edOutputFiles if edFile is not None])
        
        raise FileExistsError(f"Euclidean distance files {formattedFiles} already have .ok files; " +
                              "move, rename, or delete these before re-running psQTL_proc.py!")
    print("Variant call ED file generation complete!")

def depth_ed(args, metadataDict, locations):
    # Format names for CNV ED files
    edOutputFiles = [locations.depthEdFile, None, None] # no alternatives for CNVs
    
    # Generate ED files for CNV variants
    if not all([ os.path.isfile(edFile + ".ok") for edFile in edOutputFiles if edFile != None ]) or \
    not all([ os.path.isfile(edFile) for edFile in edOutputFiles if edFile != None ]):
        generate_ed_files(args.depthFile, metadataDict, edOutputFiles,
                          parentSamples=[], # no parents used for CNVs
                          isCNV=True,
                          ignoreIdentical=False) # don't ignore identical
        for edFile in edOutputFiles:
            if edFile != None:
                open(edFile + ".ok", "w").close() # touch a .ok file to indicate success
    else:
        raise FileExistsError(f"Euclidean distance file '{locations.depthEdFile}' already has a .ok file; " +
                              "move, rename, or delete it before re-running psQTL_proc.py!")
    print("CNV variant ED file generation complete!")

def smain(args, metadataDict, locations):
    # Validate sPLS-DA arguments
    validate_s(args)
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
        print("# Encoding 'call' variants for sPLS-DA analysis ...")
        recode_vcf(args.vcfFile, locations.variantRecodedFile, metadataDict, isCNV=False)
        open(locations.variantRecodedFile + ".ok", "w").close() # touch a .ok file to indicate success
    else:
        print("# 'call' variants already encoded for sPLS-DA analysis; skipping ...")
    
    # Run windowed sPLS-DA for variant calls
    if (not os.path.isfile(locations.variantSplsdaSelectedFile) or \
        not os.path.isfile(locations.variantSplsdaSelectedFile + ".ok")) or \
        (not os.path.isfile(locations.variantSplsdaBerFile) or \
        not os.path.isfile(locations.variantSplsdaBerFile + ".ok")) or \
        (not os.path.isfile(locations.variantSplsdaRdataFile) or \
        not os.path.isfile(locations.variantSplsdaRdataFile + ".ok")):
            print("# Running windowed sPLS-DA for 'call' variants ...")
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
            print("'call' sPLS-DA analysis complete!")
    else:
        print("# 'call' variants already processed for sPLS-DA analysis; skipping ...")

def depth_splsda(args, metadataDict, locations):
    # Encode CNVs for sPLS-DA analysis
    if (not os.path.isfile(locations.depthRecodedFile)) or (not os.path.isfile(locations.depthRecodedFile + ".ok")):
        print("# Encoding 'depth' CNVs for sPLS-DA analysis ...")
        recode_vcf(args.depthFile, locations.depthRecodedFile, metadataDict, isCNV=True)
        open(locations.depthRecodedFile + ".ok", "w").close() # touch a .ok file to indicate success
    else:
        print("# 'depth' CNVs already encoded for sPLS-DA analysis; skipping ...")
    
    # Run windowed sPLS-DA for CNV variants
    if (not os.path.isfile(locations.depthSplsdaSelectedFile) or \
        not os.path.isfile(locations.depthSplsdaSelectedFile + ".ok")) or \
        (not os.path.isfile(locations.depthSplsdaBerFile) or \
        not os.path.isfile(locations.depthSplsdaBerFile + ".ok")) or \
        (not os.path.isfile(locations.depthSplsdaRdataFile) or \
        not os.path.isfile(locations.depthSplsdaRdataFile + ".ok")):
            print("# Running windowed sPLS-DA for 'depth' CNVs ...")
            run_windowed_splsda(args.metadataFile, locations.depthRecodedFile,
                                locations.depthSplsdaSelectedFile,
                                locations.depthSplsdaBerFile,
                                locations.depthSplsdaRdataFile,
                                locations.windowedSplsdaRscript,
                                args.threads, args.splsdaWindowSize, args.berFilter,
                                args.mafFilter, args.numRepeats, args.maxIterations)
            open(locations.depthSplsdaSelectedFile + ".ok", "w").close()
            open(locations.depthSplsdaBerFile + ".ok", "w").close()
            open(locations.depthSplsdaRdataFile + ".ok", "w").close()
            print("'depth' sPLS-DA analysis complete!")
    else:
        print("# 'depth' CNVs already processed for sPLS-DA analysis; skipping ...")

def integrative_splsda(args, metadataDict, locations):
    # Check if it is possible to run integrative sPLS-DA
    if (not os.path.isfile(locations.variantSplsdaRdataFile) or \
        not os.path.isfile(locations.depthSplsdaRdataFile)):
        raise FileNotFoundError(f"Cannot run integrative sPLS-DA without both 'call' ({locations.variantSplsdaRdataFile}) " + 
                                f"and 'depth' ({locations.depthSplsdaRdataFile}) sPLS-DA RData files!")
    
    # Run integrative sPLS-DA for variant calls
    if (not os.path.isfile(locations.integrativeSplsdaSelectedFile) or \
        not os.path.isfile(locations.integrativeSplsdaSelectedFile + ".ok")):
            print("# Running integration of sPLS-DA for 'call' and 'depth' variants ...")
            run_integrative_splsda(locations.variantSplsdaRdataFile, locations.depthSplsdaRdataFile,
                                   locations.integrativeSplsdaSelectedFile,
                                   locations.integrativeSplsdaRscript,
                                   args.threads, args.numRepeats, args.maxIterations)
            open(locations.integrativeSplsdaSelectedFile + ".ok", "w").close()
            print("Integrative sPLS-DA analysis complete!")
    else:
        print("# Integrative sPLS-DA of variants already processed; skipping ...")

if __name__ == "__main__":
    main()
