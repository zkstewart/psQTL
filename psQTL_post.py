#! python3
# psQTL_post.py
# Represents step 3 of the psQTL pipeline, which is to 'post'-process the data
# generated from psQTL_proc.py. It can plot segregation statistics with a plug-and-play
# approach of several different plot types (line, scatter, histogram, genes) or it
# can report on genes that are proximal to or contained within potential QTLs.

import os, argparse, re, sys, pickle
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.parameters import ParameterCache
from modules.parsing import parse_metadata
from modules.depth import parse_bins_as_dict, normalise_coverage_dict, convert_dict_to_depthncls
from modules.ed import parse_ed_as_dict, convert_dict_to_edncls
from modules.plotting import linescatter, histogram, genes, coverage, scalebar, NUM_SAMPLE_LINES
from modules.reporting import report_genes, report_depth
from modules.gff3 import GFF3
from _version import __version__

def validate_args(args):
    # Validate working directory
    if not os.path.exists(args.workingDirectory):
        raise FileNotFoundError(f"-d working directory '{args.workingDirectory}' does not exist!")
    
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
    
    # Validate input file
    if args.inputType == "call":
        args.inputFile = os.path.join(args.workingDirectory, "psQTL_variants.ed.tsv.gz")
    else:
        args.inputFile = os.path.join(args.workingDirectory, "psQTL_depth.ed.tsv.gz")
    
    if not os.path.isfile(args.inputFile):
        raise FileNotFoundError(f"Euclidean distance file '{args.inputFile}' does not exist!")
    elif not os.path.isfile(args.inputFile + ".ok"):
        raise FileNotFoundError(f"Euclidean distance file '{args.inputFile}' does not have a '.ok' flag!")
    
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
            args.gff3Obj = GFF3(args.annotationGFF3) # parsing now to raise errors early
            args.gff3Obj.create_ncls_index("gene")
    
    # Validate output file
    if os.path.exists(args.outputFileName):
        raise FileExistsError(f"-o output file '{args.outputFileName}' already exists!")
    args.outputFileName = os.path.abspath(args.outputFileName)
    
    if not os.path.isdir(os.path.dirname(args.outputFileName)):
        raise FileNotFoundError(f"-o parent directory '{os.path.dirname(args.outputFileName)}' does not exist!")

def validate_regions(args, lengthsDict):
    # Validate regions
    regions = []
    regionsRegex = re.compile(r"^([^:]+):(\d+)-(\d+)$")
    for region in args.regions:
        reMatch = regionsRegex.match(region)
        
        # Handle chr:start-end format
        if reMatch != None:
            contigID, start, end = reMatch.groups()
            start = int(start)
            end = int(end)
            
            # Validate contig ID
            if not contigID in lengthsDict:
                raise ValueError(f"--region contig ID '{contigID}' not found in the -f FASTA!")
            # Validate start and end positions
            if start < 0:
                raise ValueError(f"--region start position '{start}' is < 0!")
            reverse = False
            if start >= end:
                start, end = end, start
                reverse = True
            # Store region
            regions.append([contigID, start, end, reverse])
        # Handle invalid format
        elif ":" in region:
            raise ValueError(f"Invalid region input '{region}'; you included a ':' but did " + 
                             "not format the region as 'chr:start-end'!")
        # Handle chr format
        else:
            if not region in lengthsDict:
                raise ValueError(f"--region contig ID '{region}' not found in the -f FASTA!")
            regions.append([region, 0, lengthsDict[region], False])
    
    # Handle empty regions
    if regions == []:
        regions = [[contigID, 0, lengthsDict[contigID], False] for contigID in lengthsDict]
    
    args.regions = regions

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
    if args.binSize < 2:
        raise ValueError(f"--bin value '{args.binSize}' must be >= 2!")
    if args.binThreshold < 0:
        raise ValueError(f"--threshold value '{args.binThreshold}' must be >= 0!")
    
    # Validate plot types
    if len(set(args.plotTypes)) != len(args.plotTypes):
        raise ValueError(f"-p must not contain duplicate plot types!")
    if "genes" in args.plotTypes and args.annotationGFF3 == None:
        raise ValueError(f"Cannot plot gene locations without providing an --annotation GFF3 file!")
    if "histogram" in args.plotTypes and args.inputType == "depth":
        raise ValueError(f"Cannot plot histogram for -i depth data!")
    
    # Validate samples for coverage plot
    if "coverage" in args.plotTypes:
        for sampleID in args.sampleCoverage:
            if not sampleID in args.metadataDict["bulk1"] + args.metadataDict["bulk2"]:
                raise ValueError(f"Sample '{sampleID}' specified in --sampleCoverage not found in metadata!")
        if len(args.sampleCoverage) > NUM_SAMPLE_LINES:
            raise ValueError(f"Cannot plot more than {NUM_SAMPLE_LINES} samples using --sampleCoverage for clarity")
    
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

def derive_window_size(args, edDict):
    windowSize = None
    if args.inputType == "depth":
        # Derive window size from dictionary
        for key, posEDpairs in edDict.items():
            if len(posEDpairs[0]) > 1:
                windowSize = posEDpairs[0][1] - posEDpairs[0][0]
                break
        if windowSize != None:
            print(f"# Window size detected as {windowSize} bp from parsing of Euclidean distance file")
    if windowSize == None:
        if hasattr(args, "plotTypes") and "coverage" in args.plotTypes:
            if args.windowSize != None:
                print(f"# Window size specified as {args.windowSize} in cached parameters")
                windowSize = args.windowSize
            else:
                raise ValueError("Could not determine window size from Euclidean distance file or cached parameters!")
        else:
            windowSize = 0
    args.windowSize = windowSize

def raise_to_power(edDict, power):
    if power > 1:
        for key, posEDpairs in edDict.items():
            for i in range(len(posEDpairs[1])):
                posEDpairs[1][i] = posEDpairs[1][i] ** power

def locate_depth_files(args):
    DEPTH_DIR = os.path.join(args.workingDirectory, "depth")
    if not os.path.exists(DEPTH_DIR):
        raise FileNotFoundError(f"Depth directory '{DEPTH_DIR}' does not exist!")
    
    # Locate and validate depth files
    notFound = []
    depthFileDict = {"bulk1": [], "bulk2": []}
    for bulk, sampleList in args.metadataDict.items():
        for sample in sampleList:
            depthFile = os.path.join(DEPTH_DIR, f"{sample}.binned.{args.windowSize}.tsv")
            if not os.path.isfile(depthFile):
                notFound.append(sample)
            else:
                depthFileDict[bulk].append([sample, depthFile])
    if notFound != []: # for testing and development
        raise FileNotFoundError(f"Could not find depth files with bin size of {args.windowSize} " +
                                f"for samples: {', '.join(notFound)}")
    return depthFileDict

def main():
    usage = """%(prog)s processes the output of psQTL_proc.py to either 1) plot segregation
    statistics or 2) report on gene proximity to potential QTLs. The segregation statistics
    can be plotted in as a combination of line plots, scatter plots, histograms,
    and/or gene locations. The gene proximity report will identify genes that are proximal to or
    contained within potential QTLs (in the case of deletions). The input directory is expected
    to have been 'initialise'd by psQTL_prep.py and 'process'ed by psQTL_proc.py.
    """
    # Establish main parser
    p = argparse.ArgumentParser()
    
    # Set arguments shared by subparsers
    p.add_argument("-d", dest="workingDirectory",
                   required=True,
                   help="Specify the location where the analysis is being performed")
    p.add_argument("-f", dest="genomeFasta",
                   required=True,
                   help="Specify the location of the genome FASTA file")
    p.add_argument("-i", dest="inputType",
                   required=True,
                   choices=["depth", "call"],
                   help="""Specify whether you are analysing the Euclidean distance calculation
                   from variant 'call's or 'depth' predictions of deleted regions""")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="""Specify the location to write the output file; for 'plot', this must
                   end with '.pdf', '.png', or '.svg'; for 'report', this must end with
                   '.tsv' or '.csv'""")
    p.add_argument("--power", dest="power",
                   required=False,
                   type=int,
                   help="""Optionally, specify the power to raise Euclidean distances to
                   reduce noise (default: 4)""",
                   default=4)
    p.add_argument("--regions", dest="regions",
                   required=False,
                   nargs="+",
                   help="""Optionally, specify which regions to plot. Providing no input
                   to this argument will plot all chromosomes; otherwise, specify the
                   chromosome(s) to plot by their individual genome contig identifiers
                   (e.g., 'chr1') with the option to specify a range within the chromosome
                   with chr:start-end format (e.g., 'chr1:1000000-2000000').
                   """,
                   default=[])
    p.add_argument("--missing", dest="missingFilter",
                   type=float,
                   required=False,
                   help="""Optionally, specify the proportion of missing data that is
                   tolerated in both bulk populations before a variant is filtered out
                   (recommended: 0.5)""",
                   default=0.5)
    p.add_argument("-v", "--version",
                   action="version",
                   version="psQTL_post.py {version}".format(version=__version__))
    
    # Establish subparsers
    subParentParser = argparse.ArgumentParser(description=usage)
    subParentParser.add_argument("-v", "--version",
                                 action="version",
                                 version="psQTL_post.py {version}".format(version=__version__))
    
    subparsers = subParentParser.add_subparsers(dest="mode",
                                                required=True)
    
    pparser = subparsers.add_parser("plot",
                                    parents=[p],
                                    add_help=False,
                                    help="Plot segregation statistics")
    pparser.set_defaults(func=pmain)
    
    rparser = subparsers.add_parser("report",
                                    parents=[p],
                                    add_help=False,
                                    help="Report gene proximity to potential QTLs")
    rparser.set_defaults(func=rmain)
    
    # Plot-subparser arguments
    ## Required arguments
    pparser.add_argument("-p", dest="plotTypes",
                         required=True,
                         nargs="+",
                         choices=["line", "scatter", "histogram", "coverage", "genes"],
                         help="Specify one or more plot types to generate")
    ## Optional file arguments
    pparser.add_argument("--annotation", dest="annotationGFF3",
                         required=False,
                         help="""Optionally, specify the location of the genome annotation
                         GFF3 file if you want to plot gene locations""")
    ## Data arguments
    pparser.add_argument("--wma", dest="wmaSize",
                         type=int,
                         required=False,
                         help="""LINE PLOT: optionally, specify the number of previous
                         values to consider during weighted moving average
                         calculation (default: 5)""",
                         default=5)
    pparser.add_argument("--bin", dest="binSize",
                         type=int,
                         required=False,
                         help="""HISTOGRAM PLOT: Optionally, specify the bin size to 
                         count variants within (default: 100000)""",
                         default=100000)
    pparser.add_argument("--threshold", dest="binThreshold",
                         type=float,
                         required=False,
                         help="""HISTOGRAM PLOT: Optionally, specify the Euclidean
                         distance threshold for counting a variant within each bin
                         (default: 0.4)""",
                         default=0.4)
    pparser.add_argument("--sampleCoverage", dest="sampleCoverage",
                         required=False,
                         nargs="+",
                         help="""COVERAGE PLOT: Optionally, specify one or more samples
                         to plot coverage data as individual lines; these samples will
                         be omitted from the bulk group values""",
                         default=[])
    ## Style arguments
    pparser.add_argument("--width", dest="width",
                         type=int,
                         required=False,
                         help="""Optionally, specify the total output plot width
                         (default: calculated internally with 5 per region)""",
                         default=None)
    pparser.add_argument("--height", dest="height",
                         type=int,
                         required=False,
                         help="""Optionally, specify the output plot height
                         (default: calculated internally with 5 per plot type)""",
                         default=None)
    
    # Report-subparser arguments
    rparser.add_argument("-a", dest="annotationGFF3",
                         required=True,
                         help="Specify the location of the genome annotation GFF3 file")
    rparser.add_argument("--radius", dest="radiusSize",
                         type=int,
                         required=False,
                         help="""Optionally, specify the radius (in bp) surrounding a variant or
                         deletion window that you want to consider as being 'proximal' to
                         a gene (default: 50000)""",
                         default=50000)
    
    args = subParentParser.parse_args()
    validate_args(args)
    
    # Perform mode-specific validation
    "Validate upfront before we get into time-consuming parsing to frontload the error checking"
    if args.mode == "plot":
        print("## psQTL_post.py - Plot Euclidean Statistics ##")
        validate_p(args)
    elif args.mode == "report":
        print("## psQTL_post.py - Report Gene Proximity ##")
        validate_r(args)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }
    
    # Raise error if no contigs are found in genome FASTA
    if lengthsDict == {}:
        raise ValueError(f"No contigs found in genome FASTA '{args.genomeFasta}'; is it actually a FASTA file?")
    
    # Validate and impute regions
    validate_regions(args, lengthsDict)
    for contigID, start, end, reverse in args.regions:
        if end > lengthsDict[contigID]:
            raise ValueError(f"--region '{contigID, start, end}' end position is > contig length '{lengthsDict[contigID]}'!")
    
    # Parse input file into dictionary data structure or load pre-existing pickle
    PICKLE_FILE = args.inputFile.replace(".tsv.gz", f".{args.missingFilter}.pkl")
    if os.path.isfile(PICKLE_FILE):
        with open(PICKLE_FILE, "rb") as fileIn:
            edDict = pickle.load(fileIn)
    else:
        edDict = parse_ed_as_dict(args.inputFile, args.metadataDict, args.missingFilter)
        with open(PICKLE_FILE, "wb") as fileOut:
            pickle.dump(edDict, fileOut)
    
    # Obtain window size
    derive_window_size(args, edDict)
    
    # Raise Euclidean distances to the power specified by the user
    "Raising to power after pickling lets us reuse the pickled data with different power values"
    raise_to_power(edDict, args.power)
    
    # Convert dictionary to Euclidean distance NCLS data structure
    "EDNCLS cannot be pickled so we need to do it like file->dict->EDNCLS"
    edNCLS = convert_dict_to_edncls(edDict, args.windowSize)
    
    # Make sure all contigs in Euclidean distance data are in genome FASTA
    for contigID in edNCLS.contigs:
        if contigID not in lengthsDict:
            raise ValueError(f"Contig ID '{contigID}' from Euclidean distance data not found in -f FASTA!")
    
    # Drop any regions that are not in the Euclidean distance data
    passedRegions = []
    for contigID, start, end, reverse in args.regions:
        if contigID in edNCLS.contigs:
            passedRegions.append([contigID, start, end, reverse])
        else:
            print(f"WARNING: {contigID} not found in Euclidean distance data; it will be skipped")
    args.regions = passedRegions
    
    # Exit if no regions are left
    if args.regions == []:
        raise ValueError("No regions remain after filtering; set --regions to one or more valid regions " + 
                         "and ensure your Euclidean distance data and genome FASTA are compatible!")
    
    # Split into mode-specific functions
    if args.mode == "plot":
        pmain(args, edNCLS, lengthsDict)
    elif args.mode == "report":
        rmain(args, edNCLS)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def pmain(args, edNCLS, lengthsDict):
    STANDARD_DIMENSION = 5
    PLOT_DIR = os.path.join(args.workingDirectory, "plots")
    os.makedirs(PLOT_DIR, exist_ok=True)
    
    # Locate depth files if necessary
    if "coverage" in args.plotTypes:
        depthFileDict = locate_depth_files(args)
    else:
        depthFileDict = None
    
    # Get our labels for the plots
    rowLabels = [
        f"$ED^{args.power}$" if "scatter" in args.plotTypes or "line" in args.plotTypes else None,
        f"Num. variants with $ED^{args.power}$ ≥ {args.binThreshold}\n" + \
            f"in {args.binSize} bp windows" if "histogram" in args.plotTypes else None,
        "Median-normalised coverage" if "coverage" in args.plotTypes else None,
        "Representative models" if "genes" in args.plotTypes else None
    ]
    rowLabels = [label for label in rowLabels if label != None]
    colLabels = [f"{region[0]}:{region[1]}-{region[2]}" if region[3] == False
                 else f"{region[0]}:{region[2]}-{region[1]}" # if reversed
                 for region in args.regions]
    
    # Derive plot dimensions
    if args.width == None:
        args.width = STANDARD_DIMENSION * len(colLabels)
        print(f"# Calculated width as num. regions ({len(colLabels)}) multiplied by 5 = {args.width}")
    if args.height == None:
        args.height = STANDARD_DIMENSION * len(rowLabels)
        print(f"# Calculated height as num. plot types ({len(rowLabels)}) multiplied by 5 = {args.height}")
    
    # Set up the overall figure object
    fig, axs = plt.subplots(nrows=len(rowLabels), ncols=len(colLabels),
                            figsize=(args.width, args.height))
    axs = np.reshape(axs, (len(rowLabels), len(colLabels))) # ensure shape is as expected
    fig.tight_layout()
    
    # Set titles and labels
    for ax, label in zip(axs[0], colLabels):
        ax.set_title(label, fontweight="bold")
    for ax, label in zip(axs[:,0], rowLabels):
        ax.set_ylabel(label)
    for ax in axs[-1]:
        ax.set_xlabel(f"Chromosomal position", fontweight="bold")
    
    # Plot a line and/or scatter plot
    rowNum = 0
    if "line" in args.plotTypes or "scatter" in args.plotTypes:        
        linescatter(axs, rowNum, edNCLS, args.regions, args.wmaSize,
                    True if "line" in args.plotTypes else False,
                    True if "scatter" in args.plotTypes else False,
                    PLOT_DIR, rowNum+1 == len(rowLabels),
                    args.inputType)
        rowNum += 1
    
    # Plot a histogram
    if "histogram" in args.plotTypes:
        histogram(axs, rowNum, edNCLS, args.regions, args.binSize, args.binThreshold,
                  PLOT_DIR, rowNum+1 == len(rowLabels))
        rowNum += 1
    
    # Plot coverage data
    if "coverage" in args.plotTypes:
        # Parse coverage data
        coverageDict = parse_bins_as_dict(depthFileDict, args.windowSize)
        normalise_coverage_dict(coverageDict)
        depthNCLSDict = convert_dict_to_depthncls(coverageDict, args.windowSize)
        
        # Plot the parsed data
        coverage(axs, rowNum, depthNCLSDict, args.regions,
                 args.sampleCoverage, rowNum+1 == len(rowLabels))
        rowNum += 1
    
    # Plot gene locations
    if "genes" in args.plotTypes:
        genes(fig, axs, rowNum, args.gff3Obj, args.regions,
              rowNum+1 == len(rowLabels))
        rowNum += 1
    
    # Write plot to file
    fig.savefig(args.outputFileName, bbox_inches="tight")
    
    print("Plotting complete!")

def rmain(args, edNCLS):
    if args.inputType == "depth":
        report_depth(edNCLS, args.gff3Obj, args.regions, args.outputFileName, args.radiusSize)
    else:
        report_genes(edNCLS, args.gff3Obj, args.regions, args.outputFileName, args.radiusSize)
    
    print("Reporting complete!")

if __name__ == "__main__":
    main()
