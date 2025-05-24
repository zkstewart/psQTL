#! python3
# psQTL_post.py
# Represents step 3 of the psQTL pipeline, which is to 'post'-process the data
# generated from psQTL_proc.py. It can plot segregation statistics with a plug-and-play
# approach of several different plot types (line, scatter, coverage, genes) or it
# can report on genes that are proximal to or contained within potential QTLs.

import os, argparse, sys, pickle
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from modules.validation import validate_post_args, validate_regions, validate_depth_files, \
                               validate_p, validate_r
from modules.depth import parse_bins_as_dict, normalise_coverage_dict, convert_dict_to_depthncls
from modules.ed import parse_ed_as_dict, convert_dict_to_windowed_ncls
from modules.splsda import parse_selected_to_windowed_ncls, parse_ber_to_windowed_ncls
from modules.plot import HorizontalPlot, CircosPlot
from modules.reporting import report_genes, report_depth
from _version import __version__

def derive_window_size(args, edDict):
    windowSize = None
    for key, posEDpairs in edDict.items():
        if len(posEDpairs[0]) > 1:
            windowSize = posEDpairs[0][1] - posEDpairs[0][0]
            break
    return windowSize

def raise_to_power(edDict, power):
    if power > 1:
        for key, posEDpairs in edDict.items():
            for i in range(len(posEDpairs[1])):
                posEDpairs[1][i] = posEDpairs[1][i] ** power

def main():
    usage = """%(prog)s processes the output of psQTL_proc.py to either 1) plot segregation
    statistics or 2) report on gene proximity to potential QTLs. The segregation statistics
    can be plotted in as a combination of line plots, scatter plots, alignment coverage plots,
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
    p.add_argument("-r", dest="resultTypes",
                   required=True,
                   nargs="+",
                   choices=["call", "depth"],
                   help="""Specify whether you are analysing the results from variant 'call's
                   and/or 'depth' predictions of deleted regions""")
    p.add_argument("-m", dest="measurementTypes",
                   required=True,
                   nargs="+",
                   choices=["ed", "splsda"],
                   help="""Specify whether you are analysing 'ed' (Euclidean distance)
                   and/or 'splsda' (Sparse Partial Least Squares Discriminant Analysis)
                   measurements""")
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
                   tolerated in either bulk population before a variant is filtered out
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
                         choices=["line", "scatter", "coverage", "genes"],
                         help="Specify one or more plot types to generate")
    pparser.add_argument("-s", dest="plotStyle",
                         required=True,
                         choices=["horizontal", "circos"],
                         help="Specify the style of plot to generate")
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
    locations = validate_post_args(args) # sets args.metadataDict; args.gff3Obj if relevant
    
    # Perform mode-specific validation
    "Validate upfront before we get into time-consuming parsing to frontload the error checking"
    if args.mode == "plot":
        print("## psQTL_post.py - Plot QTL Statistics ##")
        validate_p(args) # sets args.depthFileDict if relevant
    elif args.mode == "report":
        print("## psQTL_post.py - Report Gene Proximity ##")
        validate_r(args)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }
    if lengthsDict == {}:
        raise ValueError(f"No contigs found in genome FASTA '{args.genomeFasta}'; is it actually a FASTA file?")
    
    # Validate and impute regions
    args.regions = validate_regions(args, lengthsDict)
    
    # Parse result and measurement type data
    dataDict = {}
    if "call" in args.resultTypes:
        dataDict["call"] = {}
        if "ed" in args.measurementTypes:
            # Parse the Euclidean distance data
            pickleFile = locations.variantEdPickleFile(args.missingFilter)
            if os.path.isfile(pickleFile) and os.path.isfile(pickleFile + ".ok"):
                with open(pickleFile, "rb") as fileIn:
                    dataDict["call"]["ed"] = pickle.load(fileIn)
            else:
                dataDict["call"]["ed"] = parse_ed_as_dict(locations.variantEdFile, args.metadataDict, args.missingFilter)
                with open(pickleFile, "wb") as fileOut:
                    pickle.dump(dataDict["call"]["ed"], fileOut)
                open(pickleFile + ".ok", "w").close()
            
            # Raise Euclidean distances to the power specified by the user
            "Raising to power after pickling lets us reuse the pickled data with different power values"
            raise_to_power(dataDict["call"]["ed"], args.power)
            
            # Convert dictionary to Euclidean distance NCLS data structure
            "WindowedNCLS cannot be pickled so we need to do it like file->dict->WindowedNCLS"
            dataDict["call"]["ed"] = convert_dict_to_windowed_ncls(dataDict["call"]["ed"], 0) # windowSize = 0
        
        if "splsda" in args.measurementTypes:
            # Parse the Sparse Partial Least Squares Discriminant Analysis data
            dataDict["call"]["selected"] = parse_selected_to_windowed_ncls(locations.variantSplsdaSelectedFile)
            dataDict["call"]["ber"], dataDict["call"]["ber_windowSize"] = parse_ber_to_windowed_ncls(locations.variantSplsdaBerFile)
    
    if "depth" in args.resultTypes:
        dataDict["depth"] = {}
        if "ed" in args.measurementTypes:
            # Parse the Euclidean distance data
            pickleFile = locations.deletionEdPickleFile(args.missingFilter)
            if os.path.isfile(pickleFile) and os.path.isfile(pickleFile + ".ok"):
                with open(pickleFile, "rb") as fileIn:
                    dataDict["depth"]["ed"] = pickle.load(fileIn)
            else:
                dataDict["depth"]["ed"] = parse_ed_as_dict(locations.deletionEdFile, args.metadataDict, args.missingFilter)
                with open(pickleFile, "wb") as fileOut:
                    pickle.dump(dataDict["depth"]["ed"], fileOut)
                open(pickleFile + ".ok", "w").close()
            
            # Obtain window size
            args.windowSize = derive_window_size(args, dataDict["depth"]["ed"])
            
            # Raise Euclidean distances to the power specified by the user
            raise_to_power(dataDict["depth"]["ed"], args.power)
            
            # Convert dictionary to Euclidean distance NCLS data structure
            dataDict["depth"]["ed"] = convert_dict_to_windowed_ncls(dataDict["depth"]["ed"], args.windowSize)
        
        if "splsda" in args.measurementTypes:
            # Parse the Sparse Partial Least Squares Discriminant Analysis data
            dataDict["depth"]["selected"] = parse_selected_to_windowed_ncls(locations.deletionSplsdaSelectedFile)
            dataDict["depth"]["ber"], dataDict["depth"]["ber_windowSize"] = parse_ber_to_windowed_ncls(locations.deletionSplsdaBerFile)
    
    # Parse depth data if necessary
    if "coverage" in args.plotTypes and "depth" in args.resultTypes:
        depthFileDict = validate_depth_files(locations.depthDir, args.metadataDict, args.windowSize)
        coverageDict = parse_bins_as_dict(depthFileDict, args.windowSize)
        normalise_coverage_dict(coverageDict)
        dataDict["depth"]["ncls"] = convert_dict_to_depthncls(coverageDict, args.windowSize)
    
    # Split into mode-specific functions
    if args.mode == "plot":
        pmain(args, locations, dataDict)
    elif args.mode == "report":
        rmain(args, locations, dataDict)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def pmain(args, locations, dataDict):
    # Establish plotting object
    if args.plotStyle == "horizontal":
        plotter = HorizontalPlot(args.regions,
            callED=dataDict["call"]["ed"] if "call" in dataDict and "ed" in dataDict["call"] else None,
            depthED=dataDict["depth"]["ed"] if "depth" in dataDict and "ed" in dataDict["depth"] else None,
            callSPLSDA=(dataDict["call"]["selected"], dataDict["call"]["ber"]) \
                if "call" in dataDict and "selected" in dataDict["call"] else None,
            depthSPLSDA=(dataDict["depth"]["selected"], dataDict["depth"]["ber"]) \
                if "depth" in dataDict and "selected" in dataDict["depth"] else None,
            coverageNCLSDict=dataDict["depth"]["ncls"] if "depth" in dataDict and "ncls" in dataDict["depth"] else None,
            annotationGFF3=args.gff3Obj if "genes" in args.plotTypes else None,
            power=args.power, wmaSize=args.wmaSize, width=args.width, height=args.height)
        plotter.plot(args.plotTypes, args.outputFileName)
    elif args.plotStyle == "circos":
        plotter = CircosPlot(args.regions,
            callED=dataDict["call"]["ed"] if "call" in dataDict and "ed" in dataDict["call"] else None,
            depthED=dataDict["depth"]["ed"] if "depth" in dataDict and "ed" in dataDict["depth"] else None,
            callSPLSDA=(dataDict["call"]["selected"], dataDict["call"]["ber"]) \
                if "call" in dataDict and "selected" in dataDict["call"] else None,
            depthSPLSDA=(dataDict["depth"]["selected"], dataDict["depth"]["ber"]) \
                if "depth" in dataDict and "selected" in dataDict["depth"] else None,
            coverageNCLSDict=dataDict["depth"]["ncls"] if "depth" in dataDict and "ncls" in dataDict["depth"] else None,
            annotationGFF3=args.gff3Obj if "genes" in args.plotTypes else None,
            power=args.power, wmaSize=args.wmaSize, width=args.width, height=args.height)
        plotter.plot(args.plotTypes, args.outputFileName)
    
    print("Plotting complete!")

def rmain(args, locations, dataDict):
    if args.resultTypes == "depth":
        report_depth(edNCLS, args.gff3Obj, args.regions, args.outputFileName, args.radiusSize)
    else:
        report_genes(edNCLS, args.gff3Obj, args.regions, args.outputFileName, args.radiusSize)
    
    print("Reporting complete!")

if __name__ == "__main__":
    main()
