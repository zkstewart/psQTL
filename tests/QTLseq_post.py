#! python3
# QTLseq_post.py
# Allow comparable plotting of QTL-seq statistics to allow for comparison to psQTL

import os, argparse, re, sys, pickle
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from snpindex import DeltaNCLS, parse_qtlseq_as_dict, convert_dict_to_deltancls

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.plot import HorizontalPlot
from modules.gff3 import GFF3
from modules.validation import validate_regions

def validate_args(args):
    # Validate input file
    if not os.path.isfile(args.inputFile):
        raise FileNotFoundError(f"QTL-seq file '{args.inputFile}' does not exist!")
    
    # Validate genome FASTA file
    if not os.path.isfile(args.genomeFasta):
        raise FileNotFoundError(f"-f '{args.genomeFasta}' is not a file!")
    
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
        raise ValueError(f"--threshold value '{args.binThreshold}' must be >= 0! (note that we use ABSOLUTE values)")
    
    # Validate plot types
    if len(set(args.plotTypes)) != len(args.plotTypes):
        raise ValueError(f"-p must not contain duplicate plot types!")
    if "genes" in args.plotTypes and args.annotationGFF3 == None:
        raise ValueError(f"Cannot plot gene locations without providing an --annotation GFF3 file!")
    
    # Validate output file suffix
    if not args.outputFileName.endswith(".pdf") and not args.outputFileName.endswith(".png"):
        raise ValueError(f"-o output file '{args.outputFileName}' must end with '.pdf' or '.png'!")

def main():
    usage = """%(prog)s performs post-processing of QTL-seq output files to generate plots for comparison
    with psQTL output files. Borrows all relevant functions of psQTL_post.py with changes made to accommodate
    the different statistics and file formats used by QTL-seq.
    """
    # Establish main parser
    p = argparse.ArgumentParser()
    
    # Set required arguments
    p.add_argument("-i", dest="inputFile",
                   required=True,
                   help="Specify the QTL-seq output file to process")
    p.add_argument("-f", dest="genomeFasta",
                   required=True,
                   help="Specify the location of the genome FASTA file")
    p.add_argument("-p", dest="plotTypes",
                   required=True,
                   nargs="+",
                   choices=["line", "scatter", "histogram", "genes"],
                   help="Specify one or more plot types to generate")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Specify the location to write the output file")
    # Set optional arguments
    ## File arguments
    p.add_argument("--annotation", dest="annotationGFF3",
                   required=False,
                   help="""Optionally, specify the location of the genome annotation
                   GFF3 file if you want to plot gene locations""")
    ## Data arguments
    p.add_argument("--wma", dest="wmaSize",
                   type=int,
                   required=False,
                   help="""LINE PLOT: optionally, specify the number of previous
                   values to consider during weighted moving average
                   calculation (default: 5)""",
                   default=5)
    p.add_argument("--bin", dest="binSize",
                   type=int,
                   required=False,
                   help="""HISTOGRAM PLOT: Optionally, specify the bin size to 
                   count variants within (default: 100000)""",
                   default=100000)
    p.add_argument("--threshold", dest="binThreshold",
                   type=float,
                   required=False,
                   help="""HISTOGRAM PLOT: Optionally, specify the ABSOLUTE delta
                   SNP-index threshold for counting a variant within each bin
                   (default: 0.4)""",
                   default=0.4)
    ## Style arguments
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
    p.add_argument("--width", dest="width",
                   type=int,
                   required=False,
                   help="""Optionally, specify the total output plot width
                   (default: calculated internally with 5 per region)""",
                   default=None)
    p.add_argument("--height", dest="height",
                   type=int,
                   required=False,
                   help="""Optionally, specify the output plot height
                   (default: calculated internally with 5 per plot type)""",
                   default=None)
    
    args = p.parse_args()
    validate_args(args)
    
    # Perform mode-specific validation
    print("## QTLseq_post.py - Plot SNP-index Statistics ##")
    validate_p(args)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }
    
    # Validate and impute regions
    args.regions = validate_regions(args, lengthsDict)
    
    # Parse input file into NCLS data structure
    deltaDict = parse_qtlseq_as_dict(args.inputFile)
    deltaNCLS = convert_dict_to_deltancls(deltaDict)
    
    # Make sure all contigs in QTLseq data are in genome FASTA
    for contigID in deltaNCLS.contigs:
        if contigID not in lengthsDict:
            raise ValueError(f"Contig ID '{contigID}' from QTL-seq data not found in -f FASTA!")
    
    # Split into mode-specific functions
    "There are no alternate modes with this script, just leaving it like this for consistency"
    pmain(args, deltaNCLS, lengthsDict)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def pmain(args, deltaNCLS, lengthsDict):
    args.inputType = ["call"]
    args.measurementTypes = ["ed"]
    
    # Establish plotting object
    plotter = HorizontalPlot(args.inputType, args.measurementTypes, args.plotTypes, args.regions,
                             args.wmaSize, args.binSize, args.binThreshold,
                             args.width, args.height)
    plotter.start_plotting()
    
    rowLabels = []
    colLabels = []
    
    # Generate plots grouped by plotTypes > inputType > measurementTypes
    alreadyPlottedScatter = False
    for pType in args.plotTypes:
        # Plot a line and/or scatter plot
        if pType in ["line", "scatter"] and not alreadyPlottedScatter:
            alreadyPlottedScatter = True
            rowLabels.append(f"abs. ΔSNP-index")
            plotter.plot_linescatter(deltaNCLS, deltaNCLS)
        # Plot a histogram
        elif pType == "histogram":
            rowLabels.append(f"Num. variants with abs. ΔSNP-index ≥ {args.binThreshold}\n" + \
                f"in {args.binSize} bp windows")
            plotter.plot_histogram(deltaNCLS)
        # Plot gene locations
        elif pType == "genes":
            plotter.plot_genes(args.gff3Obj)
            rowLabels.append("Representative models")
    
    # Format column labels
    colLabels = [f"{region[0]}:{region[1]}-{region[2]}" if region[3] == False
                 else f"{region[0]}:{region[2]}-{region[1]}" # if reversed
                 for region in args.regions]
    
    # Set the row and column labels
    plotter.set_row_labels(rowLabels)
    plotter.set_col_labels(colLabels)
    
    # Write plot to file
    plotter.savefig(args.outputFileName)
    
    print("Plotting complete!")

if __name__ == "__main__":
    main()
