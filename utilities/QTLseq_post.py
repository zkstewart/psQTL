#! python3
# QTLseq_post.py
# Allow comparable plotting of QTL-seq statistics to allow for comparison to psQTL

import os, argparse, re, sys, pickle
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from pycirclize.parser import Gff

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.plot import HorizontalPlot, CircosPlot
from modules.gff3 import GFF3
from modules.validation import validate_regions
from modules.parsing import read_gz_file
from modules.ed import convert_dict_to_windowed_ncls

# Plotting and parsing functions
def parse_qtlseq_as_dict(qtlseqFile):
    '''
    Parameters:
        qtlseqFile -- a string indicating the path to QTL-seq file
    Returns:
        deltaDict -- a dictionary with structure like:
                     {
                         "chr1": [[pos1, pos2, ...], [delta1, delta2, ...]],
                         "chr2": [[pos1, pos2, ...], [delta1, delta2, ...]],
                         ...
                     }
    '''
    HEADER = ["CHROM", "POSI", "variant", "bulk1_alleles", "bulk2_alleles", "euclideanDist"]
    
    # Parse the ED file
    deltaDict = {}
    starts, ends = [], []
    with read_gz_file(qtlseqFile) as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            # Optionally skip header
            "Not every QTL-seq index file has a header"
            if firstLine:
                firstLine = False
                if sl[0] == "CHROM":
                    continue
            
            # Parse relevant details and validate format
            chrom, posi, variant, bulk1_depth, bulk2_depth, p99, p95, \
                bulk1_SNPindex, bulk2_SNPindex, delta_SNPindex = sl
            try:
                posi = int(posi)
            except:
                raise ValueError(f"Position '{posi}' is not an integer; offending line is '{line}'")
            try:
                deltaValue = float(delta_SNPindex)
            except:
                raise ValueError(f"Delta SNP-index '{delta_SNPindex}' is not a float; offending line is '{line}'")
            
            # Store in dictionary
            deltaDict.setdefault(chrom, [[], []])
            deltaDict[chrom][0].append(posi)
            deltaDict[chrom][1].append(abs(deltaValue))
    return deltaDict

# Validation functions
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
            raise FileNotFoundError(f"--annotation file '{args.annotationGFF3}' is not a file!")
        else:
            if args.plotStyle == "horizontal":
                args.gff3Obj = GFF3(args.annotationGFF3) # parsing now to raise errors early
                args.gff3Obj.create_ncls_index("gene")
            elif args.plotStyle == "circos":
                parser = Gff(args.annotationGFF3)
                args.gff3Obj = parser.get_seqid2features(feature_type=None)
    
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
    if args.axisSpace < 0:
        raise ValueError(f"--axisSpace value '{args.axisSpace}' must be >= 0!")
    if args.axisSpace >= 360:
        raise ValueError(f"--axisSpace value '{args.axisSpace}' must be < 360!")
    
    # Alert user to Circos dimensions logic
    if args.plotStyle == "circos":
        if (args.width != None and args.height != None):
            if args.width != 8 or args.height != 8:
                print("# Note: psQTL is not yet equipped to handle non-default width or height while maintaining " +
                      "clarity for Circos plots, so width and height will be set to 8.")
                args.width = 8
                args.height = 8
            if args.width != args.height:
                print(f"# Note: Circos plotting enforces square dimensions, so width and height " +
                    f"will be set to the smallest value provided i.e., '{min(args.width, args.height)}'.")
    
    # Validate plot types
    if len(set(args.plotTypes)) != len(args.plotTypes):
        raise ValueError(f"-p must not contain duplicate plot types!")
    if "genes" in args.plotTypes and args.annotationGFF3 == None:
        raise ValueError(f"Cannot plot gene locations without providing an --annotation GFF3 file!")
    
    # Validate output file suffix
    if not (args.outputFileName.endswith(".pdf") or args.outputFileName.endswith(".png") or args.outputFileName.endswith(".svg")):
        raise ValueError(f"-o output file '{args.outputFileName}' must end with '.pdf', '.png', or '.svg'!")

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
                   choices=["line", "scatter", "genes"],
                   help="Specify one or more plot types to generate")
    p.add_argument("-s", dest="plotStyle",
                   required=True,
                   choices=["horizontal", "circos"],
                   help="Specify the style of plot to generate")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="""Specify the location to write the output file; this must
                   end with '.pdf', '.png', or '.svg'""")
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
    p.add_argument("--space", dest="axisSpace",
                   type=int,
                   required=False,
                   help="""CIRCOS PLOT: Optionally, specify the space (in degrees)
                   to allow for the Y-axis labels in the top centre of the plot
                   (default: 10)""",
                   default=10)
    
    args = p.parse_args()
    validate_args(args)
    
    # Perform mode-specific validation
    print("## QTLseq_post.py - Plot SNP-index Statistics ##")
    validate_p(args)
    
    # Get contig lengths from genome FASTA
    genomeRecords = SeqIO.parse(open(args.genomeFasta, 'r'), "fasta")
    lengthsDict = { record.id:len(record) for record in genomeRecords }
    if lengthsDict == {}:
        raise ValueError(f"No contigs found in genome FASTA '{args.genomeFasta}'; is it actually a FASTA file?")
    
    # Validate and impute regions
    args.regions = validate_regions(args, lengthsDict)
    
    # Parse input file into NCLS data structure
    deltaDict = parse_qtlseq_as_dict(args.inputFile)
    deltaNCLS = convert_dict_to_windowed_ncls(deltaDict)
    
    # Split into mode-specific functions
    "There are no alternate modes with this script, just leaving it like this for consistency with psQTL_post.py"
    pmain(args, deltaNCLS, lengthsDict)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def pmain(args, deltaNCLS, lengthsDict):
    # Establish plotting object
    if args.plotStyle == "horizontal":
        plotter = HorizontalPlot(args.regions, callED=deltaNCLS, depthED=None,
                                 callSPLSDA=None, depthSPLSDA=None, integratedSPLSDA=None,
                                 annotationGFF3=args.gff3Obj if "genes" in args.plotTypes else None,
                                 power=1, wmaSize=args.wmaSize, width=args.width, height=args.height,
                                 coverageSamples=None)
        plotter.plot(args.plotTypes, args.outputFileName)
    elif args.plotStyle == "circos":
        plotter = CircosPlot(args.regions, callED=deltaNCLS, depthED=None,
                                 callSPLSDA=None, depthSPLSDA=None, integratedSPLSDA=None,
                                 annotationGFF3=args.gff3Obj if "genes" in args.plotTypes else None,
                                 power=1, wmaSize=args.wmaSize, width=args.width, height=args.height,
                                 coverageSamples=None)
        plotter.axisSpace = args.axisSpace
        plotter.plot(args.plotTypes, args.outputFileName)
    print("Plotting complete!")

if __name__ == "__main__":
    main()
