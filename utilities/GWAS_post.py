#! python3
# GWAS_post.py
# Allow comparable plotting of GWAS statistics to allow for comparison to psQTL.
# Also just provides access to the prettiness and functionality of psQTL's plotting functions.

import os, argparse, re, sys, math
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from modules.plot import HorizontalPlot, CircosPlot
from modules.gff3 import GFF3Graph
from modules.validation import validate_regions
from modules.parsing import read_gz_file
from modules.ed import convert_dict_to_windowed_ncls

# Plotting and parsing functions
def parse_2col_tsv_as_dict(tsvFile, keyIsLeft=True):
    '''
    General purpose parser for 2-column TSV to render a dictionary pairing
    one column (key) to the other (value) according to the keyIsLeft behavioural boolean.
    
    Parameters:
        tsvFile -- a string indicating the location of the file to parse
        keyIsLeft -- (OPTIONAL) a boolean controlling whether the left column should be
                     the dictionary key (default; ==True) or if this should be inversed
                     (==False)
    Returns:
        tsvDict -- a dictionary pairing keys and values corresponding to TSV file column contents.
    '''
    tsvDict = {}
    with read_gz_file(tsvFile) as fileIn:
        for line in fileIn:
            delim = "\t" if "\t" in line else ","
            sl = line.strip().split(delim)
            if len(sl) != 2:
                raise ValueError(f"--chromRelabel file line should have two columns, but has {len(sl)}; offending line is '{line.rstrip()}'")
            
            if keyIsLeft:
                left, right = sl
            else:
                right, left = sl
            
            tsvDict[left] = right
    return tsvDict

def parse_gwas_as_dict(gwasFiles, chromHead, posHead, statHead, relabelDict=None):
    '''
    Parameters:
        gwasFile -- a list containing strings that indicate the path to a GWAS results file
        chromHead -- a string indicating the (case sensitive) header value
                     listing the chromosome identifiers
        posHead -- a string indicating the (case sensitive) header value listing
                   the positions of each variant
        statHead -- a string indicating the (case sensitive) header value listing
                    the statistic (e.g., the P-value) associated with each variant
        relabelDict -- (OPTIONAL) None if chromosomes should not be relabeled during
                       parsing, OR a dictionary with structure like:
                       {
                           "
                       }
    Returns:
        statDict -- a dictionary with structure like:
                    {
                        "chr1": [[pos1, pos2, ...], [stat1, stat2, ...]],
                        "chr2": [[pos1, pos2, ...], [stat1, stat2, ...]],
                        ...
                    }
    '''
    statDict = {}
    starts, ends = [], []
    for gwasFile in gwasFiles:
        with read_gz_file(gwasFile) as fileIn:
            firstLine = True
            for line in fileIn:
                sl = line.rstrip("\r\n").split("\t")
                
                # Parse header
                if firstLine:
                    firstLine = False
                    
                    # Find each column
                    try:
                        chromIndex = sl.index(chromHead)
                    except KeyError:
                        raise KeyError(f"Chromosome column header '{chromHead}' not found within GWAS file; actual header == '{sl}'")
                    
                    try:
                        posIndex = sl.index(posHead)
                    except KeyError:
                        raise KeyError(f"Position column header '{posHead}' not found within GWAS file; actual header == '{sl}'")
                    
                    try:
                        statIndex = sl.index(statHead)
                    except KeyError:
                        raise KeyError(f"Statistic column header '{statHead}' not found within GWAS file; actual header == '{sl}'")
                    
                    continue
                
                # Parse relevant details and validate format
                chrom, posi, stat = sl[chromIndex], sl[posIndex], sl[statIndex]
                if relabelDict != None:
                    try:
                        chrom = relabelDict[chrom]
                    except:
                        raise ValueError(f"Chromosome '{chrom}' has no match within the --chromRelabel file")
                try:
                    posi = int(posi)
                except:
                    raise ValueError(f"Position '{posi}' is not an integer; offending line is '{line}'")
                try:
                    statValue = float(stat)
                    statValue = -math.log10(statValue)
                except:
                    raise ValueError(f"Statistic '{stat}' is not a float; offending line is '{line}'")
                if not statValue >= 0:
                    raise ValueError(f"Statistic '{stat}' is a negative value (we expect something like a P-value which has a minimum of 0); offending line is '{line}'")
                
                # Store in dictionary
                statDict.setdefault(chrom, [[], []])
                statDict[chrom][0].append(posi)
                statDict[chrom][1].append(statValue)
    return statDict

# Validation functions
def validate_args(args):
    # Validate input file
    for inputFile in args.inputFiles:
        if not os.path.isfile(inputFile):
            raise FileNotFoundError(f"GWAS file '{inputFile}' does not exist!")
    
    # Validate genome FASTA file
    if not os.path.isfile(args.genomeFasta):
        raise FileNotFoundError(f"-f '{args.genomeFasta}' is not a file!")
    
    # Validate annotation GFF3 file
    if args.annotationGFF3 != None:
        if not os.path.isfile(args.annotationGFF3):
            raise FileNotFoundError(f"--annotation file '{args.annotationGFF3}' is not a file!")
        else:
            args.gff3Obj = GFF3Graph(args.annotationGFF3) # parsing now to raise errors early
            
            # Validate that GFF3 is not empty
            if args.gff3Obj.features == {}:
                raise ValueError(f"--annotation file '{args.annotationGFF3}' is empty?")
            
            # Validate that GFF3 contains genes and mRNAs
            if not "gene" in args.gff3Obj.ftypes or not "mRNA" in args.gff3Obj.ftypes:
                raise ValueError(f"--annotation file '{args.annotationGFF3}' does not contain 'gene' and/or 'mRNA' features; " +
                                 "psQTL requires these features to be present to make use of a GFF3 annotation file.")
            
            # Create NCLS index and perform QC
            args.gff3Obj.create_ncls_index("gene")
            args.gff3Obj.qc(typesToCheck=["gene", "mRNA"]) # prints out warnings if any issues found
    
    # Validate metadata file
    if args.chromosomesMetadataTsv != None:
        if not os.path.isfile(args.chromosomesMetadataTsv):
            raise FileNotFoundError(f"--chromRelabel file '{args.chromosomesMetadataTsv}' is not a file!")
        else:
            args.relabelDict = parse_2col_tsv_as_dict(args.chromosomesMetadataTsv)
            # the tsv is expected with old:new layout since that will have been used during chromosome relabelling;
            # however, this script will instead want the relabelling dict to have new:old layout so we can restore
            # contig IDs back to their original form. Hence, keyIsLeft == False above.
    else:
        args.relabelDict = None
    
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
    usage = """%(prog)s performs post-processing of GWAS output files to generate psQTL-styled plots. This might be useful
    for enabling comparison to psQTL, or simply as a pretty method of displaying the GWAS statistics. This utility script
    borrows all relevant functions of psQTL_post.py with changes made to accommodate the different statistics and file format.
    Default values are set for --chr and --pos and --stat expecting outputs as derived from SAIGE; you should change these if
    your software of choice produces an output table with different header values. Note that the output plot will list
    the statistic as the "ED^1" which reflects a quirk of the underlying psQTL code; you should probably manually edit
    the resulting plot to change this to "-log10(P)".
    """
    # Establish main parser
    p = argparse.ArgumentParser(description=usage)
    
    # Set required arguments
    p.add_argument("-i", dest="inputFiles",
                   required=True,
                   nargs="+",
                   help="""Specify the GWAS output file(s) to process; multiple can be given
                   if e.g., results are separated into one file per chromosome""")
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
    ## Input parsing arguments
    p.add_argument("--chr", dest="chromColumnHead",
                   required=False,
                   help="""Optionally, specify the header value for the column listing
                   chromosome identifiers (default=='CHR')""",
                   default="CHR")
    p.add_argument("--pos", dest="posColumnHead",
                   required=False,
                   help="""Optionally, specify the header value for the column listing
                   variant positions (default=='POS')""",
                   default="POS")
    p.add_argument("--stat", dest="statColumnHead",
                   required=False,
                   help="""Optionally, specify the header value for the column listing
                   statistical values (default=='p.value')""",
                   default="p.value")
    ## File arguments
    p.add_argument("--annotation", dest="annotationGFF3",
                   required=False,
                   help="""Optionally, specify the location of the genome annotation
                   GFF3 file if you want to plot gene locations""")
    p.add_argument("--chromRelabel", dest="chromosomesMetadataTsv",
                   required=False,
                   help="""Optionally, if genomic contigs were renamed when performing GWAS analysis
                   (e.g., PLINK2 has enforced chromosomes as integers like '1', '2', and so on)
                   specify the location of a two-column TSV listing the GWAS contig labels (left)
                   against the original contig IDs (right). This will allow for contig ID consistency between
                   the GWAS file (-i), the genome file (-f), and the annotation file (--annotation)
                   if applicable.""")
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
    validate_args(args) # sets args.relabelDict and args.gff3Obj
    
    # Perform mode-specific validation
    print("## GWAS_post.py - Plot GWAS Statistics ##")
    validate_p(args)
    
    # Get contig lengths from genome FASTA
    with read_gz_file(args.genomeFasta) as fileIn:
        genomeRecords = SeqIO.parse(fileIn, "fasta")
        lengthsDict = { record.id:len(record) for record in genomeRecords }
    if lengthsDict == {}:
        raise ValueError(f"No contigs found in genome FASTA '{args.genomeFasta}'; is it actually a FASTA file?")
    
    # Validate and impute regions
    args.regions = validate_regions(args.regions, "plot", args.plotStyle, lengthsDict)
    
    # Parse input file into NCLS data structure
    statDict = parse_gwas_as_dict(args.inputFiles, args.chromColumnHead, args.posColumnHead, args.statColumnHead, args.relabelDict)
    statNCLS = convert_dict_to_windowed_ncls(statDict)
    
    # Split into mode-specific functions
    "There are no alternate modes with this script, just leaving it like this for consistency with psQTL_post.py"
    pmain(args, statNCLS, lengthsDict)
    
    # Print completion flag if we reach this point
    print("Program completed successfully!")

def pmain(args, statNCLS, lengthsDict):
    # Establish plotting object
    if args.plotStyle == "horizontal":
        plotter = HorizontalPlot(args.regions, callED=statNCLS, depthED=None,
                                 callSPLSDA=None, depthSPLSDA=None, integratedSPLSDA=None,
                                 annotationGFF3=args.gff3Obj if "genes" in args.plotTypes else None,
                                 power=1, wmaSize=args.wmaSize, width=args.width, height=args.height,
                                 coverageSamples=None)
        plotter.plot(args.plotTypes, args.outputFileName)
    elif args.plotStyle == "circos":
        plotter = CircosPlot(args.regions, callED=statNCLS, depthED=None,
                                 callSPLSDA=None, depthSPLSDA=None, integratedSPLSDA=None,
                                 annotationGFF3=args.gff3Obj if "genes" in args.plotTypes else None,
                                 power=1, wmaSize=args.wmaSize, width=args.width, height=args.height,
                                 coverageSamples=None)
        plotter.axisSpace = args.axisSpace
        plotter.plot(args.plotTypes, args.outputFileName)
    print("Plotting complete!")

if __name__ == "__main__":
    main()
