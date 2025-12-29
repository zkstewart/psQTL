#! python3
# heuristic_assessment_of_lines.py
# Automatically assesses the plot data (as found in the *.call_line.tsv files) for
# data features which were determined by empirical evaluation of plotted outputs.
# These data features should contribute to visual clarity when identifying the QTL
# at the known, simulated location.

import os, argparse, math

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import concurrent.futures
from PIL import Image

# Image constants
NUM_ROWS = 2
NUM_COLS = 6
IMG_WIDTH = 518
IMG_HEIGHT = 524

# Heuristic assessment constants
WEAK_RSQ = 0.25
MID_RSQ = 0.5
STRONG_RSQ = 0.75

def validate_args(args):
    # Validate numeric args
    if args.threads < 1:
        raise ValueError("-t must be a positive integer")
    
    # Create output locations
    if not os.path.exists(args.plotsDirectory):
        os.makedirs(args.plotsDirectory, exist_ok=True)
        print(f"# Created '{os.path.abspath(args.plotsDirectory)}' as part of argument validation")

# Functions for parallel processing of R^2 metric
def calc_r_squared(y, ypred):
    '''
    Calculates R-squared for a line/curve fitting.
    
    Parameters:
        y -- a numpy array of measured data values
        ypred -- a numpy array of predicted data values
    '''
    residuals = y - ypred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    return r_squared

def triangle_fit(x, y):
    '''
    Fits a triangle shape to the x (chromosome position) and y (ED^4) data points. This is
    done as 1) a triangle with points at the minimum y value (left and right) with the maximum
    at the centre. It is also done as 2) the same concept but with a plateau at a local minimum
    on the left and right borders.
    
    The goal is to find a distinct and noticeable peak in the statistics occurring at the site
    where the simulated QTL exists i.e., in the centre of the chromosome. A simple line is
    used to avoid potential overfitting, and to conform to an intuitive sense of how a QTL
    should manifest visually. This trend is measured with R-squared.
    
    Parameters:
        x -- a numpy array of integer coordinates, indicating the location of SNPs
        y -- a numpy array of numeric values for the ED^4 segregation of the SNPs
    '''
    # Get the triangle points
    minY = np.min(y)
    maxY = np.max(y)
    midY = (minY + maxY) / 2
    
    centreIndex = len(x) / 2
    quarterIndex = centreIndex / 2
    
    # Handle flat lines
    if (minY+0.1) >= maxY: # we need a noticeable difference between min and max for QTL detection
        maxY += 0.1
    
    # Triangle 1: full range peak (^)
    slopeUp = np.linspace(minY, maxY, num=math.floor(centreIndex))
    slopeDown = np.linspace(maxY, minY, num=math.ceil(centreIndex))
    fullTriangleY = np.concatenate((slopeUp, slopeDown))
    fullRsquared = calc_r_squared(y, fullTriangleY)
    
    # Triangle 2: subrange peak (_^_)
    leftFlatIndex = 0
    while leftFlatIndex < quarterIndex: # only plateau up to 1/4 into the 'plot'
        if y[leftFlatIndex] > midY:
            break
        leftFlatIndex += 1 # check the next position in this quadrant
    
    if leftFlatIndex != 0: # this is 0 if the starting position is >= midY
        leftAvgY = np.mean(y[0:leftFlatIndex])
    else:
        leftAvgY = minY
    slopeUp = np.concatenate((
        np.array([ leftAvgY for _ in range(leftFlatIndex)]), # plateau
        np.linspace(leftAvgY, maxY, num=math.floor(centreIndex) - leftFlatIndex) # peak (incline)
    ))
    
    rightFlatIndex = len(x)-1
    while rightFlatIndex > (math.ceil(centreIndex) + quarterIndex): # only plateau for the last 1/4 of the 'plot'
        if y[rightFlatIndex] > midY:
            break
        rightFlatIndex -= 1 # crawl back into this quadrant
    
    if rightFlatIndex != len(x)-1: # this is the final index if the ending position is >= midY
        rightAvgY = np.mean(y[rightFlatIndex:])
    else:
        rightAvgY = minY
    slopeDown = np.concatenate((
        np.linspace(maxY, rightAvgY, num=rightFlatIndex - math.floor(centreIndex)), # peak (decline)
        np.array([ rightAvgY for _ in range(len(x) - rightFlatIndex)]) # plateau
    ))
    
    subrangeTriangleY = np.concatenate((slopeUp, slopeDown))
    subrangeRsquared = calc_r_squared(y, subrangeTriangleY)
    
    # Return the best R^2 value
    if (fullRsquared >= subrangeRsquared) or pd.isna(subrangeRsquared):
        return fullTriangleY, fullRsquared
    else:
        return subrangeTriangleY, subrangeRsquared

def fit_process_handler(lineFileName):
    '''
    Multiprocessing function to handle a line file and output the R-squared value obtained
    through simple line fitting.
    
    Parameters:
        lineFileName -- a string indicating the location of text file listing line plot data.
    '''
    # Parse data from line file
    x, y = [], []
    with open(lineFileName, "r") as fileIn:
        for line in fileIn:
            if not line.startswith("contigID"): # skip header
                contigID, position, ed, smoothedEd = line.strip().split("\t")
                position = int(position)
                smoothedEd = float(smoothedEd)
                x.append(position)
                y.append(smoothedEd)
    
    # Fit a triangle to the data and obtain the R^2 value
    x = np.array(x)
    y = np.array(y)
    try:
        triangleY, r_squared = triangle_fit(x, y)
    except ValueError:
        raise ValueError(f"Error when running triangle_fit() for '{lineFileName}'")
    return r_squared

def plot_triangle_fit(x, y, triangleY):
    '''
    Method left in for manual testing and visualisation of how the line fits
    the ED^4 curve.
    '''
    plt.plot(x, triangleY, "g--")
    plt.scatter(x, y, label="SNP segregation")
    plt.xlabel("Chromosome position")
    plt.ylabel("$ED^4$")
    plt.title('Curve Fit')
    plt.legend()
    plt.show()

def main():
    usage = """%(prog)s assesses the results of auto_psqtl_runner.py
    to roughly assess whether the simulated data analysis gave results that could
    enable QTL identification. This script was made for use in psQTL
    validation as depicted in the publication associated with this software.
    
    Note that this script assumes the simulation results are nested within
    the directory this script is being called from ($PWD).
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-t", dest="threads",
                   required=True,
                   type=int,
                   help="""Specify the number of threads to use with parallel steps;
                   actual threads used is -t+1 for the file output writing thread""")
    p.add_argument("--raw", dest="rawFileName",
                   required=False,
                   help="Output file name for each raw (replicated) result",
                   default="qtl_test.tsv")
    p.add_argument("--summary", dest="summaryFileName",
                   required=False,
                   help="Output file name for summarised results for parameter combinations",
                   default="qtl_summary.tsv")
    p.add_argument("--plots", dest="plotsDirectory",
                   required=False,
                   help="Directory to write summarised results for parameter combinations",
                   default="param_plots")
    
    args = p.parse_args()
    validate_args(args)
    
    # Locate files containing plot data
    fileGroups = {}
    for root, dirs, files in os.walk(os.getcwd()):
        if len(files) == 1 and files[0] == "chr1.0-10000000.call_line.tsv":
            callLineFile = os.path.join(root, files[0])
            params = os.path.dirname(root).split("/", maxsplit=8)[-1].replace("/", "_")
            fileGroups.setdefault(params, [])
            seed = root.split("/")[-5]
            fileGroups[params].append([callLineFile, seed])
    
    # Parallel processing for R^2 measure
    if not os.path.exists(args.rawFileName):
        with open(args.rawFileName, "w") as fileOut:
            header = ["pop_size", "pop_balance", "phenotype_error", "seed", "r_squared"]
            fileOut.write("\t".join(header) + "\n")
            
            with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
                for params, files in fileGroups.items():
                    paramsString = "\t".join(params.split("_"))
                    
                    callLineFiles, seeds = map(list, zip(*files)) # pythonic unzip
                    rSquaredList = executor.map(fit_process_handler, callLineFiles)
                    for r_squared, seed in zip(rSquaredList, seeds):
                        fileOut.write(f"{paramsString}\t{seed}\t{r_squared}\n")
    
    # Summarise results
    df = pd.read_csv(args.rawFileName, sep="\t")
    pops = sorted(df["pop_size"].unique())
    balances = sorted(df["pop_balance"].unique())
    errors = sorted(df["phenotype_error"].unique())
    
    toPlot = {} # dict to retain exemplar plots for manual visualisation of R^2 cutoff points
    with open(args.summaryFileName, "w") as fileOut:
        fileOut.write("pop_size\tpop_balance\tphenotype_error\tnone\tweak\tmid\tstrong\n")
        for pop_size in pops:
            for pop_balance in balances:
                for phenotype_error in errors:
                    paramString = f"{pop_size}_{pop_balance}_{phenotype_error}"
                    
                    # Parse out simulation results for this param combination
                    paramsDf = df[(df["pop_size"] == int(pop_size)) & (df["pop_balance"] == float(pop_balance)) & (df["phenotype_error"] == float(phenotype_error))]
                    if len(paramsDf) == 0:
                        continue
                    
                    # Calculate signal strength values
                    noSignal = sum((paramsDf["r_squared"] < WEAK_RSQ) | (pd.isna(paramsDf["r_squared"])))
                    weakSignal = sum((paramsDf["r_squared"] >= WEAK_RSQ) & (paramsDf["r_squared"] < MID_RSQ))
                    midSignal = sum((paramsDf["r_squared"] >= MID_RSQ) & (paramsDf["r_squared"] < STRONG_RSQ))
                    strongSignal = sum(paramsDf["r_squared"] >= STRONG_RSQ)
                    
                    # Output summarised line
                    fileOut.write(f"{pop_size}\t{pop_balance}\t{phenotype_error}\t{noSignal}\t{weakSignal}\t{midSignal}\t{strongSignal}\n")
                    
                    # Get a range of simulations from this param combination with different R^2 values
                    "This lets us visualise the distribution of R^2 values and validate that the results make sense"
                    toPlot[paramString] = {}
                    found = paramsDf[pd.isna(paramsDf["r_squared"])]
                    if len(found) == 0:
                        toPlot[paramString]["nan"] = None
                    else:
                        toPlot[paramString]["nan"] = found.iloc[0]["seed"]
                    
                    lastCutoff = -math.inf
                    for cutoff in range(0, 11, 1):
                        cutoff = cutoff*0.1
                        found = paramsDf[(paramsDf["r_squared"] <= cutoff) & (paramsDf["r_squared"] > lastCutoff)]
                        
                        if len(found) == 0:
                            toPlot[paramString]["{:.1f}".format(lastCutoff)] = None
                        else:
                            toPlot[paramString]["{:.1f}".format(lastCutoff)] = found.iloc[0]["seed"]
                        
                        lastCutoff = cutoff
    
    # Visualise plots for manual assessment of R^2 measurement cutoffs
    for paramKey, cutoffDict in toPlot.items():
        pop_size, pop_balance, phenotype_error = paramKey.split("_")
        
        # Obtain images for this parameter combination
        paramImages = []
        for cutoff, seed in cutoffDict.items():
            if seed == None: # i.e., if no replicate obtained an R^2 within this cutoff range
                cutoffImage = Image.new("RGB", (IMG_WIDTH, IMG_HEIGHT)) # this will be a blank/black square
            else:
                imageFile = os.path.join(os.getcwd(), seed, pop_size, pop_balance, phenotype_error, f"{seed}_{paramKey}.png")
                cutoffImage = Image.open(imageFile)
            paramImages.append(cutoffImage)
        
        # Join all plots together into a grid
        concatImage = Image.new("RGB", (IMG_WIDTH*NUM_COLS, IMG_HEIGHT*NUM_ROWS))
        
        x_offset = 0
        y_offset = 0
        ongoingCount = 0
        for image in paramImages:
            concatImage.paste(image, (x_offset, y_offset))
            x_offset += IMG_WIDTH
            
            ongoingCount += 1
            if ongoingCount % NUM_COLS == 0:
                x_offset = 0
                y_offset += IMG_HEIGHT
        
        # Output concatenated plot
        concatImage.save(os.path.join(args.plotsDirectory, f"{paramKey}.png"))
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
