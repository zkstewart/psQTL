#! python3
# heuristic_assessment_of_lines.py
# Automatically assesses the plot data (as found in the *.call_line.tsv files) for
# data features which were determined by empirical evaluation of plotted outputs.
# These data features should contribute to visual clarity when identifying the QTL
# at the known, simulated location.

import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from PIL import Image

# File location constants
RAW_FILE = "qtl_test.tsv"
SUMMARY_FILE = "qtl_summary.tsv"
PARENT_DIR = os.getcwd()
PLOTS_DIR = "param_plots"

# Metric constants
SNP_SPACING = 1000 # 1Kbp per SNP position
CHECK1_THRESHOLD = 2
CHECK2_RADIUS = 50 # 50 Kbp since using an index on the y list returns the position at y*SNP_SPACING
CHECK3_RADIUS = 500000 # 500 Kbp
CHECK4_THRESHOLD = 0.5

# Image constants
NUM_ROWS = 2
NUM_COLS = 6
IMG_WIDTH = 518
IMG_HEIGHT = 524

# Functions for heuristic metrics
def heuristic_1(y, centreY):
    if centreY > (y[0] * CHECK1_THRESHOLD) and centreY > (y[-1] * CHECK1_THRESHOLD):
        return True
    else:
        return False

def heuristic_2(y, centrePeak):
    if centrePeak == max(y):
        return True
    else:
        return False

def isOverlapping(start1, end1, start2, end2):
    """Does the range (start1, end1) overlap with (start2, end2)?"""
    return end1 >= start2 and end2 >= start1

def heuristic_3(y, centreIndex, centrePeak):
    maxPositions = [i*SNP_SPACING for i, j in enumerate(y) if j == centrePeak] # i*SNP_SPACING gives the x position
    centrePosition = centreIndex*SNP_SPACING
    if any([ isOverlapping(xCoord, xCoord, centrePosition-CHECK3_RADIUS, centrePosition+CHECK3_RADIUS) for xCoord in maxPositions ]):
        return False
    else:
        return True

def heuristic_4(centreY):
    if centreY > CHECK4_THRESHOLD:
        return True
    else:
        return False

# Functions for R^2 metric
def quadratic_func(x, a, b, c):
    return a*x**2 + b*x + c

def constrained_fit(x, y):
    # Establish weights to force the curve to go down to 0 at the start and end
    y[[0, -1]] = 0 # set first and last points artificially to 0
    sigma = np.ones(len(x))
    sigma[[0, -1]] = 0.0001 # weight first and last points heavily to force curve fitting through these points
    
    # Fit the curve
    popt, pcov = curve_fit(quadratic_func, x, y,
                           p0=(-0.01, 0, -0.01), # initial a, b, c values
                           sigma=sigma, # ensure line goes through y=0 at start and end
                           bounds=([-1, -1, -np.inf], [0, 1, np.inf])) # make sure a is negative
    
    # Calculate R^2
    residuals = y - quadratic_func(x, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    return x, y, popt, r_squared

def unconstrained_fit(x, y):
    # Fit the curve
    popt, pcov = curve_fit(quadratic_func, x, y,
                           p0=(-0.01, 0, -0.01), # initial a, b, c values
                           bounds=([-1, -1, -np.inf], [0, 1, np.inf])) # make sure a is negative
    
    # Calculate R^2
    residuals = y - quadratic_func(x, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    return x, y, popt, r_squared

def triangle_fit(x, y, centreIndex):
    # Get the triangle points
    minY = np.min(y)
    maxY = np.max(y)
    
    # Handle flat lines
    if minY == maxY:
        maxY += 0.1
    
    # Obtain the triangle data points
    slopeUp = np.linspace(minY, maxY, num=centreIndex)
    slopeDown = np.linspace(maxY, minY, num=centreIndex)
    triangleY = np.concatenate((slopeUp, slopeDown))
    
    # Calculate R^2
    residuals = y - triangleY
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    return triangleY, r_squared

def plot_fit(x, y, popt):
    plt.plot(x, quadratic_func(x, *popt), 'g--',
                label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
    plt.scatter(x, y, label='data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Curve Fit')
    plt.legend()
    plt.show()

def plot_triangle_fit(x, y, triangleY):
    plt.plot(x, triangleY, 'g--')
    plt.scatter(x, y, label='data')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Curve Fit')
    plt.legend()
    plt.show()

####

# Locate files containing plot data
fileGroups = {}
for root, dirs, files in os.walk(os.getcwd()):
    if len(files) == 1 and files[0] == "chr1.0-10000000.call_line.tsv":
        callLineFile = os.path.join(root, files[0])
        params = os.path.dirname(root).split("/", maxsplit=8)[-1].replace("/", "_")
        fileGroups.setdefault(params, [])
        seed = root.split("/")[-5]
        fileGroups[params].append([callLineFile, seed])

# Assess each parameter set
if not os.path.exists(RAW_FILE):
    header = ["pop_size", "pop_balance", "phenotype_error", "seed", "centre_exceeds_edges", "centre_is_peak",
              "centre_is_only_peak", "centre_exceeds_cutoff", "r_squared"]
    with open(RAW_FILE, "w") as fileOut:
        fileOut.write("\t".join(header) + "\n")
        for params, files in fileGroups.items():
            for file, seed in files:
                # Parse data from call line file
                x, y = [], []
                with open(file, "r") as fileIn:
                    for line in fileIn:
                        if not line.startswith("contigID"):
                            contigID, position, ed, smoothedEd = line.strip().split("\t")
                            position = int(position)
                            smoothedEd = float(smoothedEd)
                            
                            x.append(position)
                            y.append(smoothedEd)
                
                # Run heuristic assessments
                centreIndex = len(x) // 2
                centreY = y[centreIndex]
                centrePeak = max(y[centreIndex - CHECK2_RADIUS:centreIndex + CHECK2_RADIUS])
                
                check1 = heuristic_1(y, centreY) # centre exceeds the edges
                check2 = heuristic_2(y, centrePeak) # centre is a peak
                check3 = heuristic_3(y, centreIndex, centrePeak) # centre is the only peak
                check4 = heuristic_4(centreY) # centre ED exceeds a cutoff
                
                # Fit a curve to the data and check the R^2 value
                x = np.array(x)
                y = np.array(y)
                #x, y, popt, r_squared = unconstrained_fit(x, y)
                #x, y, popt, r_squared = constrained_fit(x, y)
                #plot_fit(x, y, popt)
                triangleY, r_squared = triangle_fit(x, y, centreIndex)
                #plot_triangle_fit(x, y, triangleY)
                
                # Output results
                fileOut.write("\t".join(params.split("_")))
                fileOut.write(f"\t{seed}\t{check1}\t{check2}\t{check3}\t{check4}\t{r_squared}\n")

# Assess results
WEAK = 0.25
MID = 0.5
STRONG = 0.75

df = pd.read_csv(RAW_FILE, sep="\t")
pops = sorted(df["pop_size"].unique())
balances = sorted(df["pop_balance"].unique())
errors = sorted(df["phenotype_error"].unique())

toCheck = {}
with open(SUMMARY_FILE, "w") as fileOut:
    fileOut.write("pop_size\tpop_balance\tphenotype_error\tnone\tweak\tmid\tstrong\n")
    for pop_size in pops:
        for pop_balance in balances:
            for phenotype_error in errors:
                paramString = f"{pop_size}_{pop_balance}_{phenotype_error}"
                
                # Parse out individual simulation result values
                paramsDf = df[(df["pop_size"] == int(pop_size)) & (df["pop_balance"] == float(pop_balance)) & (df["phenotype_error"] == float(phenotype_error))]
                if len(paramsDf) == 0:
                    continue
                
                # Calculate strength values
                noSignal = sum(paramsDf["r_squared"] < WEAK)
                weakSignal = sum((paramsDf["r_squared"] >= WEAK) & (paramsDf["r_squared"] < MID))
                midSignal = sum((paramsDf["r_squared"] >= MID) & (paramsDf["r_squared"] < STRONG))
                strongSignal = sum(paramsDf["r_squared"] >= STRONG)
                
                # Output summarised line
                fileOut.write(f"{pop_size}\t{pop_balance}\t{phenotype_error}\t{noSignal}\t{weakSignal}\t{midSignal}\t{strongSignal}\n")
                
                # Get a range of simulations from this param combination with different R^2 values
                toCheck[paramString] = {}
                lastCutoff = -1
                for cutoff in range(0, 11, 1):
                    cutoff = cutoff*0.1
                    found = paramsDf[(paramsDf["r_squared"] <= cutoff) & (paramsDf["r_squared"] > lastCutoff)]
                    
                    if len(found) == 0:
                        toCheck[paramString]["{:.1f}".format(lastCutoff)] = None
                    else:
                        toCheck[paramString]["{:.1f}".format(lastCutoff)] = found.iloc[0]["seed"]
                    
                    lastCutoff = cutoff

# Visualise plots
os.makedirs(PLOTS_DIR, exist_ok=True)
for paramKey, cutoffDict in toCheck.items():
    pop_size, pop_balance, phenotype_error = paramKey.split("_")
        
    # Obtain images for this parameter combination
    paramImages = []
    for cutoff, seed in cutoffDict.items():
        if seed == None:
            cutoffImage = Image.new('RGB', (IMG_WIDTH, IMG_HEIGHT))
        else:
            imageFile = os.path.join(PARENT_DIR, seed, pop_size, pop_balance, phenotype_error, f"{seed}_{paramKey}.png")
            cutoffImage = Image.open(imageFile)
        paramImages.append(cutoffImage)
    
    # Join all plots together into a grid
    concatImage = Image.new('RGB', (IMG_WIDTH*NUM_COLS, IMG_HEIGHT*NUM_ROWS))
    
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
    concatImage.save(os.path.join(PLOTS_DIR, f"{paramKey}.png"))

print("Program completed successfully!")
