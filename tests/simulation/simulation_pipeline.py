#! python3
# simulation_pipeline.py
# Runs a pipeline involving three key steps: 1) simulate genotype calls for a population that
# segregates over a single QTL SNP, then 2) randomly sample data in a way that simulates population
# parameters including size, phenotype balance, and error, then 3) evaluate the results to determine
# whether successful QTL prediction would be possible using psQTL.

import os, sys, random, argparse, pickle, math
import numpy as np
import pandas as pd

from itertools import product
from math import sqrt, ceil

import matplotlib.pyplot as plt
from PIL import Image
import matplotlib.ticker as ticker

from concurrent.futures import ProcessPoolExecutor
from chromax import Simulator
from chromax.sample_data import genetic_map, genome

# Set global variables
NREP = 1000 # 1000 replicates for each parameter combination
MAIN_SEED = 1 # seed number is arbitrary; set to ensure reproducibility
random.seed(MAIN_SEED)
np.random.seed(MAIN_SEED)

# Image constants
NUM_ROWS = 2
NUM_COLS = 6
IMG_WIDTH = 525
IMG_HEIGHT = 525
IMG_DPI = 150
LINE_STYLES = ["solid", "dotted", "dashed"]
LINE_COLOURS = ["#348ABD", "#7A68A6", "#A60628", "#467821", "#CF4457"]
plt.style.use("ggplot")

# Heuristic assessment constants
WEAK_RSQ = 0.25
MID_RSQ = 0.5
STRONG_RSQ = 0.75

# Define functions
def validate_args(args):
    if args.threads < 1:
        raise ValueError("-t must be a positive integer")

def establish_directories(size, balance, phePct, symmetry):
    sizeOutDir = os.path.join(args.outputDirectory, str(size))
    os.makedirs(sizeOutDir, exist_ok=True)
    
    balanceOutDir = os.path.join(sizeOutDir, str(balance))
    os.makedirs(balanceOutDir, exist_ok=True)
    
    errorOutDir = os.path.join(balanceOutDir, str(phePct))
    os.makedirs(errorOutDir, exist_ok=True)
    
    finalOutDir = os.path.join(errorOutDir, symmetry)
    os.makedirs(finalOutDir, exist_ok=True)
    
    return os.path.join(str(size), str(balance), str(phePct), symmetry)

def random_sample(arr: np.array, size: int = 1) -> np.array:
    if size < 1:
        return []
    else:
        return arr[np.random.choice(len(arr), size=size, replace=False)]

def experiment_segregation(bulk1, bulk2, numSNPs=10000):
    '''
    Handles the calculation of ED^4 and subsequent line fitting to get R^2 for each
    process in the ProcessPool.
    
    Parameters:
        bulk1 / bulk2 -- a numpy array with shape like (num_samples, nump_snps, 2) where
                         2 refers to the alleles which are diploid
        numSNPs -- the number of SNPs that were simulated
    Returns:
        numAllelesG1 -- the number of genotyped alleles in group 1
        numAllelesG2 -- the number of genotyped alleles in group 2
        edist -- a float of the the Euclidean distance between the two groups
    '''
    # Reshape the data to have each 'row' be a SNP containing allele values (n=num_samples*2)
    numAlleles1 = len(bulk1)*2
    numAlleles2 = len(bulk2)*2
    
    b1Alleles = np.moveaxis(bulk1, 0, -1)
    b1Alleles = b1Alleles.reshape(numSNPs, numAlleles1)
    
    b2Alleles = np.moveaxis(bulk2, 0, -1)
    b2Alleles = b2Alleles.reshape(numSNPs, numAlleles2)
    
    # Calculate the ED^4 statistic for each SNP
    '''
    This uses modified logic of the main psQTL ed.py module calculate_allele_frequency_ed()
    function to receive numpy arrays of alleles for the entire chromosome rather than the
    list of lists for a single SNP it normally expects.
    '''
    distances = []
    for alleles1, alleles2 in zip(b1Alleles, b2Alleles):
        b1Count = { 1: int(np.sum(alleles1)), 0: len(alleles1) - int(np.sum(alleles1)) }
        b2Count = { 1: int(np.sum(alleles2)), 0: len(alleles2) - int(np.sum(alleles2)) }
        
        edist = sqrt(sum([
            ((b1Count[allele] / numAlleles1) - (b2Count[allele] / numAlleles2))**2
            for allele in [0, 1]
        ]))
        distances.append(edist**4)
    
    return distances

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
    diffY = max(maxY*0.1, minY*0.50) # account for flat lines by enforcing some difference between min and max
    if diffY < 1e-3:
        diffY = 1e-3 # mitigate issues with extremely low ED^4 values
    
    if (minY+diffY) >= maxY: # we need a noticeable difference between min and max for QTL detection
        maxY += diffY
    
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

def plot_triangle_fit(x, y, triangleY, r_squared):
    '''
    Obtains a plot (as a Pillow Image object) for visualising the line fit to the ED^4 data.
    
    Credit to https://www.icare.univ-lille.fr/how-to-convert-a-matplotlib-figure-to-a-numpy-array-or-a-pil-image/
    '''
    # Plot the image with line fit
    fig = plt.figure(figsize=(IMG_WIDTH/IMG_DPI, IMG_HEIGHT/IMG_DPI), dpi=IMG_DPI,
                     num=1, clear=True) # prevent memory leak
    ax = fig.add_subplot()
    
    ax.plot(x, triangleY, "g--")
    ax.scatter(x, y, label="SNP segregation")
    ax.set_xlabel("Chromosome position")
    ax.set_ylabel("$ED^4$")
    ax.set_title(round(r_squared, 4))
    
    # Convert to PIL Image object
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    buf = np.frombuffer(fig.canvas.tostring_argb(), dtype=np.uint8)
    buf.shape = (w, h, 4)
    buf = np.roll(buf, 3, axis = 2)
    
    plotImg = Image.frombytes("RGBA", (w , h), buf.tobytes())
    return plotImg

def simulate_f1(geneticMap, toKeep, seed, PARENTS_NPY, numOffspring=10000, numSNPs=10000):
    '''
    Parameters:
        geneticMap -- 
        toKeep -- a string of "good" or "bad" to indicate which progeny we want to accumulate
        seed -- an integer indicating the seed to use during chromax simulation
        PARENTS_NPY -- a string indicating the location of the parental genomes to use
                       stored as a .npy file
    '''
    # Centre the QTL position in the middle of the genome/chromosome
    qtlSNP = numSNPs//2
    if qtlSNP % 1 != 0:
        raise ValueError(f"The QTL SNP index ({qtlSNP}) must be a whole number. Please check the genome length and SNP MB values " +
                         "to make sure they divide evenly.")
    
    # Init the simulator
    simulator = Simulator(genetic_map=geneticMap, seed=seed)#, backend="cpu")
    
    # Generate the F1 population
    f0 = simulator.load_population(PARENTS_NPY)
    f1, _ = simulator.random_crosses(f0, 1, n_offspring=numOffspring) # returns (n_crosses, n_individuals, n_loci, n_alleles)
    f1 = f1.reshape(numOffspring, numSNPs, 2) # reshape to (n_individuals, n_loci, n_alleles)
    
    # Separate out progeny based on qtlSNP occurrence and save to file
    goodF1 = np.array(f1[(f1[:, qtlSNP, 0] == 1) & (f1[:, qtlSNP, 1] == 1)])
    badF1 = np.array(f1[(f1[:, qtlSNP, 0] == 0) | (f1[:, qtlSNP, 1] == 0)]) # this prevents errors during ProcessPoolExecutor below
    if toKeep == "good":
        return goodF1
    else:
        return badF1

def interpolate_popsize_milestones(plotDf):
    '''
    Assesses data being plotted to locate milestones to mark. These milestones correspond to the lowest
    population size necessary to achieve a certain confidence or likelihood that an experiment would have
    at least as much "success" with QTL prediction. as defined by the signal strength.
    
    Interpolation occurs since not all parameter combinations are evenly dispersed over a range, and we would
    like to give the most granular estimate possible.
    
    Returns:
        weak95 -- an integer of the pop. size which first has >=950 samples with at least a weak signal
        mid95 -- an integer of the pop. size which first has >=950 samples with at least a mid/intermediate signal
        strong95 -- an integer of the pop. size which first has >=950 samples with at least a strong signal
    '''
    CUTOFF = 950
    
    # Get interpolated milestone counts
    x = plotDf["pop_size"]
    milestoneInterpolations = []
    for countCol in ["atleast_weak", "atleast_mid", "strong"]:
        y = plotDf[countCol]
        interpolatedX, interpolatedY = [], []
        for i in range(len(x)-1):
            thisX, nextX = x.iloc[i], x.iloc[i+1]
            thisY, nextY = y.iloc[i], y.iloc[i+1]
            slope = (nextY - thisY) / (nextX - thisX)
            
            interpolatedX += [thisX] + [thisX + pos for pos in range(1, 10)] + [nextX]
            interpolatedY += [int(thisY)] + [int(thisY + (slope * pos)) for pos in range(1, 10)] + [int(nextY)]
        milestoneInterpolations.append((interpolatedX, interpolatedY))
        
    # Find the earliest interpolated cutoff that meets a milestone
    milestones = []
    for interpolations in milestoneInterpolations:
        foundMilestone = None
        for xVal, yVal in zip(*interpolations):
            if yVal >= CUTOFF:
                foundMilestone = xVal
                break
        milestones.append(foundMilestone)
    
    return milestones

def main():
    usage = """%(prog)s simulates genotype calls for a population that segregates over
    a single QTL SNP. It assumes a simple scenario with diploid heterozygous parents
    and a recessive trait. It will output a VCF and metadata file suitable for
    psQTL analysis.
    """
    
    # Parse command line arguments
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-t", dest="threads",
                   required=True,
                   type=int,
                   help="Specify the number of threads to run for program speedup")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify the location to write outputs")
    p.add_argument("--snpMB", dest="snpMB",
                   required=False,
                   type=int,
                   help="Specify the number of SNPs per megabase (default==1000)",
                   default=1000)
    p.add_argument("--genomeLength", dest="genomeLength",
                   required=False,
                   type=int,
                   help="Specify the genome length in base pairs (default==10000000)",
                   default=10000000)
    p.add_argument("--centimorgans", dest="cmMB",
                   required=False,
                   type=float,
                   help="Specify the centiMorgan per megabase (default==3.0)",
                   default=3.0)
    p.add_argument("--offspring", dest="numOffspring",
                   required=False,
                   type=int,
                   help="Specify the number of offspring to simulate (default=10000)",
                   default=10000)
    
    args = p.parse_args()
    os.makedirs(args.outputDirectory, exist_ok=True)
    
    ########################
    # SIMULATE POPULATIONS #
    ########################
    
    # Simulate diploid parent genomes
    PARENTS_NPY = os.path.join(args.outputDirectory, "parents.npy")
    if not os.path.isfile(PARENTS_NPY):
        ploidy1 = np.ones(args.genomeLength//args.snpMB).reshape(1, args.genomeLength//args.snpMB) # ones indicate the positive trait; must be homozygous
        ploidy2 = np.zeros(args.genomeLength//args.snpMB).reshape(1, args.genomeLength//args.snpMB) # zeros indicate the negative trait; any presence means no good trait
        parentGenome = np.dstack((ploidy1, ploidy2))
        parentGenome = parentGenome.astype(bool)
        
        # Duplicate the parent to enable crossing and save to file
        parentGenomes = np.vstack((parentGenome, parentGenome))
        np.save(PARENTS_NPY, parentGenomes)
    
    # Generate genetic map
    mapArray = [["CHR.PHYS", "cM", "Trait"]]
    for i in range(0, args.genomeLength//args.snpMB):
        physicalPosition = i*args.snpMB
        cMPosition = (physicalPosition / 1000000) * args.cmMB
        mapArray.append(["chr1", cMPosition, 0.01])
    geneticMap = pd.DataFrame(mapArray[1:], columns=mapArray[0])
    
    # Simulate progeny
    "Run as batches to limit memory consumption"
    GOOD_PICKLE = os.path.join(args.outputDirectory, f"goodF1.pkl")
    if not os.path.isfile(GOOD_PICKLE):
        goodF1 = None
        seed = -1
        while goodF1 is None or len(goodF1) < args.numOffspring:
            seed += 1
            simulatedF1 = simulate_f1(geneticMap, "good", seed, PARENTS_NPY, args.numOffspring, args.genomeLength//args.snpMB)
            if goodF1 is None:
                goodF1 = simulatedF1
            else:
                goodF1 = np.concatenate((goodF1, simulatedF1), 0)
        goodF1 = np.array(goodF1) # stop this being a jaxlib object as it doesn't play nice with the ProcessPoolExecutor
        
        with open(GOOD_PICKLE, "wb") as pickleOut:
            pickle.dump(goodF1, pickleOut)
    else:
        with open(GOOD_PICKLE, "rb") as pickleIn:
            goodF1 = pickle.load(pickleIn)
    
    BAD_PICKLE = os.path.join(args.outputDirectory, f"badF1.pkl")
    if not os.path.isfile(BAD_PICKLE):
        badF1 = None
        seed = -1
        while badF1 is None or len(badF1) < args.numOffspring:
            seed += 1
            simulatedF1 = simulate_f1(geneticMap, "bad", seed, PARENTS_NPY, args.numOffspring, args.genomeLength//args.snpMB)
            if badF1 is None:
                badF1 = simulatedF1
            else:
                badF1 = np.concatenate((badF1, simulatedF1), 0)
        badF1 = np.array(badF1) # stop this being a jaxlib object as it doesn't play nice with the ProcessPoolExecutor
        
        with open(BAD_PICKLE, "wb") as pickleOut:
            pickle.dump(badF1, pickleOut)
    else:
        with open(BAD_PICKLE, "rb") as pickleIn:
            badF1 = pickle.load(pickleIn)
    
    #############################
    #    SIMULATE EXPERIMENT    #
    #############################
    
    # Establish intermediate and output file locations
    OKAY_DIR = os.path.join(args.outputDirectory, "okay_flags")
    os.makedirs(OKAY_DIR, exist_ok=True)
    R2_PLOTS_DIR = os.path.join(args.outputDirectory, "param_plots")
    os.makedirs(R2_PLOTS_DIR, exist_ok=True)
    RESULT_DATA_DIR = os.path.join(args.outputDirectory, "result_data")
    os.makedirs(RESULT_DATA_DIR, exist_ok=True)
    
    # Set combination variables
    populationSize = list(range(10, 200, 10)) + list(range(200, 420, 20))
    populationBalance = [ x/100 for x in range(10, 60, 10) ]
    phenotypeErrorPct = [ x/100 for x in range(0, 55, 5) ]
    errorSymmetry = [ "bulk1", "both", "bulk2" ]
    
    # Obtain evaluation results for all combinations of variables
    productList = list(product(populationSize, populationBalance, phenotypeErrorPct, errorSymmetry))
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        for i, (size, balance, phePct, symmetry) in enumerate(productList):
            paramString = f"{size}_{balance}_{phePct}_{symmetry}"
            
            OKAY_FILE = os.path.join(OKAY_DIR, f"{paramString}.okay")
            PLOT_FILE = os.path.join(R2_PLOTS_DIR, f"{paramString}.png")
            DATA_FILE = os.path.join(RESULT_DATA_DIR, f"{paramString}.txt")
            if os.path.exists(OKAY_FILE) and os.path.exists(PLOT_FILE) and os.path.exists(DATA_FILE):
                continue
            
            # Determine bulk1 (good) size
            numBulk1 = size * balance
            if numBulk1 % 1 != 0: # skip if the bulk size is not a whole number
                continue
            
            # Convert to integer & derive the bulk2 (bad) size
            numBulk1 = int(numBulk1)
            numBulk2 = size - numBulk1
            if (numBulk1 * phePct % 1 != 0) or (numBulk2 * phePct % 1 != 0): # skip if we can't simulate phenotype error with whole numbers
                continue
            
            # Identify the number of samples to have errors in each bulk
            numBulk1Error = int(numBulk1 * phePct)
            numBulk2Error = int(numBulk2 * phePct)
            
            # Parallel processing for ED^4 calculation
            replicatedED4 = []
            for repNum in range(NREP):
                # Randomly select individuals for each bulk with phenotype error modeling
                "We want to pull out individuals here to avoid duplicating the full F1 data across each subprocess"
                if symmetry in ["bulk1", "both"] and numBulk1Error != 0:
                    bulk1 =  np.concatenate((
                        random_sample(goodF1, numBulk1-numBulk1Error),
                        random_sample(badF1, numBulk1Error),
                    ))
                else:
                    bulk1 = random_sample(goodF1, numBulk1)
                
                if symmetry in ["bulk2", "both"] and numBulk2Error != 0:
                    bulk2 =  np.concatenate((
                        random_sample(badF1, numBulk2-numBulk2Error),
                        random_sample(goodF1, numBulk2Error),
                    ))
                else:
                    bulk2 = random_sample(badF1, numBulk2)
                
                # Submit the job to the process executor
                futureResult = executor.submit(experiment_segregation, bulk1, bulk2, args.genomeLength//args.snpMB)
                replicatedED4.append(futureResult)
            
            # Establish container for plotting of R^2 measurement cutoffs
            paramImages = {}
            lastCutoff = -math.inf
            for cutoff in list(range(0, 11, 1)):
                cutoff = cutoff*0.1
                paramImages[cutoff] = [lastCutoff, None]
                lastCutoff = cutoff
            
            # Get R^2 calculation for each replicated simulation result
            replicatedR2 = []
            for result in replicatedED4:
                x = np.array(range(0, args.genomeLength, int(1000000 / args.snpMB))) # 1Mbp divide by num SNPs per Mbp gives the space between each SNP
                y = np.array(result.result())
                triangleY, r_squared = triangle_fit(x, y)
                replicatedR2.append(r_squared)
                
                # Store this as an image for validation if we don't have one in this cutoff range yet
                for upperCutoff, (lowerCutoff, image) in paramImages.items():
                    if (lowerCutoff < r_squared <= upperCutoff) and image == None:
                        plotImg = plot_triangle_fit(x, y, triangleY, r_squared)
                        paramImages[upperCutoff] = [lowerCutoff, plotImg]
            replicatedR2 = np.array(replicatedR2)
            
            # Calculate signal strength values and write to file
            "Writing to file rather than storing in memory allows for program resumption if interrupted"
            noSignal = sum((replicatedR2 < WEAK_RSQ) | (np.isnan(replicatedR2)))
            weakSignal = sum((replicatedR2 >= WEAK_RSQ) & (replicatedR2 < MID_RSQ))
            moderateSignal = sum((replicatedR2 >= MID_RSQ) & (replicatedR2 < STRONG_RSQ))
            strongSignal = sum(replicatedR2 >= STRONG_RSQ)
            with open(DATA_FILE, "w") as fileOut:
                fileOut.write(f"pop_size\tpop_balance\tphenotype_error\tsymmetry\tnone\tweak\tmid\tstrong\n")
                fileOut.write(f"{size}\t{balance}\t{phePct}\t{symmetry}\t{noSignal}\t{weakSignal}\t{moderateSignal}\t{strongSignal}\n")
            
            # Format and output the R^2 measurement cutoffs plot
            concatImage = Image.new("RGB", (IMG_WIDTH*NUM_COLS, IMG_HEIGHT*NUM_ROWS))
            x_offset = 0
            y_offset = 0
            ongoingCount = 0
            for upperCutoff, (lowerCutoff, image) in paramImages.items():
                if image == None:
                    image = Image.new("RGB", (IMG_WIDTH, IMG_HEIGHT)) # this will be a blank/black square
                
                concatImage.paste(image, (x_offset, y_offset))
                x_offset += IMG_WIDTH
                
                ongoingCount += 1
                if ongoingCount % NUM_COLS == 0:
                    x_offset = 0
                    y_offset += IMG_HEIGHT
            concatImage.save(PLOT_FILE)
            
            # Touch an .okay file to denote successful completion
            "This allows for the program to be resumed correctly if interrupted"
            with open(OKAY_FILE, "w") as fileOut:
                pass
    
    ########################
    #    FINAL PLOTTING    #
    ########################
    
    # Establish output file location
    RESULT_PLOTS_DIR = os.path.join(args.outputDirectory, "result_plots")
    os.makedirs(RESULT_PLOTS_DIR, exist_ok=True)
    
    # Load in all evaluation results to a common pandas DataFrame
    summaryDF = None
    for i, (size, balance, phePct, symmetry) in enumerate(productList):
        paramString = f"{size}_{balance}_{phePct}_{symmetry}"
        DATA_FILE = os.path.join(RESULT_DATA_DIR, f"{paramString}.txt")
        if not os.path.isfile(DATA_FILE): # TBD: remove after testing
            continue # TBD: remove after testing
        thisDF = pd.read_csv(DATA_FILE, sep="\t")
        if summaryDF is None:
            summaryDF = thisDF
        else:
            summaryDF = pd.concat([summaryDF, thisDF], axis=0)
    
    # Create plot for depiction of overall parameter results
    pops = sorted(summaryDF["pop_size"].unique())
    balances = sorted(summaryDF["pop_balance"].unique(), reverse=True)
    errors = sorted(summaryDF["phenotype_error"].unique())
    symmetries = sorted(summaryDF["symmetry"].unique())
    cutoffs = []
    for fignum, phenotype_error in enumerate(errors):
        errorDF = summaryDF[(summaryDF["phenotype_error"] == phenotype_error)]
        for symmetry in symmetries:
            symmetryDF = errorDF[(errorDF["symmetry"] == symmetry)]
            
            # Establish the matplotlib figure and axes
            fig, ax = plt.subplots(nrows=5, ncols=1, tight_layout=True)
            fig.supxlabel("Population Size")
            fig.supylabel("Proportion", x=0)
            fig.suptitle("Phenotype Error: {:.0f}%".format(phenotype_error*100))
            
            # Plot the stacked data for each balance interval
            for i, pop_balance in enumerate(balances):
                plotDf = symmetryDF[(symmetryDF["pop_balance"] == pop_balance)]
                if len(plotDf) == 0:
                    continue
                
                # Identify signal milestones
                plotDf["atleast_weak"] = plotDf["weak"] + plotDf["mid"] + plotDf["strong"]
                plotDf["atleast_mid"] = plotDf["mid"] + plotDf["strong"]
                weak95, mid95, strong95 = interpolate_popsize_milestones(plotDf)
                cutoffs.append([phenotype_error, pop_balance, symmetry, weak95, mid95, strong95])
                
                # Stackplot
                ax[i].stackplot(plotDf["pop_size"],
                                (plotDf["strong"], plotDf["mid"], plotDf["weak"], plotDf["none"]),
                                labels=["Strong", "Mid", "Weak", "None"])
                ax[i].margins(0,0) # fill all of the rectangular space
                
                # Mark the 95% cutoffs
                if weak95:
                    ax[i].vlines(weak95, ymin=0, ymax=1000, linestyle="--", color="white")
                if mid95:
                    ax[i].vlines(mid95, ymin=0, ymax=1000, linestyle="--", color="grey")
                if strong95:
                    ax[i].vlines(strong95, ymin=0, ymax=1000, linestyle="--", color="black")
                
                # Set y axis aesthetics
                ax[i].set_yticks([0, 250, 500, 750, 1000])
                
                ax[i].yaxis.set_minor_locator(ticker.FixedLocator([0, 250, 500, 750, 1000]))
                ax[i].yaxis.set_major_locator(ticker.FixedLocator([0, 500, 1000]))
                ax[i].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{x/10:.0f}%"))
                
                # Set x axis aesthetics
                ax[i].set_xlim(10, 400)
                if i < (len(balances) - 1):
                    ax[i].set_xticklabels([]) # only show tick labels on final (bottom) plot
            
            #handles, labels = ax[0].get_legend_handles_labels() # turn this on to manually obtain a legend to reformat in PDF editor
            #fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1, 0.5))
            
            # Save the figure
            outputFileName = os.path.join(RESULT_PLOTS_DIR, f"{fignum}.{symmetry}.pdf")
            fig.savefig(outputFileName, dpi=300, pad_inches=1.5)
    
    # Create a plot to depict the trend in cutoffs with varying error and pop balance
    errors = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5] # only showing 0.1 increments
    
    df = pd.DataFrame(cutoffs) # format the data stored during results plotting
    df.columns = ["phenotype_error", "pop_balance", "symmetry", "weak95", "mid95", "strong95"]
    
    fig = plt.figure()
    ax = fig.add_subplot()
    for i, balance in enumerate(balances):
        colour = LINE_COLOURS[i]
        ax.plot([], [], color=colour, label=balance) # dummy line for legend
        
        data = df[(df["symmetry"] == "both") & (df["pop_balance"] == balance) & (df["phenotype_error"].isin(errors))]
        for j, cutoff in enumerate(["weak95", "mid95", "strong95"]):
            style = LINE_STYLES[j]
            lineData = data[pd.notna(data[cutoff])]
            x = lineData[cutoff]
            y = lineData["phenotype_error"] * 100
            
            ax.plot(x, y, linestyle=style, color=colour)
            if i == 0:
                ax.plot([], [], linestyle=style, label=cutoff) # dummy line for legend
    
    # Handle axis labels and aesthetics
    ax.invert_yaxis() # higher Y values at the top
    ax.set_xlabel("Population Size")
    ax.set_ylabel("Phenotype Error (%)")
    
    # Handle legend plotting
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1, 0.5))
    
    # Save the figure
    outputFileName = os.path.join(RESULT_PLOTS_DIR, f"trends.pdf")
    fig.savefig(outputFileName, dpi=300, pad_inches=1.5)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
