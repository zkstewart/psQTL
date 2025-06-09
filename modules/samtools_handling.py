import os, shutil, subprocess, math
import concurrent.futures
import numpy as np
from itertools import repeat

from .parsing import parse_samtools_depth_tsv

# Validation functions
def validate_samtools_exists():
    '''
    Will check if 'samtools' is in the PATH and raise an exception if it isn't.
    '''
    if not shutil.which("samtools"):
        raise FileNotFoundError("samtools not found in PATH")

# Single threaded operations
def run_samtools_faidx(fastaFile):
    '''
    Indexes a FASTA file using samtools faidx.
    
    Parameters:
        fastaFile -- the input FASTA file.
    '''
    # Format command
    cmd = ["samtools", "faidx", fastaFile]
    
    # Run samtools depth
    run_faidx = subprocess.Popen(" ".join(cmd), shell=True,
                                 stdout = subprocess.DEVNULL,
                                 stderr = subprocess.PIPE)
    faidxout, faidxerr = run_faidx.communicate()
    
    # Check for errors
    if run_faidx.returncode == 0:
        return None
    else:
        errorMsg = faidxerr.decode("utf-8").rstrip("\r\n ")
        raise Exception(("run_samtools_faidx encountered an unhandled situation when processing " + 
                         f"'{fastaFile}'; have a look at the stderr to make sense of this:\n'{errorMsg}'"))

# Threaded operations
def depth_task(ioPair):
    '''
    Partner function for run_samtools_depth. Will run samtools depth on a single BAM file.
    
    Parameters:
        ioPair -- a list of lists containing the input BAM file and the output file name.
    '''
    inputBamFile, outputFileName = ioPair
    
    # Format command
    cmd = ["samtools", "depth", "-a", "-q", "13", inputBamFile]
    
    # Run samtools depth
    with open(outputFileName, "w") as fileOut:
        run_depth = subprocess.Popen(" ".join(cmd), shell=True,
                                     stdout = fileOut,
                                     stderr = subprocess.PIPE)
        depthout, deptherr = run_depth.communicate()
    
    # Check for errors
    if run_depth.returncode == 0:
        open(outputFileName + ".ok", "w").close() # touch a .ok file to indicate success
        return None
    else:
        errorMsg = deptherr.decode("utf-8").rstrip("\r\n ")
        raise Exception(errorMsg)

def run_samtools_depth(ioList, threads):
    '''
    Will run samtools depth on a list of BAM files in parallel.
    
    Parameters:
        ioList -- a list of lists containing paired strings for the input BAM
                  file and the output file name.
        threads -- an integer indicating how many threads to run.
    '''
    if len(ioList) == 0:
        print("# All depth files already exist; skipping...")
        return
    else:
        print(f"# Generating depth files for {len(ioList)} BAM file{'s' if len(ioList) > 1 else ''} ...")
    
    futures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for ioPair in ioList:
            futures.append(executor.submit(depth_task, ioPair))
    for f in futures:
        try:
            result = f.result()
        except Exception as e:
            raise Exception(f"run_samtools_depth encountered the following error:\n{e}")

##

def depth_to_histoDict(depthFile, lengthsDict, binSize):
    '''
    Parameters:
        depthFile -- a string containing the path to the input TSV file output by samtools depth.
        lengthsDict -- a dictionary pairing contig IDs (keys) with their lengths (values).
        binSize -- an integer indicating the size of the bin to sum depth values within.
    Returns:
        histoDict -- a dictionary pairing contig IDs (keys) with a numpy array of binned depth values.
    '''
    histoDict = {}
    contigZeroBase = {}
    for contigID, pos, depth in parse_samtools_depth_tsv(depthFile):
        if not contigID in histoDict:
            # Establish storage data structure for this contig
            histoDict[contigID] = np.array([ 0
                for windowChunk in range(math.ceil(lengthsDict[contigID] / binSize))
            ])
            
            # Check if the contig is 0-based or 1-based
            if pos == 0:
                contigZeroBase[contigID] = 0
            elif pos == 1:
                contigZeroBase[contigID] = 1 # results in subtracting 1 from the position
            else:
                raise ValueError(f"First line in contig {contigID} is not 0 or 1; please check " +
                                 "that the input depth file is valid and produced by samtools depth.")
        
        binIndex = (pos-contigZeroBase[contigID]) // binSize # adjust position to be 0-based if necessary
        histoDict[contigID][binIndex] += depth
    return histoDict

def bin_task(ioPair, lengthsDict, binSize):
    '''
    Partner function for bin_samtools_depth. Will receive samtools depth outputs and bin them
    into a new TSV file.
    
    Parameters:
        ioPair -- a list of lists containing the input depth file and the output file name.
        lengthsDict -- a dictionary pairing contig IDs (keys) with their lengths (values).
    '''
    depthFile, outputFileName = ioPair
    if binSize < 1:
        raise ValueError("binSize must be a positive integer")
    
    # Parse the input depth file
    histoDict = depth_to_histoDict(depthFile, lengthsDict, binSize)
    
    # Write the binned depth values to a file
    with open(outputFileName, "w") as fileOut:
        for contigID, depthList in histoDict.items():
            for binIndex, depthValue in enumerate(depthList):
                fileOut.write(f"{contigID}\t{binIndex * binSize}\t{depthValue}\n")
    open(outputFileName + ".ok", "w").close() # touch a .ok file to indicate success

def bin_samtools_depth(ioList, lengthsDict, threads, binSize):
    '''
    Will receive TSV files output by samtools depth to bin them into new TSV files in parallel.
    
    Parameters:
        ioList -- a list of lists containing paired strings for the input depth
                  file and the output file name.
        lengthsDict -- a dictionary pairing contig IDs (keys) with their lengths (values).
        threads -- an integer indicating how many threads to run.
        binSize --  an integer indicating the size of the bin to sum depth values within.
    '''
    if len(ioList) == 0:
        print("# All depth files have already been binned; skipping...")
        return
    else:
        print(f"# Binning {len(ioList)} depth file{'s' if len(ioList) > 1 else ''} ...")
    
    futures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for ioPair in ioList:
            futures.append(executor.submit(bin_task, ioPair, lengthsDict, binSize))
    for f in futures:
        try:
            result = f.result()
        except Exception as e:
            raise Exception(f"bin_samtools_depth encountered the following error:\n{e}")
