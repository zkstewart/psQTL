import os, shutil, subprocess
import concurrent.futures

# Validation functions
def validate_samtools_exists():
    '''
    Will check if 'samtools' is in the PATH and raise an exception if it isn't.
    '''
    if not shutil.which("samtools"):
        raise FileNotFoundError("samtools not found in PATH")

# Single threaded operations
def run_samtools_index(fastaFile):
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
        raise Exception(("run_samtools_index encountered an unhandled situation when processing " + 
                         f"'{fastaFile}'; have a look at the stderr '{errorMsg}' " +
                         "to make sense of this."))

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
        raise Exception(("run_samtools_depth encountered an unhandled situation when processing " + 
                         f"'{inputBamFile}'; have a look at the stderr '{errorMsg}' " +
                         "to make sense of this."))

def run_samtools_depth(ioList, threads):
    '''
    Will run samtools depth on a list of BAM files in parallel.
    
    Parameters:
        ioList -- a list of lists containing paired strings for the input BAM
                  file and the output file name.
        threads -- an integer indicating how many threads to run GMAP with.
    '''
    if len(ioList) == 0:
        print("# All depth files already exist; skipping...")
    else:
        print(f"# Generating depth files for {len(ioList)} unprocessed BAM file{'s' if len(ioList) > 1 else ''} ...")
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        executor.map(depth_task, ioList)
