import os, shutil, subprocess
import concurrent.futures
from itertools import repeat

# Validation functions
def validate_bcftools_exists():
    '''
    Will check if 'bcftools' is in the PATH and raise an exception if it isn't.
    '''
    if not shutil.which("bcftools"):
        raise FileNotFoundError("bcftools not found in PATH")

# Threaded operations
def call_task(bamListFile, genomeFasta, contigID, outputFileName):
    '''
    Partner function for run_bcftools_call. Will pipeline bcftools mpileup->call to
    call variants on a single BAM file.
    
    Parameters:
        bamListFile -- a string pointing to a file containing a list of BAM files to process.
        genomeFasta -- a string pointing to the genome FASTA that the BAM files were aligned to.
        contigID -- the contig ID to process.
        outputFileName -- the output file name.
    '''
    # Format command
    cmd = ["bcftools", "mpileup", "-Ou", "-f", genomeFasta,
           "-r", contigID, "--bam-list", bamListFile,
           "-q", "10", "-Q", "20", "-a", "AD",
           "|",
           "bcftools", "call", "-m", "-v",
           "-Oz", "-o", outputFileName]
    
    # Run bcftools mpileup->call
    run_call = subprocess.Popen(" ".join(cmd), shell=True,
                                 stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE)
    callout, callerr = run_call.communicate()
    
    # Check for errors
    if run_call.returncode == 0:
        open(outputFileName + ".ok", "w").close() # touch a .ok file to indicate success
        return None
    else:
        errorMsg = callerr.decode("utf-8").rstrip("\r\n ")
        raise Exception(("run_bcftools_call encountered an unhandled situation; have a " + 
                         f"look at the stderr '{errorMsg}' to make sense of this."))

def run_bcftools_call(bamListFile, genomeFasta, outputDirectory, threads):
    '''
    Will run samtools depth on a list of BAM files in parallel.
    
    Parameters:
        bamListFile -- a string pointing to a file containing a list of BAM files to process.
        genomeFasta -- a string pointing to the genome FASTA that the BAM files were aligned to.
        outputDirectory -- a string pointing to the directory to write output files to.
        threads -- an integer indicating how many threads to run GMAP with.
    '''
    # Get contig IDs from the genome FASTA
    contigIDs = set()
    with open(genomeFasta, "r") as fastaFile:
        for line in fastaFile:
            if line.startswith(">"):
                contigID = line[1:].split(" ")[0].rstrip("\r\n")
                if contigID not in contigIDs:
                    contigIDs.add(contigID)
                else:
                    raise ValueError(f"Contig ID '{contigID}' is duplicated in the genome FASTA")
    
    # Drop any contig IDs that have already been processed
    contigsToProcess = []
    outputFileNames = []
    for contigID in contigIDs:
        outputFileName = os.path.join(outputDirectory, f"{contigID}.vcf.gz")
        if not os.path.exists(outputFileName):
            contigsToProcess.append(contigID)
            outputFileNames.append(outputFileName)
        else:
            print(f"# '{contigID}' has already had variants called on it; skipping...")
    
    # Skip if there's nothing to do
    if len(contigsToProcess) == 0:
        print("# No contigs to call variants on; skipping...")
        return
    
    # Plug data into the threaded function
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        executor.map(
            call_task,
            repeat(bamListFile, len(contigsToProcess)),
            repeat(genomeFasta, len(contigsToProcess)),
            contigsToProcess,
            outputFileNames
        )
