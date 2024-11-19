import os, shutil, subprocess
import concurrent.futures
from itertools import repeat

from .parsing import read_gz_file

# Validation functions
def validate_bcftools_exists():
    if not shutil.which("bcftools"):
        raise FileNotFoundError("bcftools not found in PATH")

def validate_bgzip_exists():
    if not shutil.which("bgzip"):
        raise FileNotFoundError("bgzip not found in PATH")

def validate_vt_exists():
    if not shutil.which("vt"):
        raise FileNotFoundError("vt not found in PATH")

# Helper functions
def get_contig_ids(genomeFasta):
    '''
    Parses a genome FASTA file and returns a set of contig IDs.
    
    Parameters:
        genomeFasta -- a string pointing to the genome FASTA file.
    Returns:
        contigIDs -- a set of contig IDs.
    '''
    # Get contig IDs from the genome FASTA
    contigIDs = set()
    with read_gz_file(genomeFasta) as fastaFile:
        for line in fastaFile:
            if line.startswith(">"):
                contigID = line[1:].split(" ")[0].rstrip("\r\n")
                if contigID not in contigIDs:
                    contigIDs.add(contigID)
                else:
                    raise ValueError(f"Contig ID '{contigID}' is duplicated in the genome FASTA")
    return contigIDs

# Single threaded operations
def run_bcftools_index(bamFile):
    '''
    Indexes a BAM file using bcftools index.
    
    Parameters:
        bamFile -- the input BAM file.
    '''
    # Format command
    cmd = ["bcftools", "index", bamFile]
    
    # Run samtools depth
    run_index = subprocess.Popen(" ".join(cmd), shell=True,
                                 stdout = subprocess.DEVNULL,
                                 stderr = subprocess.PIPE)
    indexout, indexerr = run_index.communicate()
    
    # Check for errors
    if run_index.returncode == 0:
        return None
    else:
        errorMsg = indexerr.decode("utf-8").rstrip("\r\n ")
        raise Exception(("run_bcftools_index encountered an unhandled situation when processing " + 
                         f"'{bamFile}'; have a look at the stderr '{errorMsg}' " +
                         "to make sense of this."))

def run_bcftools_concat(genomeFasta, workingDirectory, outputFileName):
    '''
    Concatenates multiple VCF files using bcftools concat.
    
    Parameters:
        genomeFasta -- a string pointing to the genome FASTA that the BAM files were aligned to.
        workingDirectory -- a string pointing to the directory to read/write files to/from.
        outputFileName -- a string indicating the output file name to write to.
    '''
    # Get contig IDs from the genome FASTA
    contigIDs = get_contig_ids(genomeFasta)
    
    # Derive VCF file names for concatenation
    vcfFileNames = []
    for contigID in contigIDs:
        vcfFileName = os.path.join(workingDirectory, f"{contigID}.decomposed.vcf.gz")
        
        if not os.path.exists(vcfFileName + ".csi"):
            raise FileNotFoundError(f"Could not find a .csi index for '{vcfFileName}'!")
        
        if os.path.exists(vcfFileName + ".ok"):
            vcfFileNames.append(vcfFileName)
        else:
            raise FileNotFoundError(f"Could not find an OK '{vcfFileName}' file to concatenate!")
    
    # Format command
    cmd = ["bcftools", "concat", "-Oz", "-o", outputFileName, *vcfFileNames]
    
    # Run samtools depth
    print(f"# Concatenating {len(vcfFileNames)} VCF file{'s' if len(vcfFileNames) > 1 else ''} ...")
    run_concat = subprocess.Popen(" ".join(cmd), shell=True,
                                  stdout = subprocess.DEVNULL,
                                  stderr = subprocess.PIPE)
    concatout, concaterr = run_concat.communicate()
    
    # Check for errors
    if run_concat.returncode == 0:
        open(outputFileName + ".ok", "w").close() # touch a .ok file to indicate success
        return None
    else:
        errorMsg = concaterr.decode("utf-8").rstrip("\r\n ")
        raise Exception(("run_bcftools_concat encountered an unhandled situation; " + 
                         f"have a look at the stderr '{errorMsg}' " +
                         "to make sense of this."))

# Threaded operations
def call_task(bamListFile, genomeFasta, contigID, outputFileName):
    '''
    Partner function for run_bcftools_call. Will pipeline bcftools mpileup->call to
    call variants on a single BAM file.
    
    Parameters:
        bamListFile -- a string pointing to a file containing a list of BAM files to process.
        genomeFasta -- a string pointing to the genome FASTA that the BAM files were aligned to.
        contigID -- the contig ID to process.
        outputFileName -- a string indicating the output file name to write to.
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
    Will run bcftools mpileup->call on a list of BAM files in parallel.
    
    Parameters:
        bamListFile -- a string pointing to a file containing a list of BAM files to process.
        genomeFasta -- a string pointing to the genome FASTA that the BAM files were aligned to.
        outputDirectory -- a string pointing to the directory to write output files to.
        threads -- an integer indicating how many threads to run GMAP with.
    '''
    # Get contig IDs from the genome FASTA
    contigIDs = get_contig_ids(genomeFasta)
    
    # Drop any contig IDs that have already been processed
    contigsToProcess = []
    outputFileNames = []
    skipMessages = []
    for contigID in contigIDs:
        outputFileName = os.path.join(outputDirectory, f"{contigID}.vcf.gz")
        if not os.path.exists(outputFileName + ".ok"):
            contigsToProcess.append(contigID)
            outputFileNames.append(outputFileName)
        else:
            skipMessages.append(f"# '{contigID}' has already had variants called on it; skipping...")
    
    # Skip if there's nothing to do
    if len(contigsToProcess) == 0:
        print("# No contigs to call variants on; skipping...")
        return
    else:
        for skipMessage in skipMessages:
            print(skipMessage)
    
    # Plug data into the threaded function
    print(f"# Calling variants on {len(contigsToProcess)} contig{'s' if len(contigsToProcess) > 1 else ''} ...")
    futures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for contigID, outputFileName in zip(contigsToProcess, outputFileNames):
            futures.append(executor.submit(call_task, bamListFile, genomeFasta, contigID, outputFileName))
    for f in futures:
        try:
            result = f.result()
        except Exception as e:
            raise Exception(f"run_bcftools_call encountered the following error: {e}")

##

def normalisation_task(inputFileName, outputFileName, genomeFasta):
    '''
    Partner function for run_bcftools_call. Will pipeline bcftools mpileup->call to
    call variants on a single BAM file.
    
    Parameters:
        inputFileName -- a string pointing to a VCF file to normalise.
        outputFileName -- a string indicating the output file name to write to.
        genomeFasta -- a string pointing to the genome FASTA that the BAM files were aligned to.
    '''    
    # Format command
    cmd = ["bcftools", "norm", "-m-", "-Oz", "-N", inputFileName,
           "|",
           "bcftools", "norm", "-m+", "-Oz", "-N",
           "|",
           "bcftools", "norm", "-f", genomeFasta,
           "|",
           "vt", "decompose_blocksub", "-o", outputFileName, "-"]
    
    # Run bcftools mpileup->call
    run_norm = subprocess.Popen(" ".join(cmd), shell=True,
                                 stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE)
    normout, normerr = run_norm.communicate()
    
    # Check for errors
    if run_norm.returncode == 0:
        open(outputFileName + ".ok", "w").close() # touch a .ok file to indicate success
        return None
    else:
        errorMsg = normerr.decode("utf-8").rstrip("\r\n ")
        raise Exception(("run_normalisation encountered an unhandled situation; have a " + 
                         f"look at the stderr '{errorMsg}' to make sense of this."))

def run_normalisation(genomeFasta, workingDirectory, threads):
    '''
    Will run bcftools/vt normalisation procedure on a list of BAM files in parallel.
    
    Parameters:
        genomeFasta -- a string pointing to the genome FASTA that the BAM files were aligned to.
        workingDirectory -- a string pointing to the directory to read/write files to/from.
        threads -- an integer indicating how many threads to run GMAP with.
    '''
    # Get contig IDs from the genome FASTA
    contigIDs = get_contig_ids(genomeFasta)
    
    # Drop any contig IDs that have already been processed
    inputFileNames = []
    outputFileNames = []
    skipMessages = []
    for contigID in contigIDs:
        inputFileName = os.path.join(workingDirectory, f"{contigID}.vcf.gz")
        if not os.path.exists(inputFileName + ".ok"):
            raise FileNotFoundError(f"Could not find an OK '{inputFileName}' file to normalise!")
        
        outputFileName = os.path.join(workingDirectory, f"{contigID}.decomposed.vcf.gz")
        if not os.path.exists(outputFileName + ".ok"):
            inputFileNames.append(inputFileName)
            outputFileNames.append(outputFileName)
        else:
            skipMessages.append(f"# '{contigID}' has already been normalised; skipping...")
    
    # Skip if there's nothing to do
    if len(inputFileNames) == 0:
        print("# No contigs to normalise; skipping...")
        return
    else:
        for skipMessage in skipMessages:
            print(skipMessage)
    
    # Plug data into the threaded function
    print(f"# Normalising variants for {len(inputFileNames)} VCF{'s' if len(inputFileNames) > 1 else ''} ...")
    futures = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for inputFileName, outputFileName in zip(inputFileNames, outputFileNames):
            futures.append(executor.submit(normalisation_task, inputFileName, outputFileName, genomeFasta))
    for f in futures:
        try:
            result = f.result()
        except Exception as e:
            raise Exception(f"run_normalisation encountered the following error: {e}")
