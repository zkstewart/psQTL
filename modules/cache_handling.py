import os, json, gzip
from contextlib import contextmanager

from .parsing import open_gz_file, parse_metadata

CACHEABLE_PARAMS = ["workingDirectory", "metadataFile",
                    "vcfFile", "deletionFile",
                    "bamFiles", "bamSuffix"]

# Parameter cache functions
def merge_cache_into_args(args):
    '''
    Takes the values in the parameter cache and merges them into the argparse object
    if the argparse object has empty values for the cacheable parameters. In practice,
    this allows the 'initialise' submodule of psQTL_prep to be run, with subsequent
    submodules using the cache values if the user does not specify them in the command
    line.
    '''
    paramsDict = load_param_cache(args.workingDirectory)
    for param in CACHEABLE_PARAMS:
        if param in paramsDict and param in args.__dict__: # only merge if the parameter is both the cache and the args
            if args.__dict__[param] == None or args.__dict__[param] == []: # only use existing cache values if args are empty
                if paramsDict[param] != None or paramsDict[param] != []: # only use cache values if they are _not_ empty
                    args.__dict__[param] = paramsDict[param]
        elif param in paramsDict: # if the parameter is in the cache but not in the args, add it to the args
            args.__dict__[param] = paramsDict[param]

def load_param_cache(workingDirectory):
    '''
    Loads the parameter cache file from the working directory, if it exists,
    as a dictionary with JSON parsing.
    
    Parameters:
        workingDirectory -- a string indicating the parent dir where the analysis is being
                            run.
    '''
    # Parse any existing param cache file
    paramCacheFile = os.path.join(workingDirectory, "param_cache.json")
    if os.path.exists(paramCacheFile):
        try:
            with open(paramCacheFile, "r") as fileIn:
                paramsDict = json.load(fileIn)
        except:
            raise Exception((f"'{paramCacheFile}' exists but cannot be loaded as a JSON. " + 
                             "If the file is malformed, delete it and re-initialise this directory."))
    else:
        paramsDict = {}
    
    return paramsDict

def initialise_param_cache(args):
    '''
    Initialises a new parameter cache file in the output directory, storing the
    parameters used for the current program run.
    
    Parameters:
        args -- the argparse object generated through the psQTL_prep 'initialise' submodule.
    '''
    # Detect any existing param cache file
    paramsDict = load_param_cache(args.workingDirectory)
    
    # Initialise new param cache dictionary values
    newParamsDict = {
        param : args.__dict__[param]
        for param in CACHEABLE_PARAMS
    }
    
    # Write updated param cache to file
    if paramsDict != newParamsDict:
        with open(os.path.join(args.workingDirectory, "param_cache.json"), "w") as fileOut:
            json.dump(newParamsDict, fileOut)
    else:
        print("# Parameter cache already exists and is up-to-date; skipping ...")

def update_param_cache(workingDirectory, updatesDict):
    '''
    Updates a parameter cache file in the output directory, modifying only parameters
    that are allowed to change.
    
    Parameters:
        workingDirectory -- a string indicating the parent dir where the analysis is being
                            run.
        updatesDict -- a dictionary containing the parameters to update, with the keys
                       being the parameter names and the values being the new values.
    '''
    ALLOWED_PARAMS = ["metadataFile", "bamFiles", "bamSuffix", "vcfFile", "deletionFile"]
    
    # Detect any existing param cache file
    paramsDict = load_param_cache(workingDirectory)
    if paramsDict == {}:
        raise FileExistsError("Directory has not been initialised yet.")
    
    # Update param cache dictionary
    madeChanges = False
    for key, value in updatesDict.items():
        if key not in ALLOWED_PARAMS:
            raise KeyError(f"Parameter '{key}' not allowed to be updated!")
        else:
            if paramsDict[key] == value: # skip if the value is unchanged
                pass
            else:
                paramsDict[key] = value
                print(f"# Updated '{key}' in parameters cache.")
                madeChanges = True
    
    # Write updated param cache to file
    if madeChanges:
        with open(os.path.join(workingDirectory, "param_cache.json"), "w") as fileOut:
            json.dump(paramsDict, fileOut)

# VCF cache functions
def initialise_vcf_cache(workingDirectory):
    '''
    Initialises a new VCF cache file in the output directory, storing details pertaining
    to the VCF file used in the current program run.
    
    Parameters:
        args -- the argparse object generated through the psQTL_prep 'initialise' submodule.
    '''
    # Detect any existing param cache file
    paramsDict = load_param_cache(workingDirectory)
    if paramsDict == {}:
        raise FileExistsError("Directory has not been initialised yet.")
    
    # Detect any existing VCF cache file
    vcfDict = load_vcf_cache(workingDirectory)
    if vcfDict != {}:
        if vcfDict["vcfFile"] == paramsDict["vcfFile"]:
            print("# VCF cache already exists and is up-to-date; skipping ...")
            return
    
    # Parse the VCF file for cacheable metadata values
    vcfFile = paramsDict["vcfFile"]
    if vcfFile == None:
        raise ValueError("VCF file not specified in params cache; cannot initialise VCF cache.")
    
    numVariants = 0
    contigs = {} # using as an ordered set
    with open_gz_file(vcfFile) as fileIn:
        for line in fileIn:
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
            elif line.startswith("#"):
                continue
            else:
                contig = line.split("\t")[0]
                if contig not in contigs:
                    contigs[contig] = None
                numVariants += 1
    
    # Initialise VCF cache dictionary
    vcfDict = {
        "vcfFile" : vcfFile,
        "variants" : numVariants,
        "samples" : samples,
        "contigs" : list(contigs.keys())
    }
    
    # Write VCF cache to file
    with open(os.path.join(workingDirectory, "vcf_cache.json"), "w") as fileOut:
        json.dump(vcfDict, fileOut)

def load_vcf_cache(workingDirectory):
    '''
    Loads the VCF cache file from the working directory, if it exists,
    as a dictionary with JSON parsing.
    
    Parameters:
        workingDirectory -- a string indicating the parent dir where the analysis is being
                            run.
    '''
    # Parse any existing param cache file
    vcfCacheFile = os.path.join(workingDirectory, "vcf_cache.json")
    if os.path.exists(vcfCacheFile):
        try:
            with open(vcfCacheFile, "r") as fileIn:
                vcfDict = json.load(fileIn)
        except:
            raise Exception((f"'{vcfCacheFile}' exists but cannot be loaded as a JSON. " + 
                             "If the file is malformed, delete it and re-initialise this directory."))
    else:
        vcfDict = {}
    
    return vcfDict

# Deletion cache functions
def initialise_deletion_cache(workingDirectory, windowSize):
    '''
    Initialises a new deletion cache file in the output directory, storing the
    parameters used for the current program run.
    
    Parameters:
        args -- the argparse object generated through the psQTL_prep 'initialise' submodule.
        windowSize -- the size of the window to use for depth binning.
    '''
    # Detect any existing param cache file
    paramsDict = load_param_cache(workingDirectory)
    if paramsDict == {}:
        raise FileExistsError("Directory has not been initialised yet.")
    
    # Detect any existing deletion cache file
    deletionDict = load_deletion_cache(workingDirectory)
    if deletionDict != {}:
        if deletionDict["deletionFile"] == paramsDict["deletionFile"]:
            print("# Deletion cache already exists and is up-to-date; skipping ...")
            return
    
    # Parse the deletion file for cacheable metadata values
    deletionFile = paramsDict["deletionFile"]
    if deletionFile == None:
        raise ValueError("Deletion file not specified in params cache; cannot initialise deletion cache.")
    
    totalBins = 0
    deletionBins = 0
    contigs = {} # using as an ordered set
    with open_gz_file(deletionFile) as fileIn:
        for line in fileIn:
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
            elif line.startswith("#"):
                continue
            else:
                sl = line.split("\t")
                if sl[0] not in contigs:
                    contigs[sl[0]] = None
                
                # Tally bins
                deletionBins += 1 if any([ "1" in x for x in sl[9:]]) else 0
                totalBins += 1
    
    # Initialise deletion cache dictionary
    deletionDict = {
        "deletionFile" : deletionFile,
        "totalBins" : totalBins,
        "deletionBins" : deletionBins,
        "samples" : samples,
        "contigs" : list(contigs.keys()),
        "windowSize": windowSize
    }
    
    # Write deletion cache to file
    with open(os.path.join(workingDirectory, "deletion_cache.json"), "w") as fileOut:
        json.dump(deletionDict, fileOut)

def load_deletion_cache(workingDirectory):
    '''
    Loads the deletion cache file from the working directory, if it exists,
    as a dictionary with JSON parsing.
    
    Parameters:
        workingDirectory -- a string indicating the parent dir where the analysis is being
                            run.
    '''
    # Parse any existing param cache file
    deletionCacheFile = os.path.join(workingDirectory, "deletion_cache.json")
    if os.path.exists(deletionCacheFile):
        try:
            with open(deletionCacheFile, "r") as fileIn:
                deletionDict = json.load(fileIn)
        except:
            raise Exception((f"'{deletionCacheFile}' exists but cannot be loaded as a JSON. " + 
                             "If the file is malformed, delete it and re-initialise this directory."))
    else:
        deletionDict = {}
    
    return deletionDict

# Metadata cache functions
def initialise_metadata_cache(workingDirectory):
    '''
    Initialises a new metadata cache file in the output directory, storing the
    values for presentation to a user.
    
    Parameters:
        args -- the argparse object generated through the psQTL_prep 'initialise' submodule.
    '''
    # Detect any existing param cache file
    paramsDict = load_param_cache(workingDirectory)
    if paramsDict == {}:
        raise FileExistsError("Directory has not been initialised yet.")
    
    # Detect any existing metadata cache file
    metadataDict = load_metadata_cache(workingDirectory)
    if metadataDict != {}:
        if metadataDict["metadataFile"] == metadataDict["metadataFile"]:
            print("# Metadata cache already exists and is up-to-date; skipping ...")
            return
    
    # Parse the deletion file for cacheable metadata values
    metadataFile = paramsDict["metadataFile"]
    if metadataFile == None:
        raise ValueError("Metadata file not specified in params cache; cannot initialise metadata cache.")
    metadataDict = parse_metadata(metadataFile)
    
    # Initialise metadata cache dictionary
    initDict = {
        "metadataFile" : metadataFile,
        "bulk1" : list(metadataDict["bulk1"]),
        "bulk2" : list(metadataDict["bulk2"])
    }
    
    # Write VCF cache to file
    with open(os.path.join(workingDirectory, "metadata_cache.json"), "w") as fileOut:
        json.dump(initDict, fileOut)

def load_metadata_cache(workingDirectory):
    '''
    Loads the metadata cache file from the working directory, if it exists,
    as a dictionary with JSON parsing.
    
    Parameters:
        workingDirectory -- a string indicating the parent dir where the analysis is being
                            run.
    '''
    # Parse any existing param cache file
    metadataCacheFile = os.path.join(workingDirectory, "metadata_cache.json")
    if os.path.exists(metadataCacheFile):
        try:
            with open(metadataCacheFile, "r") as fileIn:
                metadataDict = json.load(fileIn)
        except:
            raise Exception((f"'{metadataCacheFile}' exists but cannot be loaded as a JSON. " + 
                             "If the file is malformed, delete it and re-initialise this directory."))
    else:
        metadataDict = {}
    
    return metadataDict
