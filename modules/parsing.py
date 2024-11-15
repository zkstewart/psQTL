import os, json

# Parameter cache functions
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
    CACHEABLE_PARAMS = ["workingDirectory", "metadataFile", "vcfFile", "bamFiles", "depthFiles"]
    
    # Detect any existing param cache file
    paramsDict = load_param_cache(args.workingDirectory)
    if paramsDict != {}:
        raise FileExistsError("Directory has already been initialised; delete this folder or use a new one.")
    
    # Initialise param cache dictionary
    paramsDict = {
        param : args.__dict__[param]
        for param in CACHEABLE_PARAMS
    }
    
    # Write updated param cache to file
    with open(os.path.join(args.workingDirectory, "param_cache.json"), "w") as fileOut:
        json.dump(paramsDict, fileOut)

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
    ALLOWED_PARAMS = ["bamFiles", "vcf", "depth", "binSize", "binStep"]
    
    # Detect any existing param cache file
    paramsDict = load_param_cache(workingDirectory)
    if paramsDict == {}:
        raise FileExistsError("Directory has not been initialised yet.")
    
    # Update param cache dictionary
    for key, value in updatesDict.items():
        if key not in ALLOWED_PARAMS:
            raise KeyError(f"Parameter '{key}' not allowed to be updated!")
        paramsDict[key] = value
    
    # Write updated param cache to file
    with open(os.path.join(workingDirectory, "param_cache.json"), "w") as fileOut:
        json.dump(paramsDict, fileOut)

# VCF cache functions
def initialise_vcf_cache(workingDirectory):
    '''
    Initialises a new parameter cache file in the output directory, storing the
    parameters used for the current program run.
    
    Parameters:
        args -- the argparse object generated through the psQTL_prep 'initialise' submodule.
    '''
    CACHEABLE_PARAMS = ["workingDirectory", "metadataFile", "vcfFile", "bamFiles", "depthFiles"]
    
    # Detect any existing VCF cache file
    vcfDict = load_vcf_cache(workingDirectory)
    if vcfDict != {}:
        raise FileExistsError("VCF cache has already been initialised; no need to initialise a new one.")
    
    # Detect any existing param cache file
    paramsDict = load_param_cache(workingDirectory)
    if paramsDict == {}:
        raise FileExistsError("Directory has not been initialised yet.")
    
    # Parse the VCF file for cacheable metadata values
    vcfFile = paramsDict["vcfFile"]
    if vcfFile == None:
        raise ValueError("VCF file not specified in params cache; cannot initialise VCF cache.")
    
    numVariants = 0
    contigs = {} # using as an ordered set
    with open(vcfFile, "r") as fileIn:
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
