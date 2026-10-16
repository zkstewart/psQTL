import os, sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from locations import Locations
from parameters import ParameterCache, VcfCache, DepthCache, MetadataCache
from runs import list_runs, find_matching_run, find_legacy_results
from parsing import parse_metadata

def is_complete(fileName):
    return os.path.isfile(fileName) and os.path.isfile(fileName + ".ok")

def view_inputs(locations, args):
    '''
    Presents the files and caches established by psQTL_prep.py.

    Parameters:
        locations -- a Locations object for the working directory
        args -- an argparse.Namespace object that has had the ParameterCache merged into it
    '''
    paramsCache = ParameterCache(locations.workingDirectory)
    paramsCache.load()

    # Present standard parameters
    print("# Parameters:")
    print(f"Working directory: {locations.workingDirectory}")

    if args.bamFiles != None and args.bamFiles != []:
        print(f"BAM files: {args.bamFiles}")
    else:
        print("BAM files: None")

    # Present metadata cache
    print() # blank line for spacing
    print("# Default metadata details:")
    metaSamples = None
    if args.metadataFile is not None:
        metadataCache = MetadataCache(locations.workingDirectory)
        metadataCache.establish()
        if metadataCache.metadataFile == None:
            print("## Metadata cache not found; re-initialising...")
            metadataCache.metadataFile = args.metadataFile

        print(f"Metadata file: {args.metadataFile}")
        print(f"Group 1 samples (n={len(metadataCache.group1)}): {metadataCache.group1}")
        print(f"Group 2 samples (n={len(metadataCache.group2)}): {metadataCache.group2}")
        metaSamples = set(metadataCache.group1 + metadataCache.group2)
    else:
        print("Metadata file: None (provide metadata to psQTL_proc.py with --meta)")

    # Present VCF cache
    print() # blank line for spacing
    print("# 'Call' VCF details:")
    vcfCache = None
    if args.vcfFile is not None:
        vcfCache = VcfCache(locations.workingDirectory)
        vcfCache.establish()
        if vcfCache.vcfFile == None:
            print("## VCF cache not found; re-initialising...")
            vcfCache.vcfFile = args.vcfFile

        print(f"VCF file: {args.vcfFile}")
        print(f"Num. variants: {vcfCache.variants}")
        print(f"Samples (n={len(vcfCache.samples)}): {vcfCache.samples}")
        print(f"Contigs (n={len(vcfCache.contigs)}): {vcfCache.contigs}")
    else:
        print("VCF file: None")

    print() # blank line for spacing
    print("# Filtered 'call' VCF details:")
    if args.filteredVcfFile is not None:
        vcfCache = VcfCache(locations.workingDirectory)
        vcfCache.establish()
        if vcfCache.filteredVcfFile == None:
            print("## VCF cache not found; re-initialising...")
            vcfCache.filteredVcfFile = args.filteredVcfFile

        print(f"Filtered VCF file: {args.filteredVcfFile}")
        if paramsCache.qualFilter == None:
            print(f"Quality filter: unknown")
        else:
            print(f"Quality filter: {paramsCache.qualFilter}")
        print(f"Num. filtered variants: {vcfCache.filteredVariants}")
        print(f"Filtered samples (n={len(vcfCache.filteredSamples)}): {vcfCache.filteredSamples}")
        print(f"Filtered contigs (n={len(vcfCache.filteredContigs)}): {vcfCache.filteredContigs}")
    else:
        print("Filtered VCF file: None")
    print() # blank line for spacing

    # Present depth cache
    print("# Depth VCF-like details:")
    depthCache = None
    if args.depthFile is not None:
        depthCache = DepthCache(locations.workingDirectory)
        depthCache.establish()
        if depthCache.depthFile == None:
            print("## Depth cache not found; re-initialising...")
            depthCache.depthFile = args.depthFile

        print(f"Depth file: {args.depthFile}")
        if paramsCache.windowSize == None:
            print(f"Window size: unknown")
        else:
            print(f"Window size: {paramsCache.windowSize} bp")
        print(f"Total num. bins: {depthCache.bins}")
        print(f"Num. bins with CNV: {depthCache.cnvBins}")
        print(f"Samples (n={len(depthCache.samples)}): {depthCache.samples}")
        print(f"Contigs (n={len(depthCache.contigs)}): {depthCache.contigs}")
    else:
        print("Depth file: None")
    print() # blank line for spacing

    # Identify potential conflicts or issues
    issues = []
    vcfSamples, vcfContigs, depthSamples, depthContigs = None, None, None, None
    if vcfCache != None:
        vcfSamples = set(vcfCache.filteredSamples) if vcfCache.filteredSamples != None else set(vcfCache.samples)
        vcfContigs = set(vcfCache.filteredContigs) if vcfCache.filteredContigs != None else set(vcfCache.contigs)
    if depthCache != None:
        depthSamples = set(depthCache.samples)
        depthContigs = set(depthCache.contigs)

    ## Sample issues
    if metaSamples != None and len(metaSamples) == 2:
        issues.append("Metadata only indicates two samples; analysis interpretation may be limited")
    if vcfSamples == None:
        issues.append("VCF file may need to be set or generated before 'call' analysis can proceed")
    elif metaSamples != None and metaSamples != vcfSamples:
        issues.append(f"Metadata samples (n={len(metaSamples)}) do not match VCF samples (n={len(vcfSamples)})")
    if depthSamples == None:
        issues.append("Depth file may need to be set or generated before 'depth' analysis can proceed")
    elif metaSamples != None and metaSamples != depthSamples:
        issues.append(f"Metadata samples (n={len(metaSamples)}) do not match depth samples (n={len(depthSamples)})")
    if vcfSamples != None and depthSamples != None and vcfSamples != depthSamples:
        issues.append(f"VCF samples (n={len(vcfSamples)}) do not match depth samples (n={len(depthSamples)})")

    ## Contig issues
    if vcfContigs != None and depthContigs != None and vcfContigs != depthContigs:
        issues.append(f"VCF contigs (n={len(vcfContigs)}) do not match depth contigs (n={len(depthContigs)})")

    # Present potential issues
    print("# Potential issues:")
    if issues == []:
        print("No issues found!")
    else:
        print("\n".join(issues))
    print() # blank line for spacing

def find_default_run(workingDirectory, defaultMetadataFile, runs):
    '''
    Mirrors runs.resolve_run() to identify which run psQTL_proc.py and psQTL_post.py would
    use if neither --run nor --meta were specified, without modifying anything.
    '''
    if defaultMetadataFile != None and os.path.isfile(defaultMetadataFile):
        try:
            return find_matching_run(workingDirectory, parse_metadata(defaultMetadataFile))
        except ValueError:
            return None
    return runs[0] if len(runs) == 1 else None

def summarise_results(workingDirectory, run):
    '''
    Returns:
        summary -- a one-line string listing the completed statistics and number of post
                   outputs within the run
    '''
    locations = Locations(workingDirectory, run.runDir)
    edTypes = [ label for label, edFile in [["alleles", locations.allelesEdFile],
                                            ["genotypes", locations.genotypesEdFile],
                                            ["inheritance", locations.inheritanceEdFile]]
                if is_complete(edFile) ]
    edResults = [ f"call ({', '.join(edTypes)})" ] if edTypes != [] else []
    if is_complete(locations.depthEdFile):
        edResults.append("depth")
    splsdaResults = [ label for label, selectedFile in [["call", locations.variantSplsdaSelectedFile],
                                                        ["depth", locations.depthSplsdaSelectedFile],
                                                        ["integrative", locations.integrativeSplsdaSelectedFile]]
                      if is_complete(selectedFile) ]
    numPost = len(os.listdir(locations.postDir)) if os.path.isdir(locations.postDir) else 0
    return (f"ED: {', '.join(edResults) if edResults != [] else 'none'} | " +
            f"sPLS-DA: {', '.join(splsdaResults) if splsdaResults != [] else 'none'} | " +
            f"post outputs: {numPost}")

def view_runs(workingDirectory, defaultMetadataFile, scriptName):
    '''
    Presents a summary of each analysis run created by psQTL_proc.py.
    '''
    print("# Analysis runs:")
    legacyFiles, legacySplsdaDir = find_legacy_results(workingDirectory)
    if legacyFiles != [] or legacySplsdaDir != None:
        print("Results from an earlier version of psQTL were found in the working directory; these will " +
              "be moved into an analysis run the next time psQTL_proc.py or psQTL_post.py is run")

    runs = list_runs(workingDirectory)
    if runs == []:
        print("None yet; psQTL_proc.py will create one when it is run")
    else:
        defaultRun = find_default_run(workingDirectory, defaultMetadataFile, runs)
        for run in runs:
            isDefault = defaultRun != None and run.number == defaultRun.number
            print(run.describe() + (" [default]" if isDefault else ""))
            print(f"    {summarise_results(workingDirectory, run)}")
        print(f"Use '{scriptName} -d {workingDirectory} --run <number>' to see the samples and outputs of a run")
        if len(runs) > 1:
            print("Use --run or --meta with psQTL_proc.py and psQTL_post.py to choose which run to work with" +
                  ("; the [default] run is used if neither is given" if defaultRun != None else ""))
    print() # blank line for spacing

def view_run_details(workingDirectory, run, inputSamples=None):
    '''
    Presents the metadata, segregation statistics, and post-processing outputs of a run.

    Parameters:
        workingDirectory -- a string indicating the psQTL working directory
        run -- a runs.Run object
        inputSamples -- (OPTIONAL) a set of sample IDs found in the VCF and/or depth files,
                        used to report which samples this run excludes
    '''
    locations = Locations(workingDirectory, run.runDir)
    group1, group2 = run.metadataDict["group1"], run.metadataDict["group2"]

    print(f"# Analysis run '{run.name}':")
    print(f"Directory: {run.runDir}")
    print(f"Metadata file: {run.metadataFile}")
    if run.metadataSource != None:
        sourceNote = {True: " (since modified)", False: "", None: " (no longer exists)"}[run.sourceChanged]
        print(f"Copied from: {run.metadataSource}{sourceNote}")
    if run.migrated:
        print("Origin: migrated from an earlier psQTL version")
    elif run.created != None:
        print(f"Created: {run.created}")
    print(f"Group 1 samples (n={len(group1)}): {group1}")
    print(f"Group 2 samples (n={len(group2)}): {group2}")
    if inputSamples != None:
        excluded = sorted(inputSamples.difference(group1 + group2))
        print(f"Input samples not in this run (n={len(excluded)}): {excluded}")
    print() # blank line for spacing

    print("# Segregation statistics:")
    for label, statFile in [["ED 'call' alleles", locations.allelesEdFile],
                            ["ED 'call' genotypes", locations.genotypesEdFile],
                            ["ED 'call' inheritance", locations.inheritanceEdFile],
                            ["ED 'depth'", locations.depthEdFile],
                            ["sPLS-DA 'call'", locations.variantSplsdaSelectedFile],
                            ["sPLS-DA 'depth'", locations.depthSplsdaSelectedFile],
                            ["sPLS-DA integrative", locations.integrativeSplsdaSelectedFile]]:
        print(f"{label}: {statFile if is_complete(statFile) else 'None'}")
    print() # blank line for spacing

    print(f"# Post-processing outputs ({locations.postDir}):")
    postFiles = sorted(os.listdir(locations.postDir)) if os.path.isdir(locations.postDir) else []
    print("\n".join(postFiles) if postFiles != [] else "None")
    print() # blank line for spacing

def get_input_samples(workingDirectory):
    '''
    Returns:
        inputSamples -- a set of sample IDs found in the cached VCF and/or depth files, or
                        None if neither has been cached
    '''
    inputSamples = None
    for cacheClass in [VcfCache, DepthCache]:
        cache = cacheClass(workingDirectory)
        try:
            cache.load()
        except FileNotFoundError:
            continue
        for samples in [ getattr(cache, attr, None) for attr in ["filteredSamples", "samples"] ]:
            if samples != None:
                inputSamples = set(samples) if inputSamples == None else inputSamples.union(samples)
                break
    return inputSamples
