import os
import pandas as pd
import numpy as np
from .parsing import parse_binned_tsv
from .ed import EDNCLS

def get_median_value(values):
    medianValue = np.median(values)
    if medianValue == 0:
        nonZeroValues = values[values != 0]
        if len(nonZeroValues) == 0:
            medianValue = 1
        else:
            medianValue = np.median(nonZeroValues)
    return medianValue

def predict_deletions(binDict):
    '''
    Receives a histogram dictionary and predicts regions of homozygous deletion,
    homozygous presence, and hemizygous regions on the basis of depth coverage.
    Uses a simple heuristic approach.
    
    Parameters:
        binDict -- a dictionary with structure like:
                     {
                         position1: depth1,
                         position2: depth2,
                             ...
                    }
    Returns:
        alleles -- a numpy array with the same length as the input dictionary,
                   where the integer value indicates the number of allele copies
                   present (assuming diploidy). For example, 0 indicates
                   homozygous deletion, 1 indicates hemizygous deletion, and 2
                   indicates homozygous presence. Increasing values correspond to
                   increasing copy number.
    '''
    depths = np.array(list(binDict.values()))
    
    # Median-normalise the coverage
    medianDepth = get_median_value(depths)
    depths = depths / medianDepth
    
    # Round to nearest number of alleles
    alleles = np.round(depths * 2)
    
    # Return the results
    return alleles

def convert_alleles_to_gt(alleles):
    '''
    Converts the allele predictions to a genotype array, where 0 indicates
    homozygous deletion, 1 indicates heterozygous deletion, and 2 indicates
    homozygous presence.
    
    Parameters:
        alleles -- a numpy array where the integer value indicates the number of
                   allele copies present (assuming diploidy). For example, 0
                   indicates homozygous deletion, 1 indicates hemizygous deletion,
                   2 indicates homozygous presence, and increasing values correspond
                   to increasing copy number.
    Returns:
        genotypes -- a list of strings with the same length as the input array,
                     where each string is a VCF-encoded genotype where the
                     summed value corresponds to the number of allele copies
                     present (assuming diploidy). For example, "0/0" indicates
                     0 copies, "0/1" indicates 1 copy, "1/1" indicates 2 copies,
                     and so on.
    '''
    genotypes = [
        f"{int(np.floor(allele/2))}/{int(np.ceil(allele/2))}"
        for allele in alleles
    ]
    return genotypes

def call_deletions_from_depth(samplePairs, outputFileName, windowSize):
    '''
    Calls homozygous and hemizygous deletions from binned depth files and writes
    the results to a VCF-like file.
    
    Parameters:
        samplePairs -- a list of lists containing paired strings, where the first
                       string is the sample name and the second string is the path
                       to the binned depth file.
        outputFileName -- the path to the output VCF-like file.
        windowSize -- an integer specifying the size of the bins in each binned depth
                      files.
    '''
    genotypesDict = {}
    samples = []
    for sampleName, binFile in samplePairs:
        # Error out if bin file is missing
        if not os.path.isfile(binFile):
            raise FileNotFoundError(f"Binned depth file '{binFile}' not found!")
        
        # Parse the binned depth file
        histoDict = parse_binned_tsv(binFile)
        
        # Predict deletions
        genotypesDict[sampleName] = {}
        for contigID, binDict in histoDict.items():
            alleles = predict_deletions(binDict)
            genotypes = convert_alleles_to_gt(alleles)
            genotypesDict[sampleName][contigID] = genotypes
        samples.append(sampleName)
    
    # Convert to DataFrame with VCF-like format
    df = pd.DataFrame(genotypesDict)
    exploded_df = df.apply(lambda x: x.explode()).reset_index()
    exploded_df.rename(columns={"index": "#CHROM"}, inplace=True)    
    exploded_df["POS"] = exploded_df.groupby("#CHROM").cumcount() * windowSize
    exploded_df["ID"] = "."
    exploded_df["REF"] = "N"
    exploded_df["ALT"] = "N"
    exploded_df["QUAL"] = "."
    exploded_df["FILTER"] = "."
    exploded_df["INFO"] = "."
    exploded_df["FORMAT"] = "GT"
    exploded_df = exploded_df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]]
    
    # Write to output file
    with open(outputFileName, "w") as fileOut:
        fileOut.write("##fileformat=VCF-like\n")
        fileOut.write("##0/0=nodeletion\n")
        fileOut.write("##0/1=hemizygousdeletion\n")
        fileOut.write("##1/1=homozygousdeletion\n")
        fileOut.write("##psQTL_prepDeletionPrediction\n")
        exploded_df.to_csv(fileOut, sep="\t", index=False)

def parse_bins_as_dict(depthFileDict, windowSize):
    '''
    Parameters:
        depthFileDict -- a dictionary with structure like:
                         {
                             "bulk1": [
                                 ["sample1", "depthFile1"],
                                 ["sample2", "depthFile2"],
                                 ...
                            ],
                             "bulk2": [ ... ]
                         }
        windowSize -- an integer indicating the size of the windows used for binning
    Returns:
        coverageDict -- a dictionary with structure like:
                  {
                      "bulk1": {
                          "sample1": {
                              "chr1": [[pos1, pos2, ...], [coverage1, coverage2, ...]],
                              "chr2": [[pos1, pos2, ...], [coverage1, coverage2, ...]],
                              ...,
                            },
                            "sample2": { ... },
                            ...
                      },
                      "bulk2": { ... }
                  }
    '''
    coverageDict = {}
    for bulk, depthFiles in depthFileDict.items():
        coverageDict[bulk] = {}
        # Iterate through bulk files
        for sampleID, depthFile in depthFiles:
            coverageDict[bulk][sampleID] = {}
            # Parse binned depth file
            with open(depthFile, "r") as fileIn:
                for line in fileIn:
                    # Extract relevant details
                    contigID, pos, coverage = line.strip().split("\t")
                    try:
                        pos = int(pos)
                    except:
                        raise ValueError(f"Position '{pos}' is not an integer in file '{depthFile}'")
                    try:
                        coverage = int(coverage)
                    except:
                        raise ValueError(f"Coverage '{coverage}' is not an integer in file '{depthFile}'")
                    
                    # Store the coverage
                    if not contigID in coverageDict[bulk][sampleID]:
                        coverageDict[bulk][sampleID][contigID] = [[], []]
                    coverageDict[bulk][sampleID][contigID][0].append(pos)
                    coverageDict[bulk][sampleID][contigID][1].append(coverage)
    return coverageDict

def normalise_coverage_dict(coverageDict):
    '''
    Modifies the coverage dictionary to contain median-normalised coverage values.
    
    Parameters:
        coverageDict -- a dictionary with structure like:
                  {
                      "bulk1": {
                          "sample1": {
                              "chr1": [[pos1, pos2, ...], [coverage1, coverage2, ...]],
                              "chr2": [[pos1, pos2, ...], [coverage1, coverage2, ...]],
                              ...,
                            },
                            "sample2": { ... },
                            ...
                      },
                      "bulk2": { ... }
                  }
    '''
    for bulk, sampleDict in coverageDict.items():
        for sampleID, depthDict in sampleDict.items():
            for chrom, value in depthDict.items():
                positions = np.array(value[0])
                coverages = np.array(value[1])
                
                # Figure out what our median value is
                medianCoverage = get_median_value(coverages)
                
                # Normalise the coverage values and store them
                coverages = coverages / medianCoverage
                coverageDict[bulk][sampleID][chrom][1] = coverages

def convert_dict_to_depthncls(coverageDict, windowSize):
    '''
    Parameters:
        coverageDict -- a dictionary with structure like:
                  {
                      "bulk1": {
                          "sample1": {
                              "chr1": [[pos1, pos2, ...], [coverage1, coverage2, ...]],
                              "chr2": [[pos1, pos2, ...], [coverage1, coverage2, ...]],
                              ...,
                            },
                            "sample2": { ... },
                            ...
                      },
                      "bulk2": { ... }
                  }
        windowSize -- an integer indicating the size of the windows used for binning
    Returns:
        depthNCLSDict -- a dictionary with structure like:
                         {
                             "bulk1": {
                                 "sample1": EDNCLS,
                                 "sample2": EDNCLS,
                                 ...
                            },
                             "bulk2": { ... }
                         }
    '''
    depthNCLSDict = {}
    for bulk, sampleDict in coverageDict.items():
        depthNCLSDict[bulk] = {}
        for sampleID, depthDict in sampleDict.items():
            edNCLS = EDNCLS(windowSize)
            for chrom, value in depthDict.items():
                positions = np.array(value[0])
                edValues = np.array(value[1])
                edNCLS.add(chrom, positions, edValues)
                depthNCLSDict[bulk][sampleID] = edNCLS
    return depthNCLSDict
