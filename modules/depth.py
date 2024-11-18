import os
import pandas as pd
import numpy as np
from .parsing import parse_binned_tsv

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
                   where 0 indicates homozygous deletion, 1 indicates hemizygous
                   deletion, and 2 indicates homozygous presence
    '''
    HEMI_CUTOFF = 2.5
    HOMO_CUTOFF = 8
    
    # Get breakpoints for deletions and presence
    depths = np.array(list(binDict.values()))
    medianDepth = np.median(depths)
    heteroDepth = medianDepth / HEMI_CUTOFF
    homoDepth = medianDepth / HOMO_CUTOFF
    
    # Predict deletions and presence
    hetero = np.where((depths > homoDepth) & (depths <= heteroDepth), 1, 0)
    homoPresent = np.where(depths > heteroDepth, 2, 0)
    alleles = hetero + homoPresent
    
    # Return the results
    return alleles

def convert_alleles_to_gt(alleles):
    '''
    Converts the allele predictions to a genotype array, where 0 indicates
    homozygous deletion, 1 indicates heterozygous deletion, and 2 indicates
    homozygous presence.
    
    Parameters:
        alleles -- a numpy array where 0 indicates homozygous deletion,
                   1 indicates heterozygous deletion, and 2 indicates
                   homozygous presence
    Returns:
        genotypes -- a list of strings with the same length as the input array,
                     where each string is a VCF-encoded genotype where '0/0'
                     indicates homozygous presence, '0/1' indicates heterozygous
                     deletion, and '1/1' indicates homozygous deletion
    '''
    genotypes = [
        "1/1" if allele == 0 # homozygous deletion
        else "0/1" if allele == 1 # heterozygous deletion
        else "0/0" # homozygous presence
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
    open(outputFileName + ".ok", "w").close() # touch a .ok file to indicate success
