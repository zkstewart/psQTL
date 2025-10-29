import os, sys
import pandas as pd
import numpy as np

from .parsing import parse_binned_tsv, read_gz_file, WriteGzFile
from .ncls import WindowedNCLS

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _version import __version__

def get_median_value(values):
    '''
    Parameters:
        values -- a numpy array of depth values
    Returns:
        medianValue -- the median value of the input array, or 1 if all values are 0;
                       if the median would be 0, it will return the median when ignoring
                       zeros; if there are no non-zero values, this function will return 1.
    '''
    assert isinstance(values, np.ndarray), "Input values must be a numpy array."
    
    medianValue = np.median(values)
    if medianValue == 0:
        nonZeroValues = values[values != 0] # remove zeros
        # If all values are zero, make sure we return 1
        if len(nonZeroValues) == 0:
            medianValue = 1
        # Calculate median of non-zero values
        else:
            medianValue = np.median(nonZeroValues)
    return medianValue

def convert_depth_to_alleles(binDict, ploidy=2):
    '''
    Receives a histogram dictionary and normalises depth values with respect
    to the median. Then multiples the normalised depth values by the
    ploidy number to estimate the number of allele copies present.
    
    Parameters:
        binDict -- a dictionary with structure like:
                     {
                         position1: depth1,
                         position2: depth2,
                             ...
                    }
        ploidy -- (OPTIONAL) an integer indicating the ploidy number of the
                  samples being analysed e.g., '2' for diploid, '4' for tetraploid, etc.
    Returns:
        alleles -- a list with the same length as the input dictionary,
                   where the integer value indicates the number of allele copies
                   present given the indicated ploidy number.
    '''
    depths = np.array(list(binDict.values()))
    
    # Median-normalise the coverage
    medianDepth = get_median_value(depths)
    depths = depths / medianDepth
    
    # Round to nearest number of alleles
    alleles = np.round(depths * ploidy)
    
    # Return the results
    return [ int(x) for x in alleles ] # convert to a list of integers

def split_copynum_by_ploidy(n, k):
    '''
    See https://stackoverflow.com/questions/70392403/dividing-an-even-number-into-n-parts-each-part-being-a-multiple-of-2
    
    Parameters:
        n -- an integer indicating the number of gene copies present
        k -- an integer indicating the ploidy number of the samples
    '''
    d,r = divmod(n, k)
    return [d+1]*r + [d]*(k-r)

def convert_alleles_to_gt(alleles, ploidy=2):
    '''
    Converts allele copy number predictions to VCF-encoded genotypes.
    
    Parameters:
        alleles -- a numpy array where the integer value indicates the number of
                   allele copies present (assuming the ploidy value).
        ploidy -- (OPTIONAL) an integer indicating the ploidy number of the
                  samples being analysed e.g., '2' for diploid, '4' for tetraploid, etc.
    Returns:
        genotypes -- a list of strings with the same length as the input array,
                     where each string is a VCF-encoded genotype where the
                     summed value corresponds to the number of allele copies
                     present (assuming diploidy). For example, "0/0" indicates
                     0 copies, "0/1" indicates 1 copy, "1/1" indicates 2 copies,
                     and so on.
    '''
    genotypes = [
        "/".join(map(str, sorted(split_copynum_by_ploidy(allele, ploidy))))
        for allele in alleles
    ]
    
    return genotypes

def validate_df(df):
    '''
    Validates a DataFrame to ensure that all rows have consistent list lengths
    across all columns, with the exception of NaN values. This is a core assumption
    necessary for the subsequent processing of the DataFrame especially with the
    explode_dedup() function.
    '''
    for row in df.itertuples():
        rowListLength = None
        for rowColumn in row._fields[1:]:
            rowValue = getattr(row, rowColumn)
            if isinstance(rowValue, list):
                if rowListLength is None:
                    rowListLength = len(rowValue)
                elif len(rowValue) != rowListLength:
                    raise ValueError(f"Contig '{row}' has an inconsistent number of bins: {len(rowValue)} vs {rowListLength}; " + 
                                     "this indicates some sort of discrepancy in the input data, " +
                                     "possibly a different reference genome used for binned depth files " + 
                                     "or a different binning window size used for some of the samples. " +
                                     "Suggestion is to delete all depth files and re-run 'psQTL_prep.py depth'.")
            elif pd.isna(rowValue):
                continue
            else:
                raise ValueError("An unhandled error occurred when checking your binned depth files. " +
                                 "Suggestion is to delete all depth files and re-run 'psQTL_prep.py depth'. " +
                                 "If this error persists, please report it in the GitHub issues page.")

def explode_dedup(s):
    '''
    See https://stackoverflow.com/questions/77996844/how-to-explode-a-pandas-dataframe-that-has-nulls-in-some-rows-but-populated-in
    
    The assumptions of this function are checked by validate_df().
    '''
    s = s.explode()
    return s.set_axis(
        pd.MultiIndex.from_arrays([s.index, s.groupby(level=0).cumcount()])
    )

def call_cnvs_from_depth(samplePairs, outputFileName, windowSize, ploidy=2):
    '''
    Calls allele copy numbers (and hence CNVs) based on binned depth files and writes
    the results to a VCF-like file.
    
    Parameters:
        samplePairs -- a list of lists containing paired strings, where the first
                       string is the sample name and the second string is the path
                       to the binned depth file.
        outputFileName -- the path to the output VCF-like file.
        windowSize -- an integer specifying the size of the bins in each binned depth
                      files.
        ploidy -- (OPTIONAL) an integer indicating the ploidy number of the
                  samples being analysed e.g., '2' for diploid, '4' for tetraploid, etc.
    '''
    genotypesDict = {}
    samples = []
    for sampleName, binFile in samplePairs:
        # Error out if bin file is missing
        if not os.path.isfile(binFile):
            raise FileNotFoundError(f"Binned depth file '{binFile}' not found!")
        
        # Parse the binned depth file
        histoDict = parse_binned_tsv(binFile)
        
        # Predict allele copy numbers and encode as genotypes
        genotypesDict[sampleName] = {}
        for contigID, binDict in histoDict.items():
            alleles = convert_depth_to_alleles(binDict, ploidy=ploidy)
            genotypes = convert_alleles_to_gt(alleles, ploidy=ploidy)
            genotypesDict[sampleName][contigID] = genotypes
        samples.append(sampleName)
    
    # Convert to DataFrame and validate file structure
    df = pd.DataFrame(genotypesDict)
    validate_df(df)
    
    # Explode the DataFrame to have one row per window/bin
    exploded_df = pd.concat({c: explode_dedup(df[c]) for c in df}, axis=1).reset_index()
    exploded_df.drop(columns="level_1", inplace=True)
    exploded_df.fillna(value=".", inplace=True)
    
    # Convert to VCF-like format
    exploded_df.rename(columns={"level_0": "#CHROM"}, inplace=True)
    exploded_df["POS"] = (exploded_df.groupby("#CHROM").cumcount() * windowSize) + 1 # 1-based indexing for VCF
    exploded_df["ID"] = "."
    exploded_df["REF"] = "N"
    exploded_df["ALT"] = "N"
    exploded_df["QUAL"] = "."
    exploded_df["FILTER"] = "."
    exploded_df["INFO"] = "."
    exploded_df["FORMAT"] = "GT"
    exploded_df = exploded_df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]]
    
    # Write to output file
    todaysDate = pd.Timestamp.now().strftime("%d-%m-%Y")
    with WriteGzFile(outputFileName) as fileOut:
        fileOut.write("##fileformat=VCF-like\n")
        fileOut.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Sums to number of median-normalised allele copies\">\n")
        fileOut.write(f"##psQTL_prep;module=depth, version={__version__}, windowSize={windowSize}; ploidy={ploidy}; DD-MM-YYYY={todaysDate}\n")
        exploded_df.to_csv(fileOut, sep="\t", index=False)

def parse_bins_as_dict(depthFileDict):
    '''
    Parameters:
        depthFileDict -- a dictionary with structure like:
                         {
                             "group1": [
                                 ["sample1", "depthFile1"],
                                 ["sample2", "depthFile2"],
                                 ...
                            ],
                             "group2": [ ... ]
                         }
    Returns:
        coverageDict -- a dictionary with structure like:
                  {
                      "group1": {
                          "sample1": {
                              "chr1": [[pos1, pos2, ...], [coverage1, coverage2, ...]],
                              "chr2": [[pos1, pos2, ...], [coverage1, coverage2, ...]],
                              ...,
                            },
                            "sample2": { ... },
                            ...
                      },
                      "group2": { ... }
                  }
    '''
    coverageDict = {}
    for group, depthFiles in depthFileDict.items():
        coverageDict[group] = {}
        # Iterate through group files
        for sampleID, depthFile in depthFiles:
            coverageDict[group][sampleID] = {}
            # Parse binned depth file
            with read_gz_file(depthFile) as fileIn:
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
                    if not contigID in coverageDict[group][sampleID]:
                        coverageDict[group][sampleID][contigID] = [[], []]
                    coverageDict[group][sampleID][contigID][0].append(pos)
                    coverageDict[group][sampleID][contigID][1].append(coverage)
    return coverageDict

def normalise_coverage_dict(coverageDict):
    '''
    Modifies the coverage dictionary to contain median-normalised coverage values.
    
    Parameters:
        coverageDict -- a dictionary with structure like:
                  {
                      "group1": {
                          "sample1": {
                              "chr1": [[pos1, pos2, ...], [coverage1, coverage2, ...]],
                              "chr2": [[pos1, pos2, ...], [coverage1, coverage2, ...]],
                              ...,
                            },
                            "sample2": { ... },
                            ...
                      },
                      "group2": { ... }
                  }
    '''
    for group, sampleDict in coverageDict.items():
        for sampleID, depthDict in sampleDict.items():
            for chrom, value in depthDict.items():
                coverages = np.array(value[1])
                
                # Figure out what our median value is
                medianCoverage = get_median_value(coverages)
                
                # Normalise the coverage values and store them
                coverages = coverages / medianCoverage
                coverageDict[group][sampleID][chrom][1] = coverages

def convert_dict_to_depthncls(coverageDict, windowSize):
    '''
    Parameters:
        coverageDict -- a dictionary with structure like:
                  {
                      "group1": {
                          "sample1": {
                              "chr1": [[pos1, pos2, ...], [coverage1, coverage2, ...]],
                              "chr2": [[pos1, pos2, ...], [coverage1, coverage2, ...]],
                              ...,
                            },
                            "sample2": { ... },
                            ...
                      },
                      "group2": { ... }
                  }
        windowSize -- an integer indicating the size of the windows used for binning
    Returns:
        depthNCLSDict -- a dictionary with structure like:
                         {
                             "group1": {
                                 "sample1": WindowedNCLS,
                                 "sample2": WindowedNCLS,
                                 ...
                            },
                             "group2": { ... }
                         }
    '''
    depthNCLSDict = {}
    for group, sampleDict in coverageDict.items():
        depthNCLSDict[group] = {}
        for sampleID, depthDict in sampleDict.items():
            windowedNCLS = WindowedNCLS(windowSize)
            for chrom, value in depthDict.items():
                positions = np.array(value[0])
                edValues = np.array(value[1])
                windowedNCLS.add(chrom, positions, edValues)
                depthNCLSDict[group][sampleID] = windowedNCLS
    return depthNCLSDict
