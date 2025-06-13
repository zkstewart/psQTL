import numpy as np
from math import sqrt, ceil
from collections import Counter
from itertools import combinations

from .parsing import read_gz_file, vcf_header_to_metadata_validation, parse_vcf_genotypes
from .ncls import WindowedNCLS

def gt_median_adjustment(genotypeLists):
    '''
    Receives a list of genotype lists, obtains the median of all alleles found in the lists,
    and adjusts the genotypes to be relative to the median. This is specifically intended for
    use on CNV genotypes, which are often highly variable and range from 0 up to very high values.
    
    Parameters:
        genotypeLists -- a list of one or more lists containing genotype values as integers
                         with format like:
                         [
                             [
                                [0, 1],
                                [0, 0],
                                [1, 1],
                                ...
                             ],
                             [ ... ],
                             ...
                         ]
    Returns:
        adjustedGenotypes -- a list of lists with the same structure as genotypeLists
                             but with genotypes adjusted to be relative to the median of all
                             alleles found in the genotypeLists
    '''
    # Flatten the genotype lists and calculate the median allele
    alleles = [ allele for sublist in genotypeLists for gt in sublist if not "." in gt for allele in gt ]
    if len(alleles) == 0:
        medianAllele = 0 # if no alleles are found, set median to 0
    else:
        medianAllele = np.median(alleles)
    
    # Adjust the genotypes relative to the median allele
    adjustedGenotypes = []
    for sublist in genotypeLists:
        adjustedSublist = []
        for gt in sublist:
            gtList = []
            if "." in gt:
                gtList.append(".")
            else:
                for allele in gt:
                    adjustedAllele = 0 if allele <= medianAllele else 1
                    gtList.append(adjustedAllele)
            adjustedSublist.append(gtList)
        adjustedGenotypes.append(adjustedSublist)
    return adjustedGenotypes

def calculate_segregant_ed(b1Gt, b2Gt, isCNV=False, parentsGT=None):
    '''
    Parameters:
        b1Gt / b2Gt -- a list of lists containing the genotype value as integers
                       with format like:
                       [
                           [0, 1],
                           [0, 0],
                           [1, 1],
                           ...
                       ]
        isCNV -- (OPTIONAL) a boolean indicating whether the genotypes are for CNVs
                 (True) or SNPs/indels (False); default is False
        parentsGT -- a list of two lists containing the genotype value as integers
                     for the parents with format like:
                     [ [0, 1], [1, 2] ]
                     OR None if no parents are available or specified to
                     filter out non-inheritable genotypes
    Returns:
        numAllelesB1 -- the number of genotyped alleles in bulk 1
        numAllelesB2 -- the number of genotyped alleles in bulk 2
        edist -- a float of the the Euclidean distance between the two bulks
    '''
    # Adjust values if this is a CNV
    if isCNV:
        b1Gt, b2Gt = gt_median_adjustment([b1Gt, b2Gt])
    
    # Filter impossible progeny genotypes based on the parents' genotypes
    if parentsGT != None and len(parentsGT) == 2:
        possibleGTs = possible_genotypes(parentsGT[0], parentsGT[1])
        b1Gt = [ gt for gt in b1Gt if set(gt) in possibleGTs ]
        b2Gt = [ gt for gt in b2Gt if set(gt) in possibleGTs ]
    
    # Calculate Euclidean distance between the two bulks with both methods
    numAllelesB1, numAllelesB2, edist1 = calculate_allele_frequency_ed(b1Gt, b2Gt)
    numSamplesB1, numSamplesB2, edist2 = calculate_genotype_frequency_ed(b1Gt, b2Gt)
    
    # Return the values
    "Choose the method of comparison that yields the largest distance"
    return numAllelesB1, numAllelesB2, max(edist1, edist2)

def calculate_allele_frequency_ed(b1Gt, b2Gt):
    '''
    Separate function to calculate the Euclidean distance between two bulks
    based on the allele frequencies of the genotypes in each bulk.
    
    Parameters:
        b1Gt / b2Gt -- a list of lists containing the genotype value as integers
                       with format like:
                       [
                           [0, 1],
                           [0, 0],
                           [1, 1],
                           ...
                       ]
    Returns:
        numAllelesB1 -- the number of genotyped alleles in bulk 1
        numAllelesB2 -- the number of genotyped alleles in bulk 2
        edist -- a float of the the Euclidean distance between the two bulks
    '''
    # Get all the unique alleles
    alleles = [ allele for gt in b1Gt + b2Gt for allele in gt ]
    uniqueAlleles = list(set(alleles))
    
    # Tally for bulk 1
    b1Count = { allele: 0 for allele in uniqueAlleles }
    for gt in b1Gt:
        for allele in gt:
            b1Count[allele] += 1
    
    # Tally for bulk 2
    b2Count = { allele: 0 for allele in uniqueAlleles }
    for gt in b2Gt:
        for allele in gt:
            b2Count[allele] += 1
    
    # Sum the number of genotyped alleles for each bulk
    numAllelesB1 = sum(b1Count.values())
    numAllelesB2 = sum(b2Count.values())
    
    # Calculate the Euclidean distance between the two bulks if possible
    if numAllelesB1 == 0 or numAllelesB2 == 0:
        return numAllelesB1, numAllelesB2, 0 # euclidean distance cannot be calculated
    else:
        # Derive our euclidean distance value
        "Refer to 'Euclidean distance calculation' in Hill et al. 2013"
        edist = sqrt(sum([
            ((b1Count[allele] / numAllelesB1) - (b2Count[allele] / numAllelesB2))**2
            for allele in uniqueAlleles
        ]))
        # Return the values
        return numAllelesB1, numAllelesB2, edist

def calculate_genotype_frequency_ed(b1Gt, b2Gt):
    '''
    Parameters:
        b1Gt / b2Gt -- a list of lists containing the genotype value as integers
                       with format like:
                       [
                           [0, 1],
                           [0, 0],
                           [1, 1],
                           ...
                       ]
    Returns:
        numSamplesB1 -- the number of genotyped samples in bulk 1
        numSamplesB2 -- the number of genotyped samples in bulk 2
        edist -- a float of the the Euclidean distance between the two bulks
    '''
    # Get all the unique genotypes
    uniqueGts = set(tuple(sorted(gt)) for gt in b1Gt + b2Gt)
    
    # Tally for bulk 1
    b1GtCount = { gt: 0 for gt in uniqueGts }
    for gt in b1Gt:
        b1GtCount[tuple(sorted(gt))] += 1
    
    # Tally for bulk 2
    b2GtCount = { gt: 0 for gt in uniqueGts }
    for gt in b2Gt:
        b2GtCount[tuple(sorted(gt))] += 1
    
    # Sum the number of genotyped samples for each bulk
    numSamplesB1 = len(b1Gt)
    numSamplesB2 = len(b2Gt)
    
    # Calculate the Euclidean distance between the two bulks if possible
    if numSamplesB1 == 0 or numSamplesB2 == 0:
        return numSamplesB1, numSamplesB2, 0 # euclidean distance cannot be calculated
    else:
        # Derive our euclidean distance value
        "This is the same as in calculate_allele_frequency_ed(), but for genotypes instead of alleles"
        edist = sqrt(sum([
            ((b1GtCount[gt] / numSamplesB1) - (b2GtCount[gt] / numSamplesB2))**2
            for gt in uniqueGts
        ]))
        
        # Return the values
        return numSamplesB1, numSamplesB2, edist

def possible_genotypes(gt1, gt2):
    '''
    Calculate all possible inheritable genotypes from two parent genotypes.
    
    Parameters:
        gt1 / gt2 -- a list of integers representing the genotype values
                     with format like:
                     [0, 1, 2, ...]
    Returns:
        possibleGTs -- a list of sets containing all possible genotypes
                       that can be formed from the two genotypes; sets
                       mask duplicated alleles but allow for testing
                       of genotype without regard to order
    '''
    numGt1 = int(len(gt1) / 2)
    numGt2 = int(len(gt2) / 2)
    
    possibleGTs = []
    for g1 in combinations(gt1, numGt1):
        for g2 in combinations(gt2, numGt2):
            possibleGTs.append(set(g1 + g2))
    
    return possibleGTs

def calculate_inheritance_ed(b1Gt, b2Gt, parentsGT):
    '''
    DEPRECATED; use calculate_segregant_ed() instead.
    
    Parameters:
        b1Gt / b2Gt -- a list of lists containing the genotype value as integers
                       with format like:
                       [
                           [0, 1],
                           [0, 2],
                           [1, 1],
                           ...
                       ]
        parentsGT -- a list of two lists containing the genotype value as integers
                     for the parents with format like:
                     [ [0, 1], [1, 2] ]
    Returns:
        numAllelesB1 -- the number of genotyped alleles in bulk 1
        numAllelesB2 -- the number of genotyped alleles in bulk 2
        edist -- a float of the the Euclidean distance between the two bulks
    '''
    # Raise error for unmanageable ploidy values
    for i, gt in enumerate(parentsGT):
        if len(gt) % 2 != 0:
            raise ValueError(f"Odd number of alleles found in parent {i+1}'s genotype; " +
                             "cannot calculate inheritance Euclidean distance with odd ploidy")
    PARENT_FULLY_ASSIGNED = { f"p{i+1}" : int(len(gt) / 2) for i, gt in enumerate(parentsGT) }
    
    if (len(b1Gt) > 0 and len(b1Gt[0]) != sum(PARENT_FULLY_ASSIGNED.values())) or \
       (len(b2Gt) > 0 and len(b2Gt[0]) != sum(PARENT_FULLY_ASSIGNED.values())):
        raise ValueError("Progeny samples must have a number of alleles equal to the summed value of " +
                         f"half of each parent's chromosomes ({sum(PARENT_FULLY_ASSIGNED.values())}); " +
                         "cannot calculate inheritance Euclidean distance with mismatched ploidy")
    
    # Determine the possible progeny genotypes based on the parents' genotypes
    possibleGTs = possible_genotypes(parentsGT[0], parentsGT[1])
    
    # Assign alleles to parent haplotypes
    bulkSums = []
    numSamples = [0, 0] # [numSamplesB1, numSamplesB2]
    for i, bulkGt in enumerate([b1Gt, b2Gt]): # for each bulk
        # Set up data structure to hold assigned allele counts across the bulk
        bulkSum = {
            "p1": {
                gt: 0 for gt in parentsGT[0]
            },
            "p2": {
                gt: 0 for gt in parentsGT[1]
            }
        }
        # For each sample, assign alleles to parents based on the most likely inheritance
        for gt in bulkGt:
            # Detect samples that cannot be matched to parents
            if not set(gt) in possibleGTs:
                continue
            numSamples[i] += 1 # increment the sample count for this bulk since it is valid
            
            # Set up the data structure for holding this sample's assigned alleles
            sampleColumns = {
                "p1": {
                    gt: 0 for gt in parentsGT[0]
                },
                "p2": {
                    gt: 0 for gt in parentsGT[1]
                }
            }
            # Order assignable parent allele by count; ensures that least common alleles are assigned first
            Palleles = Counter([ allele for parent in parentsGT for allele in parent ])
            Palleles = Palleles.most_common()
            Palleles.sort(key=lambda x: x[1]) # sort by count for rare alleles first
            Palleles = { allele: count for allele, count in Palleles }
            
            # Order sample alleles by the parent allele counts
            Salleles = sorted(gt, key=lambda x: Palleles[x])
            
            # Assign sample alleles to most likely parent alleles
            isValidSample = True
            parentAssigned = { parent: 0 for parent in sampleColumns } # how many alleles have been assigned to each parent
            for Sallele in Salleles: # for each sample allele
                # Detect progeny mismatch with parent genotype
                count = Palleles[Sallele]
                
                # Derive the likelihood of this allele being assigned to a specific parental haplotype
                "Fractional values indicate partial assignment of the allele to a haplotype"
                likelihood = 1 / count
                
                # Partially assign the allele to all applicable haplotypes
                for parent, assignedDict in sampleColumns.items():
                    # Skip assigned parents
                    if sum([ assignedDict[_allele] for _allele in assignedDict ]) == PARENT_FULLY_ASSIGNED[parent]:
                        continue
                    
                    if Sallele in assignedDict:
                        assignedDict[Sallele] += likelihood
                        
                        # Update the parent allele counts to reflect that one allele has been confidently assigned
                        if assignedDict[Sallele] == PARENT_FULLY_ASSIGNED[parent]:
                            parentGT = list(sampleColumns[parent].keys())
                            for pGT in parentGT:
                                Palleles[pGT] -= 1
            
            # Update the column sums using this sample
            for parent, assignedDict in sampleColumns.items():
                for allele, assigned in assignedDict.items():
                    bulkSum[parent][allele] += assigned
        
        # Set aside results for this bulk
        bulkSums.append(bulkSum)
    
    # Format bulkSums to ensure all alleles are present [makes it easier to calculate ED without nested if statements]
    uniqueAlleles = list(set([ allele for gt in parentsGT for allele in gt ]))
    for allele in uniqueAlleles:
        for bSum in bulkSums:
            if allele not in bSum["p1"]:
                bSum["p1"][allele] = 0
            if allele not in bSum["p2"]:
                bSum["p2"][allele] = 0
    b1Sum, b2Sum = bulkSums
    
    # Get the number of samples and alleles in each bulk
    numSamplesB1, numSamplesB2 = numSamples # we might have skipped samples that don't match parent genotypes
    numAllelesB1 = numSamplesB1 * len(b1Gt[0]) if numSamplesB1 > 0 else 0 # num samples * ploidy of the samples
    numAllelesB2 = numSamplesB2 * len(b2Gt[0]) if numSamplesB2 > 0 else 0 # conforms to the calculate_segregant_ed() return values
    
    # Calculate the Euclidean distance between the parental inherited alleles if possible
    """Note that values are divied by the number of samples, not the number of alleles;
    doing so gives results in the same scale as the standard segregant ED calculation,
    but the rationale is difficult for me to articulate"""
    if numSamplesB1 == 0 or numSamplesB2 == 0:
        return numAllelesB1, numAllelesB2, 0 # euclidean distance cannot be calculated
    else:
        edist = sqrt(sum([
            ((b1Sum[parent][allele] / numSamplesB1) - (b2Sum[parent][allele] / numSamplesB2))**2
            for parent in ["p1", "p2"]
            for allele in uniqueAlleles
        ])/2) # divide by 2 since there are two parents and 2x the number of comparions made across bulks
    
    # Return the values
    return numAllelesB1, numAllelesB2, edist

def parse_vcf_for_ed(vcfFile, metadataDict, isCNV, parents=[], ignoreIdentical=True, quiet=False):
    '''
    Parameters:
        vcfFile -- a string pointing to the VCF or VCF-like file to parse
        metadataDict -- a dictionary with structure like:
                        {
                            "bulk1": set([ "sample1", "sample2", ... ]),
                            "bulk2": set([ "sample3", "sample4", ... ])
                        }
        isCNV -- a boolean indicating whether the VCF file is for CNVs ("depth"; True)
                 or SNPs/indels ("call"; False)
        parents -- (OPTIONAL) a list of two sample IDs to use as parents for calculating
                   haplotype inheritance ED OR an empty list if no parents are available
                   or specified; default is an empty list for standard segregant ED calculation
        ignoreIdentical -- (OPTIONAL) a boolean indicating whether to ignore
                           identical non-reference alleles shared by all samples;
                           default is True, which means that identical non-reference
                           alleles will be ignored
        quiet -- (OPTIONAL) a boolean indicating whether to suppress output messages;
                 default is False, which means that output messages will be printed
    Yields:
        contig -- the contig name for the variant
        pos -- the position of the variant
        variant -- the type of variant (snp or indel)
        numAllelesB1 -- the number of genotyped alleles in bulk 1
        numAllelesB2 -- the number of genotyped alleles in bulk 2
        euclideanDist -- the Euclidean distance between the two bulks
    '''
    # Validations
    if parents == None: # just in case
        parents = []
    
    if parents != [] and len(parents) != 2:
        raise ValueError("Parents must be a list with zero (standard ED) or two sample IDs (haplotype inheritance ED)")
        for parent in parents:
            if parent not in metadataDict["bulk1"] and parent not in metadataDict["bulk2"]:
                raise ValueError(f"Parent sample '{parent}' is not in either bulk; cannot calculate haplotype inheritance ED")
    
    # Iterate through the VCF file
    samples = None
    with read_gz_file(vcfFile) as fileIn:
        for line in fileIn:
            l = line.strip('\r\n\t "') # remove quotations to help with files opened by Excel
            sl = l.replace('"', '').split("\t") # remove any remaining quotations
            
            # Skip blank lines
            if l == "":
                continue
            
            # Handle header line
            if l.startswith("#CHROM"):
                samples = sl[9:] # This gives us the ordered sample IDs
                vcf_header_to_metadata_validation(samples, metadataDict, strict=False, quiet=quiet)
            if l.startswith("#"):
                continue
            
            # Validate line length
            if len(sl) < 11: # 9 fixed columns + minimum 2 genotype columns
                raise ValueError(f"VCF file has too few columns; offending line is '{l}'")
            
            # Validate that we've seen the header line
            if samples == None:
                raise ValueError("VCF file does not contain a #CHROM header line; cannot parse!")
            
            # Extract relevant details from line
            contig = sl[0]
            pos = int(sl[1])
            ref = sl[3]
            alt = sl[4].split(",")
            ref_alt = [ref, *alt]
            try:
                qual = float(sl[5])
            except:
                qual = 0.0 # If the quality is missing, we'll just assume it's zero
            formatField = sl[8]
            sampleFields = sl[9:]
            
            # Parse the genotypes out of the sample fields
            snpDict = parse_vcf_genotypes(formatField, sampleFields, samples)
            
            # Figure out what type of variant this is
            if any([ x == "." for x in ref_alt ]):
                variant = "indel"
            elif all([ x == "N" for x in ref_alt ]): # for parsing deletion VCF-like files
                variant = "cnv"
            elif any([ len(ref_alt[0]) != len(ref_alt[x]) for x in range(1, len(ref_alt))]):
                variant = "indel"
            else:
                variant = "snp"
            
            # Split sample genotypes into bulk1 and bulk2
            bulk1 = [ snpDict[sample] for sample in metadataDict["bulk1"] if sample in snpDict and not sample in parents ]
            bulk2 = [ snpDict[sample] for sample in metadataDict["bulk2"] if sample in snpDict and not sample in parents ]
            parentsGT = [ snpDict[parent] for parent in parents if parent in snpDict ] # if parents == [] this will always be empty
            
            # Calculate Euclidean distance
            numAllelesB1, numAllelesB2, euclideanDist = calculate_segregant_ed(bulk1, bulk2,
                                                                               isCNV=isCNV, parentsGT=parentsGT)
            
            # Skip if both bulks are identical
            if ignoreIdentical and euclideanDist == 0:
                bulk1Dedup = set(( tuple(x) for x in bulk1 ))
                bulk2Dedup = set(( tuple(x) for x in bulk2 ))
                "if both have set len==1, are the same, and have 0 Euclidean distance, they are identical non-reference alleles"
                if len(bulk1Dedup) == 1 and bulk1Dedup == bulk2Dedup:
                    continue
            
            # Yield results
            yield contig, pos, variant, numAllelesB1, numAllelesB2, euclideanDist

def parse_ed_as_dict(edFile, metadataDict, missingFilter=0.5):
    '''
    Parameters:
        edFile -- a string indicating the path to an ED file
        metadataDict -- a dictionary with structure like:
                        {
                            "bulk1": set([ "sample1", "sample2", ... ]),
                            "bulk2": set([ "sample3", "sample4", ... ])
                        } OR None if no filtering is desired
        missingFilter -- OPTIONAL; a float indicating the maximum allowed missing data
                         calculated for each bulk
    Returns:
        edDict -- a dictionary with structure like:
                  {
                      "chr1": [[pos1, pos2, ...], [ed1, ed2, ...]],
                      "chr2": [[pos1, pos2, ...], [ed1, ed2, ...]],
                      ...
                  }
    '''
    HEADER = ["CHROM", "POSI", "variant", "bulk1_alleles", "bulk2_alleles", "euclideanDist"]
    
    # Make sure the metadata is valid if missingFilter is > 0
    if missingFilter > 0:
        if metadataDict == None:
            raise ValueError("Cannot filter for missing data without metadata")
        if not all([ x in metadataDict for x in ["bulk1", "bulk2"] ]):
            raise ValueError("Metadata dictionary must contain keys 'bulk1' and 'bulk2'")
    
    # Calculate how many alleles in each bulk
    if metadataDict != None:
        BULK1_ALLELES = len(metadataDict["bulk1"]) * 2
        BULK2_ALLELES = len(metadataDict["bulk2"]) * 2
        
        # Alert user to number of samples needed to pass filtration in each bulk
        print(f"# Filtering for missing data: up to {missingFilter*100}% missing data allowed in each bulk")
        print(f"# For bulk 1: {BULK1_ALLELES} alleles are possible; {ceil(BULK1_ALLELES * missingFilter)} needed to pass")
        print(f"# For bulk 2: {BULK2_ALLELES} alleles are possible; {ceil(BULK2_ALLELES * missingFilter)} needed to pass")
    
    # Parse the ED file
    edDict = {}
    starts, ends = [], []
    with read_gz_file(edFile) as fileIn:
        firstLine = True
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            if firstLine:
                if not sl == HEADER:
                    raise ValueError(f"ED file is improperly formatted; header line '{sl}' " + 
                                     f"does not match expected header '{HEADER}'")
                firstLine = False
            else:
                # Parse relevant details and validate format
                chrom, posi, variant, bulk1_alleles, bulk2_alleles, euclideanDistance = sl
                try:
                    posi = int(posi)
                except:
                    raise ValueError(f"Position '{posi}' is not an integer; offending line is '{line}'")
                try:
                    euclideanDistance = float(euclideanDistance)
                except:
                    raise ValueError(f"Euclidean distance '{euclideanDistance}' is not a float; offending line is '{line}'")
                try:
                    bulk1_alleles = int(bulk1_alleles)
                    bulk2_alleles = int(bulk2_alleles)
                except:
                    raise ValueError(f"Bulk allele counts '{bulk1_alleles}' or '{bulk2_alleles}' are not integers; " + 
                                     f"offending line is '{line}'")
                
                # Skip if missing data exceeds threshold
                if metadataDict != None:
                    if (bulk1_alleles / BULK1_ALLELES) < missingFilter or (bulk2_alleles / BULK2_ALLELES) < missingFilter:
                        continue
                
                # Store in dictionary
                edDict.setdefault(chrom, [[], []])
                edDict[chrom][0].append(posi)
                edDict[chrom][1].append(euclideanDistance)
    return edDict

def convert_dict_to_windowed_ncls(statDict, windowSize=0):
    '''
    Parameters:
        statDict -- a dictionary with structure like:
                  {
                      "chr1": [[pos1, pos2, ...], [stat1, stat2, ...]],
                      "chr2": [[pos1, pos2, ...], [stat1, stat2, ...]],
                      ...
                  }
        windowSize -- OPTIONAL; an integer indicating the size of the window that was used
                      when generating the depth file that led to the statistics file. Default is 0
                      (no window size) which is intended for use with variant calls, whereas
                      depth deletions should use an actual window size.
    Returns:
        windowedNCLS -- a WindowedNCLS object containing statistical values indexed by chromosome
                        and position
    '''
    windowedNCLS = WindowedNCLS(windowSize)
    for chrom, value in statDict.items():
        positions = np.array(value[0])
        statsValues = np.array(value[1])
        windowedNCLS.add(chrom, positions, statsValues)
    return windowedNCLS
