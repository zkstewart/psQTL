import numpy as np
from math import sqrt, ceil
from statistics import mean
from collections import Counter
from itertools import combinations

from .parsing import read_gz_file, vcf_header_to_metadata_validation, parse_vcf_genotypes
from .ncls import WindowedNCLS

def gt_median_adjustment(genotypeLists):
    '''
    Receives a list of genotype lists, obtains the median of all alleles found in the lists,
    and adjusts the genotypes to be relative to the median. This is specifically intended for
    use on CNV genotypes, which are often highly variable and range from 0 up to very high values.
    
    Adjusted values range from 0 to 2, where:
        0 = allele is below the median
        1 = allele is equal to the median
        2 = allele is above the median
    This adjustment allows for easier comparison of CNV genotypes across samples while
    mitigating the effects of variability in CNV depth.
    
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
                    adjustedAllele = 0 if allele < medianAllele else 1 if allele == medianAllele else 2
                    gtList.append(adjustedAllele)
            adjustedSublist.append(gtList)
        adjustedGenotypes.append(adjustedSublist)
    return adjustedGenotypes

def possible_genotypes(gt1, gt2):
    '''
    Calculate all possible inheritable genotypes from two parent genotypes.
    
    Parameters:
        gt1 / gt2 -- a list of integers representing the genotype values
                     with format like:
                     [0, 1, 2, ...]
    Returns:
        possibleGTs -- a set of frozensets containing all possible genotypes
                       that can be formed from the two genotypes; sets
                       mask duplicated alleles but allow for testing
                       of genotype without regard to order
    '''
    numGt1 = int(len(gt1) / 2)
    numGt2 = int(len(gt2) / 2)
    
    possibleGTs = []
    for g1 in combinations(gt1, numGt1):
        for g2 in combinations(gt2, numGt2):
            possibleGTs.append(frozenset(g1 + g2))
    
    return set(possibleGTs)

def calculate_allele_frequency_ed(g1Gt, g2Gt):
    '''
    Calculates the Euclidean distance between two groups based on the allele
    frequencies in each group. This is conceptually similar to BSA methods like
    QTLseq, but substituting allele depth (AD) for the number of alleles
    in genotype (GT) calls.
    
    Parameters:
        g1Gt / g2Gt -- a list of lists containing the genotype value as integers
                       with format like:
                       [
                           [0, 1],
                           [0, 0],
                           [1, 1],
                           ...
                       ]
    Returns:
        numAllelesG1 -- the number of genotyped alleles in group 1
        numAllelesG2 -- the number of genotyped alleles in group 2
        edist -- a float of the the Euclidean distance between the two groups
    '''
    # Get all the unique alleles
    alleles = [ allele for gt in g1Gt + g2Gt for allele in gt ]
    uniqueAlleles = list(set(alleles))
    
    # Tally for group 1
    b1Count = { allele: 0 for allele in uniqueAlleles }
    for gt in g1Gt:
        for allele in gt:
            b1Count[allele] += 1
    
    # Tally for group 2
    b2Count = { allele: 0 for allele in uniqueAlleles }
    for gt in g2Gt:
        for allele in gt:
            b2Count[allele] += 1
    
    # Sum the number of genotyped alleles for each group
    numAllelesG1 = sum(b1Count.values())
    numAllelesG2 = sum(b2Count.values())
    
    # Calculate the Euclidean distance between the two groups if possible
    if numAllelesG1 == 0 or numAllelesG2 == 0:
        return numAllelesG1, numAllelesG2, 0 # euclidean distance cannot be calculated
    else:
        # Derive our euclidean distance value
        "Refer to 'Euclidean distance calculation' in Hill et al. 2013"
        edist = sqrt(sum([
            ((b1Count[allele] / numAllelesG1) - (b2Count[allele] / numAllelesG2))**2
            for allele in uniqueAlleles
        ]))
        # Return the values
        return numAllelesG1, numAllelesG2, edist

def calculate_genotype_frequency_ed(g1Gt, g2Gt):
    '''
    Calculates the Euclidean distance between two groups based on the
    genotype frequencies in each group. Rather than breaking apart genotypes (GT)
    into individual alleles as calculate_allele_frequency_ed() does, this method
    treats each genotype as a single entity and compares the frequency of each
    genotype in each group.
    
    However, although this sounds like a good idea, it has some problems. From
    real data, this method tends to overestimate segregation between groups. The
    cause of this is unclear. My hunch is that reads, possibly originating from
    repetitive regions in proximity to real QTLs, can contaminate those related
    repeats throughout the rest of the genome, leading to signals of segregation
    occurring throughout the whole genome rather than just at the QTL. Whatever
    it is, this method is not recommended for use in real data analysis despite
    it sounding fantastic in theory.
    
    Parameters:
        g1Gt / g2Gt -- a list of lists containing the genotype value as integers
                       with format like:
                       [
                           [0, 1],
                           [0, 0],
                           [1, 1],
                           ...
                       ]
    Returns:
        numSamplesG1 -- the number of genotyped samples in group 1
        numSamplesG2 -- the number of genotyped samples in group 2
        edist -- a float of the the Euclidean distance between the two groups
    '''
    # Get all the unique genotypes
    uniqueGts = set(tuple(sorted(gt)) for gt in g1Gt + g2Gt)
    
    # Tally for group 1
    g1GtCount = { gt: 0 for gt in uniqueGts }
    for gt in g1Gt:
        g1GtCount[tuple(sorted(gt))] += 1
    
    # Tally for group 2
    g2GtCount = { gt: 0 for gt in uniqueGts }
    for gt in g2Gt:
        g2GtCount[tuple(sorted(gt))] += 1
    
    # Sum the number of genotyped samples for each group
    numSamplesG1 = len(g1Gt)
    numSamplesG2 = len(g2Gt)
    
    # Calculate the Euclidean distance between the two groups if possible
    if numSamplesG1 == 0 or numSamplesG2 == 0:
        return numSamplesG1, numSamplesG2, 0 # euclidean distance cannot be calculated
    else:
        # Derive our euclidean distance value
        "This is the same as in calculate_allele_frequency_ed(), but for genotypes instead of alleles"
        edist = sqrt(sum([
            ((g1GtCount[gt] / numSamplesG1) - (g2GtCount[gt] / numSamplesG2))**2
            for gt in uniqueGts
        ]))
        
        # Return the values
        return numSamplesG1, numSamplesG2, edist

def calculate_inheritance_ed(g1Gt, g2Gt, parentsGT):
    '''
    Employs a method to calculate the Euclidean distance between two groups
    based on the likelihood of specific alleles/haplotypes being inherited
    from each parent. This blends the allele frequency and genotype frequency
    methods to provide a more accurate representation of _biased inheritance_
    of specific alleles from parents to progeny, rather than a simple
    comparison of genotype frequency.
    
    Parameters:
        g1Gt / g2Gt -- a list of lists containing the genotype value as integers
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
        numAllelesG1 -- the number of genotyped alleles in group 1
        numAllelesG2 -- the number of genotyped alleles in group 2
        edist -- a float of the the Euclidean distance between the two groups
    '''
    # Validate the parentsGT input
    if len(parentsGT) != 2:
        raise ValueError("parentsGT must be a list of two lists containing the genotype values for the parents")
    PARENT_FULLY_ASSIGNED = { f"p{i+1}" : int(len(gt) / 2) for i, gt in enumerate(parentsGT) }
    NUM_PARENT_ALLELES = sum(PARENT_FULLY_ASSIGNED.values())
    
    # Raise error for unmanageable ploidy values
    for i, gt in enumerate(parentsGT):
        if len(gt) % 2 != 0:
            raise ValueError(f"Odd number of alleles found in parent {i+1}'s genotype; " +
                             "cannot calculate inheritance Euclidean distance with odd ploidy")
    for groupGt in [g1Gt, g2Gt]:
        for gt in groupGt:
            if len(gt) != NUM_PARENT_ALLELES:
                raise ValueError("Progeny samples must have a number of alleles equal to the summed value of " +
                                 f"half of each parent's chromosomes ({NUM_PARENT_ALLELES}); " +
                                 "cannot calculate inheritance Euclidean distance with mismatched ploidy")
    
    # Assign alleles to parent haplotypes
    groupSums = []
    for i, groupGt in enumerate([g1Gt, g2Gt]): # for each group
        # Set up data structure to hold assigned allele counts across the group
        groupSum = {
            "p1": {
                gt: 0 for gt in parentsGT[0]
            },
            "p2": {
                gt: 0 for gt in parentsGT[1]
            }
        }
        # For each sample, assign alleles to parents based on the most likely inheritance
        for gt in groupGt:
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
                    groupSum[parent][allele] += assigned
        
        # Set aside results for this group
        groupSums.append(groupSum)
    
    # Format groupSums to ensure all alleles are present [makes it easier to calculate ED without nested if statements]
    uniqueAlleles = list(set([ allele for gt in parentsGT for allele in gt ]))
    for allele in uniqueAlleles:
        for bSum in groupSums:
            if allele not in bSum["p1"]:
                bSum["p1"][allele] = 0
            if allele not in bSum["p2"]:
                bSum["p2"][allele] = 0
    b1Sum, b2Sum = groupSums
    
    # Get the number of samples and alleles in each group
    numSamplesG1, numSamplesG2 = len(g1Gt), len(g2Gt)
    numAllelesG1 = sum([ len(gt) for gt in g1Gt ])
    numAllelesG2 = sum([ len(gt) for gt in g2Gt ])
    
    # Calculate the Euclidean distance between the parental inherited alleles if possible
    """Note that values are divied by the number of samples, not the number of alleles;
    doing so gives results in the same scale as the standard segregant ED calculation,
    but the rationale is difficult for me to articulate"""
    if numSamplesG1 == 0 or numSamplesG2 == 0:
        return numAllelesG1, numAllelesG2, 0 # euclidean distance cannot be calculated
    else:
        edist = sqrt(sum([
            ((b1Sum[parent][allele] / numSamplesG1) - (b2Sum[parent][allele] / numSamplesG2))**2
            for parent in ["p1", "p2"]
            for allele in uniqueAlleles
        ])/2) # divide by 2 since there are two parents and 2x the number of comparions made across groups
    
    # Return the values
    return numAllelesG1, numAllelesG2, edist

def filter_impossible_genotypes(g1Gt, g2Gt, parentsGT):
    '''
    Parameters:
        g1Gt / g2Gt -- a list of lists containing the genotype value as integers
                       with format like:
                          [
                            [0, 1],
                            [0, 0],
                            [1, 1],
                            ...
                          ]
        parentsGT -- a list of two lists containing the genotype value as integers
                     for the parents with format like:
                        [ [0, 1], [1, 2] ]
    Returns:
        g1Gt / g2Gt -- a modified version of the original g1 lists with
                       impossible genotypes filtered out
    '''
    possibleGTs = possible_genotypes(parentsGT[0], parentsGT[1])
    g1Gt = [ gt for gt in g1Gt if set(gt) in possibleGTs ]
    g2Gt = [ gt for gt in g2Gt if set(gt) in possibleGTs ]
    return g1Gt, g2Gt

def parse_vcf_for_ed(vcfFile, metadataDict, isCNV, parents=[],
                     ignoreIdentical=True, quiet=False):
    '''
    Parameters:
        vcfFile -- a string pointing to the VCF or VCF-like file to parse
        metadataDict -- a dictionary with structure like:
                        {
                            "group1": set([ "sample1", "sample2", ... ]),
                            "group2": set([ "sample3", "sample4", ... ])
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
        resultsDict -- a dictionary with the following keys:
            contig -- the contig name for the variant
            pos -- the position of the variant
            variant -- the type of variant (snp or indel)
            numAllelesG1 -- the number of genotyped alleles in group 1
            numAllelesG2 -- the number of genotyped alleles in group 2
            numFilteredG1 -- the number of genotyped alleles in group 1 after filtering;
                             None if parents are not provided or specified
            numFilteredG2 -- the number of genotyped alleles in group 2 after filtering;
                             None if parents are not provided or specified
            alleleED -- a float of the the Euclidean distance between the two groups
                        when using the naive allele frequency
            genotypeED -- a float of the the Euclidean distance between the two groups
                          when using the naive genotype frequency calculation
            inheritanceED -- a float of the the Euclidean distance between the two groups
                             when using the inheritance calculation;
                             None if parents are not provided or specified
            ploidy -- the ploidy of the samples in the VCF file at this position
            possibleAllelesG1 -- the number of alleles if all group 1 samples were present
            possibleAllelesG2 -- the number of alleles if all group 2 samples were present
    '''
    # Validations
    if parents == None: # just in case
        parents = []
    
    if parents != [] and len(parents) != 2:
        raise ValueError("Parents must be a list with zero (standard ED) or two sample IDs (haplotype inheritance ED)")
        for parent in parents:
            if parent not in metadataDict["group1"] and parent not in metadataDict["group2"]:
                raise ValueError(f"Parent sample '{parent}' is not in either group; cannot calculate haplotype inheritance ED")
    
    # Iterate through the VCF file
    samples = None
    with read_gz_file(vcfFile) as fileIn:
        for line in fileIn:
            l = line.strip('\r\n\t "') # remove quotations to help with files opened by Excel
            sl = l.replace('"', '').split("\t") # remove any remaining quotations
            
            # Skip blank lines
            if l == "":
                continue
            
            # Handle header lines
            if l.startswith("#CHROM"):
                # Get the ordered sample IDs from the header line
                samples = sl[9:]
                
                # Validate the metadata dictionary against the samples
                vcf_header_to_metadata_validation(samples, metadataDict, strict=False, quiet=quiet)
                
                # Check that the parents are in the samples
                if parents != []:
                    for parent in parents:
                        if parent not in samples:
                            raise ValueError(f"Parent sample '{parent}' is not in the VCF file; check your --parents argument")
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
            snpDict = parse_vcf_genotypes(formatField, sampleFields, samples) # extracts ALL sample genotypes, not just those in metadataDict
            
            # Figure out what type of variant this is
            if any([ x == "." for x in ref_alt ]):
                variant = "indel"
            elif all([ x == "N" for x in ref_alt ]): # for parsing depth VCF-like files
                variant = "cnv"
            elif any([ len(ref_alt[0]) != len(ref_alt[x]) for x in range(1, len(ref_alt))]):
                variant = "indel"
            else:
                variant = "snp"
            
            # Split sample genotypes into group1 and group2
            g1Gt, g2Gt, parentsGT = separate_genotypes_by_group(snpDict, metadataDict, parents=[]) # don't use parents yet
            
            # Skip if no genotypes are found in either group
            "This is a completely useless variant, so we skip it"
            if len(g1Gt) == 0 and len(g2Gt) == 0: # since we pass here, inferred ploidy is never == None
                continue
            
            # Adjust values if this is a CNV
            if isCNV:
                g1Gt, g2Gt = gt_median_adjustment([g1Gt, g2Gt])
            
            # Calculate naive Euclidean distances
            "These two ED measures do not utilise parent genotypes and are hence considered to be naive"
            numAllelesG1, numAllelesG2, alleleED = calculate_allele_frequency_ed(g1Gt, g2Gt)
            _, _, genotypeED = calculate_genotype_frequency_ed(g1Gt, g2Gt) # don't need the allele count again
            
            # Skip if both groups are identical
            if ignoreIdentical and alleleED == 0:
                g1Dedup = set(( tuple(x) for x in g1Gt ))
                g2Dedup = set(( tuple(x) for x in g2Gt ))
                "if both have set len==1, are the same, and have 0 Euclidean distance, they are identical non-reference alleles"
                if len(g1Dedup) == 1 and g1Dedup == g2Dedup:
                    continue
            
            # Calculate inheritance Euclidean distance if parents are provided
            numFilteredG1, numFilteredG2, inheritanceED = 0, 0, 0 # these will be set later if inheritance ED is calculated
            if parents != [] and len(parents) == 2:
                # Re-obtain our group genotypes, this time using the parents
                g1Gt, g2Gt, parentsGT = separate_genotypes_by_group(snpDict, metadataDict, parents=parents)
                
                # Skip if one or both parents are not genotyped at this position
                if len(parentsGT) != 2: # parentsGT can come back missing parents if they are not called
                    pass
                else:
                    # Adjust values if this is a CNV
                    if isCNV:
                        g1Gt, g2Gt = gt_median_adjustment([g1Gt, g2Gt])
                    
                    # Filter impossible progeny genotypes based on the parents' genotypes
                    g1Gt, g2Gt = filter_impossible_genotypes(g1Gt, g2Gt, parentsGT)
                    
                    # Calculate inheritance Euclidean distance
                    numFilteredG1, numFilteredG2, inheritanceED = calculate_inheritance_ed(g1Gt, g2Gt, parentsGT)
            
            # Infer the ploidy of the samples
            ploidy = infer_ploidy(snpDict) # this may average across samples with different ploidies, but that's okay for this purpose
            possibleAllelesG1 = int(len(metadataDict["group1"]) * ploidy)
            possibleAllelesG2 = int(len(metadataDict["group2"]) * ploidy)
            
            # Yield results
            yield {
                "contig": contig, "pos": pos, "variant": variant,
                "numAllelesG1": numAllelesG1, "numAllelesG2": numAllelesG2,
                "numFilteredG1": numFilteredG1, "numFilteredG2": numFilteredG2,
                "alleleED": alleleED, "genotypeED": genotypeED, "inheritanceED": inheritanceED,
                "ploidy": ploidy, "possibleAllelesG1": possibleAllelesG1, "possibleAllelesG2": possibleAllelesG2
            }

def separate_genotypes_by_group(snpDict, metadataDict, parents=[]):
    '''
    Parameters:
        snpDict -- a dictionary with structure like:
                   {
                       "sample1": [0, 1],
                       "sample2": [0, 0],
                       ...
                   }
        metadataDict -- a dictionary with structure like:
                        {
                            "group1": set([ "sample1", "sample2", ... ]),
                            "group2": set([ "sample3", "sample4", ... ])
                        }
        parents -- (OPTIONAL) a list of two sample IDs indicating parent samples
                   OR an empty list if no parents are available or specified;
                   default is an empty list
    Returns:
        g1Gt -- a list of lists containing the genotype values for group 1
        g2Gt -- a list of lists containing the genotype values for group 2
        parentsGT -- a list of lists containing the genotype values for the parents
    '''
    g1Gt = [ snpDict[sample] for sample in metadataDict["group1"] if sample in snpDict and not sample in parents ]
    g2Gt = [ snpDict[sample] for sample in metadataDict["group2"] if sample in snpDict and not sample in parents ]
    parentsGT = [ snpDict[parent] for parent in parents if parent in snpDict ] # if parents == [] this will always be empty
    return g1Gt, g2Gt, parentsGT

def infer_ploidy(snpDict):
    '''
    Parameters:
        snpDict -- a dictionary with structure like:
                   {
                       "sample1": [0, 1],
                       "sample2": [0, 0],
                       ...
                   }
    Returns:
        inferredPloidy -- a float indicating the inferred ploidy of the samples
    '''
    try:
        return mean([ len(gt) for gt in snpDict.values() ])
    except:
        return None

def parse_ed_as_dict(edFile, missingFilter=0.5):
    '''
    Parameters:
        edFile -- a string indicating the path to an ED file
        missingFilter -- (OPTIONAL) a float indicating the maximum percentage of
                         samples that exist in a group which can be missing
                         without filtration occurring; 0 is perfectly strict
                         and 1 is perfectly lenient; default is 1 (100% missing data allowed)
    Returns:
        edDict -- a dictionary with structure like:
                  {
                      "chr1": [[pos1, pos2, ...], [ed1, ed2, ...]],
                      "chr2": [[pos1, pos2, ...], [ed1, ed2, ...]],
                      ...
                  }
    '''
    HEADER = ["CHROM", "POSI", "variant", "group1_alleles", "group1_alleles_possible",
              "group2_alleles", "group2_alleles_possible", "euclideanDist"]
    
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
                chrom, posi, variant, group1_alleles, group1_alleles_possible, \
                    group2_alleles, group2_alleles_possible, euclideanDistance = sl
                try:
                    posi = int(posi)
                except:
                    raise ValueError(f"Position '{posi}' is not an integer; offending line is '{line}'")
                try:
                    euclideanDistance = float(euclideanDistance)
                except:
                    raise ValueError(f"Euclidean distance '{euclideanDistance}' is not a float; offending line is '{line}'")
                try:
                    group1_alleles = int(group1_alleles)
                    group2_alleles = int(group2_alleles)
                except:
                    raise ValueError(f"Group allele counts '{group1_alleles}' or '{group2_alleles}' are not integers; " + 
                                     f"offending line is '{line}'")
                try:
                    group1_alleles_possible = int(group1_alleles_possible)
                    group2_alleles_possible = int(group2_alleles_possible)
                except:
                    raise ValueError(f"Possible group allele counts '{group1_alleles_possible}' or '{group2_alleles_possible}' " + 
                                     f"are not integers; offending line is '{line}'")
                
                # Calculate the percentage of missing data in each group
                missingG1 = (group1_alleles_possible - group1_alleles) / group1_alleles_possible if group1_alleles_possible > 0 else 0
                missingG2 = (group2_alleles_possible - group2_alleles) / group2_alleles_possible if group2_alleles_possible > 0 else 0
                
                # Skip if either group has too much missing data
                if missingG1 > missingFilter or missingG2 > missingFilter:
                    continue
                
                # Store in dictionary
                edDict.setdefault(chrom, [[], []])
                edDict[chrom][0].append(posi)
                edDict[chrom][1].append(euclideanDistance)
    return edDict

def convert_dict_to_windowed_ncls(statDict, windowSize=1):
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
                      depth CNVs should use an actual window size.
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
