import numpy as np
from math import sqrt, ceil

from .parsing import read_gz_file, vcf_header_to_metadata_validation, parse_vcf_genotypes
from .ncls import WindowedNCLS

def calculate_snp_ed(b1Gt, b2Gt, isCNV=False):
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
    Returns:
        numAllelesB1 -- the number of genotyped alleles in bulk 1
        numAllelesB2 -- the number of genotyped alleles in bulk 2
        edist -- a float of the the Euclidean distance between the two bulks
    '''
    # Get all the unique alleles
    alleles = [ allele for gt in b1Gt + b2Gt for allele in gt ]
    uniqueAlleles = list(set(alleles))
    
    # Tally for bulk 1
    if isCNV:
        uniqueAlleles = [0, 1] # CNV genotypes get re-encoded to 0 and 1
        medianAllele = np.median(alleles)
        b1Count = { allele: 0 for allele in uniqueAlleles } # CNV genotypes are in relation to the median
        for allele1, allele2 in b1Gt:
            b1Count[0 if allele1 <= medianAllele else 1] += 1
            b1Count[0 if allele2 <= medianAllele else 1] += 1
    else:
        b1Count = { allele: 0 for allele in uniqueAlleles }
        for allele1, allele2 in b1Gt:
            b1Count[allele1] += 1
            b1Count[allele2] += 1
    
    # Tally for bulk 2
    if isCNV:
        b2Count = { allele: 0 for allele in uniqueAlleles } # CNV genotypes are in relation to the median
        for allele1, allele2 in b2Gt:
            b2Count[0 if allele1 <= medianAllele else 1] += 1 # medianAllele was calculated for bulk 1
            b2Count[0 if allele2 <= medianAllele else 1] += 1
    else:
        b2Count = { allele: 0 for allele in uniqueAlleles }
        for allele1, allele2 in b2Gt:
            b2Count[allele1] += 1
            b2Count[allele2] += 1
    
    # Sum the number of genotyped alleles for each bulk
    numAllelesB1 = sum(b1Count.values())
    numAllelesB2 = sum(b2Count.values())
    
    # Calculate the Euclidean distance between the two bulks if possible
    if numAllelesB1 == 0 or numAllelesB2 == 0:
        return numAllelesB1, numAllelesB2, 0 # euclidean distance cannot be calculated
    else:
        # Derive our euclidean distance value
        """Refer to "Euclidean distance calculation" in Hill et al. 2013"""
        edist = sqrt(sum([
            ((b1Count[allele] / numAllelesB1) - (b2Count[allele] / numAllelesB2))**2
            for allele in uniqueAlleles
        ]))
        
        # Return the values
        return numAllelesB1, numAllelesB2, edist

def parse_vcf_for_ed(vcfFile, metadataDict, isCNV, ignoreIdentical=True, quiet=False):
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
            bulk1 = [ snpDict[sample] for sample in metadataDict["bulk1"] if sample in snpDict ]
            bulk2 = [ snpDict[sample] for sample in metadataDict["bulk2"] if sample in snpDict ]
            
            # Calculate difference ratio
            numAllelesB1, numAllelesB2, \
                euclideanDist = calculate_snp_ed(bulk1, bulk2, isCNV=isCNV)
            
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
