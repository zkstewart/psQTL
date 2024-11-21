from math import sqrt

from .parsing import read_gz_file, vcf_header_to_metadata_validation, parse_vcf_genotypes

def calculate_snp_ed(b1Gt, b2Gt):
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
        numAllelesB1 -- the number of genotyped alleles in bulk 1
        numAllelesB2 -- the number of genotyped alleles in bulk 2
        edist -- a float of the the Euclidean distance between the two bulks
    '''
    # Get all the unique alleles
    alleles = list(set([ allele for gt in b1Gt + b2Gt for allele in gt ]))
    if 0 not in alleles:
        alleles.append(0)
    alleles.sort()
    
    # Tally for bulk 1
    b1Count = { allele: 0 for allele in alleles }
    for allele1, allele2 in b1Gt:
        b1Count[allele1] += 1
        b1Count[allele2] += 1
    
    # Tally for bulk 2
    b2Count = { allele: 0 for allele in alleles }
    for allele1, allele2 in b2Gt:
        b2Count[allele1] += 1
        b2Count[allele2] += 1
    
    # Sum the number of genotyped alleles for each bulk
    numAllelesB1 = sum(b1Count.values())
    numAllelesB2 = sum(b2Count.values())
    
    # Calculate the Euclidean distance between the two bulks if possible
    if numAllelesB1 == 0 and numAllelesB2 == 0:
        return numAllelesB1, numAllelesB2, "." # euclidean distance is null
    elif numAllelesB1 == 0 or numAllelesB2 == 0:
        return numAllelesB1, numAllelesB2, 1 # euclidean distance is 1
    else:
        # Derive our euclidean distance value
        """Refer to "Euclidean distance calculation" in Hill et al. 2013"""
        edist = sqrt(sum([
            ((b1Count[allele] / numAllelesB1) - (b2Count[allele] / numAllelesB2))**2
            for allele in alleles
        ]))
        
        # Return the values
        return numAllelesB1, numAllelesB2, edist

def parse_vcf_for_ed(vcfFile, metadataDict, ignoreIdentical=False):
    '''
    Parameters:
        vcfFile -- a string pointing to the VCF or VCF-like file to parse
        metadataDict -- a dictionary with structure like:
                        {
                            "bulk1": set([ "sample1", "sample2", ... ]),
                            "bulk2": set([ "sample3", "sample4", ... ])
                        }
        ignoreIdentical -- OPTIONAL; a boolean indicating whether to ignore
                           identical non-reference alleles shared by all samples
    Yields:
        contig -- the contig name for the variant
        pos -- the position of the variant
        variant -- the type of variant (snp or indel)
        numAllelesB1 -- the number of genotyped alleles in bulk 1
        numAllelesB2 -- the number of genotyped alleles in bulk 2
        euclideanDist -- the Euclidean distance between the two bulks
    '''
    with read_gz_file(vcfFile) as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").replace('"', '').split("\t") # remove quotations to help with files opened by Excel
            
            # Handle header line
            if line.startswith("#CHROM"):
                samples = sl[9:] # This gives us the ordered sample IDs
                vcf_header_to_metadata_validation(samples, metadataDict, strict=False)
            if line.startswith("#"):
                continue
            
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
            elif any([ len(ref_alt[0]) != len(ref_alt[x]) for x in range(1, len(ref_alt))]):
                variant = "indel"
            else:
                variant = "snp"
            
            # Split sample genotypes into bulk1 and bulk2
            bulk1 = [ snpDict[sample] for sample in metadataDict["bulk1"] if sample in snpDict ]
            bulk2 = [ snpDict[sample] for sample in metadataDict["bulk2"] if sample in snpDict ]
            
            # Calculate difference ratio
            numAllelesB1, numAllelesB2, \
                euclideanDist = calculate_snp_ed(bulk1, bulk2)
            
            # Skip if both bulks are identical
            if ignoreIdentical and euclideanDist == 0:
                bulk1Dedup = set(( tuple(x) for x in bulk1 ))
                bulk2Dedup = set(( tuple(x) for x in bulk2 ))
                "if both have set len==1, are the same, and have 0 Euclidean distance, they are identical non-reference alleles"
                if len(bulk1Dedup) == 1 and bulk1Dedup == bulk2Dedup:
                    continue
            
            # Yield results
            yield contig, pos, variant, numAllelesB1, numAllelesB2, euclideanDist
