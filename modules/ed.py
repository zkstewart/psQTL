from math import sqrt

from .parsing import SimpleGenotypeIterator

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
    # Extract metadata sample IDs from dict
    metadataSamples = set([
        sample
        for samples in metadataDict.values()
        for sample in samples
    ])
    
    # Iterate through VCF / VCF-like file
    firstYield = True
    for values in SimpleGenotypeIterator(vcfFile):
        # Grab header line containing sample IDs
        if firstYield:
            vcfSamples = set(values) # first yield gives the sample IDs
            firstYield = False
            
            # Check that the metadata file matches the VCF file
            if vcfSamples != metadataSamples:
                vcfDiff = vcfSamples.difference(metadataSamples)
                metadataDiff = metadataSamples.difference(vcfSamples)
                
                # Error out if files are incompatible
                if len(metadataDiff) == len(metadataSamples):
                    raise ValueError("Metadata file has no samples in common with VCF file; " +
                                     f"Metadata samples: {metadataSamples}\nVCF samples: {vcfSamples}")
                
                # Error out if we don't have samples from both bulks
                b1Samples = [ sample for sample in vcfSamples if sample in metadataDict["bulk1"] ]
                if len(b1Samples) == 0:
                    raise ValueError("No samples from bulk 1 are present in the VCF file.")
                
                b2Samples = [ sample for sample in vcfSamples if sample in metadataDict["bulk2"] ] 
                if len(b2Samples) == 0:
                    raise ValueError("No samples from bulk 2 are present in the VCF file.")
                
                # Warn if some samples are missing
                needsSpace = False
                if len(vcfDiff) > 0:
                    print("# WARNING: In your VCF, the following samples exist which are " + 
                        "absent from the metadata: ", ", ".join(vcfDiff))
                    print("# These will be ignored during the Euclidean distance calculation")
                    needsSpace = True
                if len(metadataDiff) > 0:
                    print("# WARNING: In your metadata, the following samples exist which are " + 
                        "absent from the VCF: ", ", ".join(metadataDiff))
                    print("# These will be ignored during the Euclidean distance calculation")
                    needsSpace = True
                if needsSpace:
                    print("")
                
                # Notify user of samples that will be used
                print(f"# Samples used as part of bulk 1 (n={len(b1Samples)}) include: " + ", ".join(b1Samples))
                print(f"# Samples used as part of bulk 2 (n={len(b2Samples)}) include: " + ", ".join(b2Samples))
        # Handle content lines
        else:
            contig, pos, ref, alt, snpDict = values
            ref_alt = [ref, *alt]
            
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
