from .parsing import read_gz_file, vcf_header_to_metadata_validation, parse_vcf_genotypes

def filter_vcf(vcfFile, outputFileName, metadataDict, bulkMissing=0.25, qualThreshold=30.0):
    '''
    Parameters:
        vcfFile -- a string pointing to the VCF or VCF-like file to parse
        outputFileName -- a string indicating the output file name to write to.
        metadataDict -- a dictionary with structure like:
                        {
                            "bulk1": set([ "sample1", "sample2", ... ]),
                            "bulk2": set([ "sample3", "sample4", ... ])
                        }
        bulkMissing -- OPTIONAL; a float indicating the maximum percentage of
                       missing data in both bulk populations (default == 0.25);
                       it is allow for one bulk to have more missing data, we
                       only filter if BOTH bulks fail this threshold.
        qualThreshold -- OPTIONAL; a float indicating the minimum QUAL score
                         to keep a variant (default == 30.0).
    Yields:
        contig -- the contig name for the variant
        pos -- the position of the variant
        variant -- the type of variant (snp or indel)
        numAllelesB1 -- the number of genotyped alleles in bulk 1
        numAllelesB2 -- the number of genotyped alleles in bulk 2
        euclideanDist -- the Euclidean distance between the two bulks
    '''
    with read_gz_file(vcfFile) as fileIn, open(outputFileName, "w") as fileOut:
        for line in fileIn:
            sl = line.rstrip("\r\n").replace('"', '').split("\t") # remove quotations to help with files opened by Excel
            
            # Handle header line
            if line.startswith("#CHROM"):
                samples = sl[9:] # This gives us the ordered sample IDs
                vcf_header_to_metadata_validation(samples, metadataDict, strict=False)
                
                b1Samples = [ sample for sample in samples if sample in metadataDict["bulk1"] ]
                b1SampleNum = len(b1Samples)
                b2Samples = [ sample for sample in samples if sample in metadataDict["bulk2"] ]
                b2SampleNum = len(b2Samples)
            
            # Handle comment lines
            if line.startswith("#"):
                fileOut.write(line)
                continue
            
            # Extract relevant details from line
            try:
                qual = float(sl[5])
            except:
                qual = 0.0 # If the quality is missing, we'll just assume it's zero
            formatField = sl[8]
            sampleFields = sl[9:]
            
            # Omit variants with QUAL score below threshold
            if qual < qualThreshold:
                continue
            
            # Parse the genotypes out of the sample fields
            snpDict = parse_vcf_genotypes(formatField, sampleFields, samples)
            
            # Count the number of genotyped samples in each bulk
            numBulk1 = sum([ 1 for sample in b1Samples if sample in snpDict ])
            numBulk2 = sum([ 1 for sample in b2Samples if sample in snpDict ])
            
            # Calculate the percentage of missing data in each bulk
            bulk1Missing = 1.0 - (numBulk1 / b1SampleNum)
            bulk2Missing = 1.0 - (numBulk2 / b2SampleNum)
            
            # Omit variants with too much missing data
            if bulk1Missing > bulkMissing and bulk2Missing > bulkMissing:
                continue
            
            # Write the line to the output file
            fileOut.write(line)
